use std::cmp::min;
use std::collections::{HashSet, VecDeque};
use std::sync::{Arc, mpsc};
use std::thread;

pub fn connected_components(sorted_vectors: Vec<Vec<usize>>, minimum_similarity: f64) -> Vec<Vec<usize>> {
    let graph = create_graph(sorted_vectors, minimum_similarity);
    let components = bfs(graph);
    components
}

/// Returns the connected components of all nodes
fn bfs(graph: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    let mut explored = HashSet::new();
    let mut connected_components = Vec::new();
    for s in 0..graph.len() {
        if explored.contains(&s) {
            continue;
        }
        connected_components.push(bfs_helper(&graph, s, &mut explored));
    }
    connected_components
}

fn bfs_helper(graph: &Vec<Vec<usize>>, start_node: usize, explored: &mut HashSet<usize>) -> Vec<usize> {
    explored.insert(start_node);
    let mut queue = VecDeque::from([start_node]);
    let mut connected_component = Vec::from([start_node]);
    while !queue.is_empty() {
        let node = queue.pop_front().unwrap();
        for adjacent_node in &graph[node] {
            if explored.contains(adjacent_node) {
                continue;
            }
            queue.push_back(*adjacent_node);
            explored.insert(*adjacent_node);
            connected_component.push(*adjacent_node);
        }
    }
    connected_component
}

fn create_graph(kmer_sets: Vec<Vec<usize>>, minimum_similarity: f64) -> Vec<Vec<usize>> {
    let mut graph = vec![vec![]; kmer_sets.len()];
    let kmer_sets_arc = Arc::new(kmer_sets);
    let (sender, receiver) = mpsc::sync_channel(64);
    for i1 in 0..kmer_sets_arc.len() {
        let sender_clone = sender.clone();
        let kmer_sets_arc_clone = kmer_sets_arc.clone();
        thread::spawn(move || {
            let mut edges = vec![];
            for i2 in 0..kmer_sets_arc_clone.len() {
                if i2 <= i1 {
                    continue;
                }
                let (sorted_vector_1, sorted_vector_2) = (&kmer_sets_arc_clone[i1], &kmer_sets_arc_clone[i2]);
                let intersect_size = intersect(sorted_vector_1, sorted_vector_2).len();
                let min_coverage = intersect_size as f64 / min(sorted_vector_1.len(), sorted_vector_2.len()) as f64;
                if min_coverage >= minimum_similarity {
                    edges.push((i1, i2));
                }
            }
            sender_clone.send(edges).unwrap();
        });
    }
    drop(sender);
    for edges in receiver {
        for edge in edges {
            graph[edge.0].push(edge.1);
            graph[edge.1].push(edge.0);
        }
    }
    graph
}

fn intersect(vector_1: &Vec<usize>, vector_2: &Vec<usize>) -> Vec<usize> {
    let mut index_1 = 0;
    let mut index_2 = 0;
    let mut intersection = vec![];
    while index_1 < vector_1.len() && index_2 < vector_2.len() {
        let first_value = vector_1[index_1];
        let second_value = vector_2[index_2];
        if first_value < second_value {
            match vector_1[index_1..vector_1.len()].binary_search(&second_value) {
                Ok(index) => {index_1 += index;}
                Err(index) => {index_1 += index;}
            }
        } else if first_value > second_value {
            match vector_2[index_2..vector_2.len()].binary_search(&first_value) {
                Ok(index) => {index_2 += index;}
                Err(index) => {index_2 += index;}
            }
        } else {
            // The values are equivalent
            intersection.push(first_value);
            index_1 += 1;
        }
    }
    intersection
}
