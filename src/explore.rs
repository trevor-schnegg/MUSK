use std::cmp::min;
use std::collections::{HashSet, VecDeque};
use std::sync::{Arc, mpsc};
use std::thread;
use vers_vecs::RsVec;

pub fn connected_components(bit_vectors: Vec<RsVec>, minimum_similarity: f64) -> Vec<Vec<usize>> {
    let graph = create_graph(bit_vectors, minimum_similarity);
    let components = bfs(graph);
    components
}

fn create_graph(bit_vectors: Vec<RsVec>, minimum_similarity: f64) -> Vec<Vec<usize>> {
    let mut graph = vec![vec![]; bit_vectors.len()];
    let kmer_sets_arc = Arc::new(bit_vectors);
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
                let intersect_size = intersect(sorted_vector_1, sorted_vector_2);
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

fn intersect(vector_1: &RsVec, vector_2: &RsVec) -> usize {
    let (mut iterator_1, mut iterator_2) = (vector_1.iter1(), vector_2.iter1());
    let (mut current_1, mut current_2) = (iterator_1.next(), iterator_2.next());
    let mut intersection_size = 0_usize;
    loop {
        match (current_1, current_2) {
            (Some(x), Some(y)) => {
                if x < y {
                    current_1 = iterator_1.next();
                } else if y > x {
                    current_2 = iterator_2.next();
                } else {
                    intersection_size += 1;
                    current_1 = iterator_1.next();
                }
            },
            _ => break,
        }
    }
    intersection_size
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
