use std::cmp::min;
use std::collections::{HashSet, VecDeque};
use std::sync::{Arc, mpsc};
use std::thread;
use roaring::RoaringBitmap;

pub fn connected_components(bit_vectors: Vec<RoaringBitmap>, minimum_similarity: f64) -> Vec<Vec<usize>> {
    let graph = create_graph(bit_vectors, minimum_similarity);
    let components = bfs(graph);
    components
}

fn create_graph(bit_vectors: Vec<RoaringBitmap>, minimum_similarity: f64) -> Vec<Vec<usize>> {
    let mut graph = vec![vec![]; bit_vectors.len()];
    let bit_vectors_arc = Arc::new(bit_vectors);
    let (sender, receiver) = mpsc::sync_channel(64);
    for i1 in 0..bit_vectors_arc.len() {
        let sender_clone = sender.clone();
        let kmer_sets_arc_clone = bit_vectors_arc.clone();
        thread::spawn(move || {
            let mut edges = vec![];
            for i2 in 0..kmer_sets_arc_clone.len() {
                if i2 <= i1 {
                    continue;
                }
                let (bit_vector_1, bit_vector_2) = (&kmer_sets_arc_clone[i1], &kmer_sets_arc_clone[i2]);
                let intersect_size = intersect(bit_vector_1, bit_vector_2);
                let minimum_containment = intersect_size as f64 / min(bit_vector_1.len(), bit_vector_2.len()) as f64;
                if minimum_containment >= minimum_similarity {
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

fn intersect(vector_1: &RoaringBitmap, vector_2: &RoaringBitmap) -> usize {
    vector_1.intersection_len(vector_2) as usize
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
