use std::collections::{HashSet, VecDeque};
use std::sync::{mpsc, Arc};
use threadpool::ThreadPool;

use crate::sorted_vector_utilities::{IntersectIterator, UnionIterator};

pub fn connected_components(
    sorted_vector_sets: Vec<Vec<u32>>,
    minimum_similarity: f64,
    thread_number: usize,
) -> Vec<Vec<usize>> {
    let graph = create_graph(sorted_vector_sets, minimum_similarity, thread_number);
    let components = bfs(graph);
    components
}

fn create_graph(
    sorted_vector_sets: Vec<Vec<u32>>,
    minimum_similarity: f64,
    thread_number: usize,
) -> Vec<Vec<usize>> {
    let mut graph = vec![vec![]; sorted_vector_sets.len()];
    let sorted_vector_sets_arc = Arc::new(sorted_vector_sets);
    let (sender, receiver) = mpsc::channel();
    let pool = ThreadPool::new(thread_number);
    for index_1 in 0..sorted_vector_sets_arc.len() {
        let sender_clone = sender.clone();
        let sorted_vector_sets_arc_clone = sorted_vector_sets_arc.clone();
        pool.execute(move || {
            for index_2 in 0..sorted_vector_sets_arc_clone.len() {
                if index_2 <= index_1 {
                    continue;
                }
                let (bit_vector_1, bit_vector_2) = (
                    &sorted_vector_sets_arc_clone[index_1],
                    &sorted_vector_sets_arc_clone[index_2],
                );
                let intersect_size = IntersectIterator::from(bit_vector_1, bit_vector_2).count();
                let union_size = UnionIterator::from(vec![bit_vector_1, bit_vector_2]).count();
                let intersection_coverage = intersect_size as f64 / union_size as f64;
                if intersection_coverage >= minimum_similarity {
                    sender_clone.send((index_1, index_2)).unwrap();
                }
            }
        });
    }
    drop(sender);
    drop(sorted_vector_sets_arc);
    for edge in receiver {
        graph[edge.0].push(edge.1);
        graph[edge.1].push(edge.0);
    }
    graph
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

fn bfs_helper(
    graph: &Vec<Vec<usize>>,
    start_node: usize,
    explored: &mut HashSet<usize>,
) -> Vec<usize> {
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
