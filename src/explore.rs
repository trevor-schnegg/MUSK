use roaring::RoaringBitmap;
use std::cmp::min;
use std::collections::{HashSet, VecDeque};
use std::sync::{mpsc, Arc};
use threadpool::ThreadPool;

pub fn connected_components(
    bit_vectors: Vec<RoaringBitmap>,
    minimum_similarity: f64,
    thread_number: usize,
) -> Vec<Vec<usize>> {
    let graph = create_graph(bit_vectors, minimum_similarity, thread_number);
    let components = bfs(graph);
    components
}

fn create_graph(
    bitmaps: Vec<RoaringBitmap>,
    minimum_similarity: f64,
    thread_number: usize,
) -> Vec<Vec<usize>> {
    let mut graph = vec![vec![]; bitmaps.len()];
    let bitmaps_arc = Arc::new(bitmaps);
    let (sender, receiver) = mpsc::channel();
    let pool = ThreadPool::new(thread_number);
    for i1 in 0..bitmaps_arc.len() {
        let sender_clone = sender.clone();
        let bitmaps_arc_clone = bitmaps_arc.clone();
        pool.execute(move || {
            for i2 in 0..bitmaps_arc_clone.len() {
                if i2 <= i1 {
                    continue;
                }
                let (bit_vector_1, bit_vector_2) = (&bitmaps_arc_clone[i1], &bitmaps_arc_clone[i2]);
                let intersect_size = intersect(bit_vector_1, bit_vector_2);
                let minimum_containment =
                    intersect_size as f64 / min(bit_vector_1.len(), bit_vector_2.len()) as f64;
                if minimum_containment >= minimum_similarity {
                    sender_clone.send((i1, i2)).unwrap();
                }
            }
        });
    }
    drop(sender);
    drop(bitmaps_arc);
    for edge in receiver {
        graph[edge.0].push(edge.1);
        graph[edge.1].push(edge.0);
    }
    graph
}

fn intersect(vector_1: &RoaringBitmap, vector_2: &RoaringBitmap) -> u64 {
    vector_1.intersection_len(vector_2)
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
