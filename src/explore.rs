use indicatif::ParallelProgressIterator;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator};
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::collections::{HashSet, VecDeque};

pub fn connected_components(
    bitmaps: Vec<RoaringBitmap>,
    minimum_similarity: f64,
) -> Vec<Vec<usize>> {
    let graph = create_graph(bitmaps);
    let components = bfs(graph, minimum_similarity);
    components
}

fn create_graph(bitmaps: Vec<RoaringBitmap>) -> Vec<Vec<f64>> {
    bitmaps
        .par_iter()
        .progress()
        .enumerate()
        .map(|(index_1, bitmap_1)| {
            let mut similarities = bitmaps[..index_1]
                .par_iter()
                .map(|bitmap_2| {
                    let intersection_size = bitmap_1.intersection_len(bitmap_2);
                    let union_size = bitmap_1.union_len(bitmap_2);
                    intersection_size as f64 / union_size as f64
                })
                .collect::<Vec<f64>>();
            similarities.push(0.0);
            similarities
        })
        .collect::<Vec<Vec<f64>>>()
}

/// Returns the connected components of all nodes
fn bfs(graph: Vec<Vec<f64>>, minimum_similarity: f64) -> Vec<Vec<usize>> {
    let mut explored = HashSet::new();
    let mut connected_components = Vec::new();
    for s in 0..graph.len() {
        if explored.contains(&s) {
            continue;
        }
        connected_components.push(bfs_helper(&graph, s, &mut explored, minimum_similarity));
    }
    connected_components
}

fn bfs_helper(
    graph: &Vec<Vec<f64>>,
    start_node: usize,
    explored: &mut HashSet<usize>,
    minimum_similarity: f64,
) -> Vec<usize> {
    explored.insert(start_node);
    let mut queue = VecDeque::from([start_node]);
    let mut connected_component = Vec::from([start_node]);
    while !queue.is_empty() {
        let node = queue.pop_front().unwrap();
        for (index, similarity) in graph[node][..node]
            .iter()
            .chain(graph[node..].iter().map(|vec| &vec[node]))
            .enumerate()
        {
            if explored.contains(&index) {
                continue;
            } else if minimum_similarity <= *similarity {
                queue.push_back(index);
                explored.insert(index);
                connected_component.push(index);
            }
        }
    }
    connected_component
}
