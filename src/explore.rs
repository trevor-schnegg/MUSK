use std::cmp::min;
use std::collections::{HashSet, VecDeque};

pub fn connected_components(kmer_sets: Vec<HashSet<usize>>, minimum_similarity: f64) -> Vec<Vec<usize>> {
    let graph = create_graph(kmer_sets, minimum_similarity);
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

fn create_graph(kmer_sets: Vec<HashSet<usize>>, minimum_similarity: f64) -> Vec<Vec<usize>> {
    let mut graph = vec![vec![]; kmer_sets.len()];
    for (i1, kmer_set_1) in kmer_sets.iter().enumerate() {
        for (i2, kmer_set_2) in kmer_sets.iter().enumerate() {
            if i2 <= i1 {
                continue;
            }
            let intersect = kmer_set_1.intersection(&kmer_set_2).map(|x| *x).collect::<Vec<usize>>().len() as f64;
            let min_coverage = intersect / min(kmer_set_1.len(), kmer_set_2.len()) as f64;
            if min_coverage >= minimum_similarity {
                graph[i1].push(i2);
                graph[i2].push(i1);
            }
        }
    }
    graph
}