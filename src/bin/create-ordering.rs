use clap::Parser;
use log::{debug, info};
use musk::io::{dump_data_to_file, load_data_from_file};
use std::collections::HashSet;
use std::path::Path;

fn find_ordering(distances: &Vec<(Vec<u32>, String, u32)>, start_index: usize) -> Vec<usize> {
    let mut connected_indices = HashSet::from([start_index]);
    let mut ordering = vec![start_index];
    let mut current_index = start_index;
    while ordering.len() < distances.len() {
        let mut next_index = 0_usize;
        let mut next_index_minimum = u32::MAX;
        for (index, distance) in distances[current_index].0.iter().enumerate() {
            if *distance < next_index_minimum && !connected_indices.contains(&index) {
                next_index = index;
                next_index_minimum = *distance;
            }
        }
        ordering.push(next_index);
        connected_indices.insert(next_index);
        current_index = next_index;
        if ordering.len() % 1000 == 0 {
            debug!("found ordering for {} bitmaps", ordering.len());
        }
    }
    ordering
}

fn average_hamming_distance(
    ordering: &Vec<usize>,
    distances: &Vec<(Vec<u32>, String, u32)>,
) -> (f64, u64) {
    let sum = ordering
        .windows(2)
        .map(|x| distances[x[0]].0[x[1]] as u64)
        .sum();
    (sum as f64 / (ordering.len() - 1) as f64, sum)
}

/// Creates an ordering of files based on distances between bitmaps
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 0)]
    /// Start index of the naive shortest path traversal
    start: usize,

    #[arg()]
    /// the distances file
    distances: String,

    #[arg()]
    /// location to output the serialized ordering
    output_file: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let distances_file = Path::new(&args.distances);
    let output_file_path = Path::new(&args.output_file);

    info!("loading distances at {}", args.distances);
    let mut distances = load_data_from_file::<Vec<(Vec<u32>, String, u32)>>(distances_file);
    debug!("length of distances: {}", distances.len());
    info!("distances loaded!, filling out matrix...");

    let mut all_distances = vec![];
    for index in 0..distances.len() {
        let mut full_distance_vector = vec![];
        for other_index in 0..index {
            full_distance_vector.push(distances[other_index].0[index - (other_index + 1)])
        }
        full_distance_vector.push(0);
        full_distance_vector.append(&mut distances[index].0);
        all_distances.push((full_distance_vector, distances[index].1.clone(), distances[index].2.clone()));
    }

    let ordering = find_ordering(&all_distances, args.start);
    let avg_dist_output = average_hamming_distance(&ordering, &all_distances);
    debug!(
        "average hamming distance of ordering: {} (total: {})",
        avg_dist_output.0, avg_dist_output.1
    );

    let ordering_output = ordering
        .into_iter()
        .map(|x| (distances[x].1.clone(), distances[x].2))
        .collect::<Vec<(String, u32)>>();
    dump_data_to_file(
        bincode::serialize(&ordering_output).unwrap(),
        output_file_path,
    )
    .unwrap();
}
