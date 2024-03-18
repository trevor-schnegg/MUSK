use clap::Parser;
use log::{debug, info};
use musk::io::{dump_data_to_file, load_data_from_file};
use rand::seq::SliceRandom;
use rand::thread_rng;
use std::collections::HashSet;
use std::path::Path;

fn find_ordering(distances: &Vec<(Vec<u64>, String, u32)>, start_index: usize) -> Vec<usize> {
    let mut connected_indices = HashSet::from([start_index]);
    let mut ordering = vec![start_index];
    let mut current_index = start_index;
    while ordering.len() < distances.len() {
        let mut next_index = 0_usize;
        let mut next_index_minimum = u64::MAX;
        for (index, distance) in distances[current_index].0.iter().enumerate() {
            if connected_indices.contains(&index) {
                continue;
            }
            if *distance < next_index_minimum {
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
    distances: &Vec<(Vec<u64>, String, u32)>,
) -> (f64, u64) {
    let sum = ordering.windows(2).map(|x| distances[x[0]].0[x[1]]).sum();
    (sum as f64 / (ordering.len() - 1) as f64, sum)
}

/// Creates an ordering of files based on distances between bitmaps
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 0)]
    /// Start index of the naive shortest walk traversal
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
    let mut distances = load_data_from_file::<Vec<(Vec<u64>, String, u32)>>(distances_file);
    debug!("length of distances: {}", distances.len());
    info!("distances loaded! filling out the matrix of distance values...");

    // Fill out the remaining distances
    for index in 0..distances.len() {
        let mut prior_distances = vec![];
        for prior_index in 0..index {
            let distance_between = distances[prior_index].0[index];
            prior_distances.push(distance_between);
        }
        prior_distances.push(0);
        prior_distances.append(&mut distances[index].0);
        distances[index].0 = prior_distances;
    }
    info!("full matrix constructed! finding ordering...");

    let ordering = find_ordering(&distances, args.start);
    assert_eq!(ordering[0], args.start);
    let avg_dist_output = average_hamming_distance(&ordering, &distances);
    debug!(
        "average hamming distance of ordering: {} (total: {})",
        avg_dist_output.0, avg_dist_output.1
    );

    let mut generic_ordering = (0..distances.len()).collect::<Vec<usize>>();
    generic_ordering.shuffle(&mut thread_rng());
    let generic_average_hamming_distance = average_hamming_distance(&generic_ordering, &distances);
    debug!(
        "average hamming distance of ordering: {} (total: {})",
        generic_average_hamming_distance.0, generic_average_hamming_distance.1
    );

    let data_to_output = ordering
        .into_iter()
        .map(|x| (distances[x].1.clone(), distances[x].2))
        .collect::<Vec<(String, u32)>>();
    dump_data_to_file(
        bincode::serialize(&data_to_output).unwrap(),
        output_file_path,
    )
    .unwrap();
}
