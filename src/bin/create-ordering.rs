use clap::Parser;
use log::{debug, info};
use musk::{io::{dump_data_to_file, load_data_from_file}, utility::{average_hamming_distance, find_ordering}};
use std::path::Path;

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

    let mut all_distances: Vec<(Vec<u32>, String, u32)> = vec![];
    for index in 0..distances.len() {
        let mut full_distance_vector = vec![];
        for prior_index in 0..index {
            full_distance_vector.push(all_distances[prior_index].0[index])
        }
        full_distance_vector.push(0);
        full_distance_vector.append(&mut distances[index].0);
        all_distances.push((
            full_distance_vector,
            distances[index].1.clone(),
            distances[index].2.clone(),
        ));
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
