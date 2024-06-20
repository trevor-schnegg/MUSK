use clap::Parser;
use musk::{
    io::load_data_from_file,
    tracing::start_musk_tracing_subscriber,
    utility::{average_hamming_distance, greedy_ordering},
};
use std::{fs::File, io::Write, path::Path};
use tracing::{debug, info};

/// Creates an ordering of files based on a pairwise distance matrix
/// This is done such that the total hamming distance of the ordering is as small as possible
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// Directory to output the file2taxid file
    output_directory: String,

    #[arg(short, long, default_value_t = 0)]
    /// Start index of the naive shortest path traversal
    start: usize,

    #[arg()]
    /// The pairwise distances file
    distances: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_musk_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let distances_file = Path::new(&args.distances);
    let output_dir_path = Path::new(&args.output_directory);

    let mut output_file = File::create(output_dir_path.join("musk.ordered.file2taxid"))
        .expect("could not create output file");

    info!("loading distances at {}", args.distances);

    let (distances, file2taxid) =
        load_data_from_file::<(Vec<Vec<u32>>, Vec<(String, usize)>)>(distances_file);

    debug!("length of distances: {}", distances.len());
    info!("distances loaded! finding ordering...");

    // Perform the greedy solution
    let greedy_ordering = greedy_ordering(&distances, args.start);
    let avg_dist_output = average_hamming_distance(&greedy_ordering, &distances);
    debug!("length of tour: {}", avg_dist_output.1);

    for index in greedy_ordering {
        let (files_string, taxid) = &file2taxid[index];
        output_file
            .write(format!("{}\t{}\n", *files_string, *taxid).as_bytes())
            .expect("could not write to output file");
    }
}
