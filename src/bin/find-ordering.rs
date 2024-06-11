use clap::Parser;
use concorde_rs::{solver, LowerDistanceMatrix};
use itertools::Itertools;
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
    #[arg(short, long, action)]
    /// Flag that specifies using the Held-Karp algorithm to solve TSP
    held_karp: bool,

    #[arg(short, long, action)]
    /// Flag that specifies using the Lin-Kernighan algorithm to solve TSP
    lin_kernighan: bool,

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

    let ordering = if args.lin_kernighan || args.held_karp {
        // Flatten the distance matrix
        let dist_mat = LowerDistanceMatrix::new(
            distances.len() as u32,
            distances.iter().flat_map(|row| row.clone()).collect_vec(),
        );

        if args.lin_kernighan {
            let solution = solver::tsp_lk(&dist_mat).unwrap();
            debug!("length of tour: {}", solution.length);
            solution.tour.iter().map(|x| *x as usize).collect_vec()
        } else {
            let solution = solver::tsp_hk(&dist_mat).unwrap();
            debug!("length of tour: {}", solution.length);
            solution.tour.iter().map(|x| *x as usize).collect_vec()
        }
    } else {
        // Perform the greedy solution
        let greedy_ordering = greedy_ordering(&distances, args.start);
        let avg_dist_output = average_hamming_distance(&greedy_ordering, &distances);
        debug!("length of tour: {}", avg_dist_output.1);
        greedy_ordering
    };

    for index in ordering {
        let (files_string, taxid) = &file2taxid[index];
        output_file
            .write(format!("{}\t{}\n", *files_string, *taxid).as_bytes())
            .expect("could not write to output file");
    }
}
