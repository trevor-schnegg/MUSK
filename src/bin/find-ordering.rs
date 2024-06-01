use clap::Parser;
use concorde_rs::{solver, LowerDistanceMatrix};
use itertools::Itertools;
use log::{debug, info};
use musk::{
    io::load_data_from_file,
    utility::{average_hamming_distance, greedy_ordering},
};
use std::path::Path;

/// Creates an ordering of files based on distances between bitmaps
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 0)]
    /// Start index of the naive shortest path traversal
    start: usize,

    #[arg(short, long, action)]
    /// Flag that specifies using the Held-Karp algorithm to solve TSP
    held_karp: bool,

    #[arg(short, long, action)]
    /// Flag that specifies using the Lin-Kernighan algorithm to solve TSP
    lin_kernighan: bool,

    #[arg()]
    /// the distances file
    distances: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let distances_file = Path::new(&args.distances);

    info!("loading distances at {}", args.distances);
    let distances = load_data_from_file::<Vec<(Vec<u32>, String, u32)>>(distances_file);
    debug!("length of distances: {}", distances.len());
    info!("distances loaded!, filling out matrix...");

    let ordering = if args.lin_kernighan {
        let dist_mat = LowerDistanceMatrix::new(
            distances.len() as u32,
            distances
                .iter()
                .flat_map(|tuple| tuple.0.clone())
                .collect_vec(),
        );
        let solution = solver::tsp_lk(&dist_mat).unwrap();
        debug!("length of tour: {}", solution.length);
        solution.tour.iter().map(|x| *x as usize).collect_vec()
    } else if args.held_karp {
        let dist_mat = LowerDistanceMatrix::new(
            distances.len() as u32,
            distances
                .iter()
                .flat_map(|tuple| tuple.0.clone())
                .collect_vec(),
        );
        let solution = solver::tsp_hk(&dist_mat).unwrap();
        debug!("length of tour: {}", solution.length);
        solution.tour.iter().map(|x| *x as usize).collect_vec()
    } else {
        let greedy_ordering = greedy_ordering(&distances, args.start);
        let avg_dist_output = average_hamming_distance(&greedy_ordering, &distances);
        debug!("length of tour: {}", avg_dist_output.1);
        greedy_ordering
    };

    for (files, taxid) in ordering
        .into_iter()
        .map(|x| (distances[x].1.clone(), distances[x].2))
    {
        println!("{}\t{}", files, taxid);
    }
}
