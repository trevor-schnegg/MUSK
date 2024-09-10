use clap::Parser;
use musk::{
    io::{create_output_file, load_data_from_file},
    tracing::start_musk_tracing_subscriber,
    utility::{greedy_ordering, ordering_statistics},
};
use std::{
    io::{BufWriter, Write},
    path::Path,
};
use tracing::{debug, info};

/// Creates an ordering of files based on a pairwise distance matrix
/// This is done such that the total hamming distance of the ordering is as small as possible
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// Where to write the output
    /// If a file, '.musk.o.f2t' is added
    /// If a directory, 'musk.o.f2t' will be the file name
    /// Name means: musk, (o)rdered, (f)ile(2)(t)axid
    output_location: String,

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
    let output_loc_path = Path::new(&args.output_location);

    // Create the output file
    let mut output_file = BufWriter::new(create_output_file(output_loc_path, "musk.o.f2t"));

    info!("loading distances at {}", args.distances);
    let (distances, file2taxid) =
        load_data_from_file::<(Vec<Vec<u32>>, Vec<(String, usize)>)>(distances_file);

    debug!("length of distances: {}", distances.len());
    info!("distances loaded! finding ordering...");

    // Perform the greedy solution
    let greedy_ordering = greedy_ordering(&distances, args.start);
    let (_avg_dist, total_dist) = ordering_statistics(&greedy_ordering, &distances);
    debug!("length of tour: {}", total_dist);

    for index in greedy_ordering {
        let (files_string, taxid) = &file2taxid[index];
        output_file
            .write(format!("{}\t{}\n", *files_string, *taxid).as_bytes())
            .expect("could not write to output file");
    }
}
