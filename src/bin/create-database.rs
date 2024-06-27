use clap::Parser;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use musk::database::Database;
use musk::io::{dump_data_to_file, load_string2taxid};
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::create_bitmap;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::ops::Neg;
use std::path::Path;
use tracing::info;

/// Creates a run length encoding database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, action)]
    /// Flag that specifies whether or not to use canonical kmers
    canonical: bool,

    #[arg(short, long, default_value_t = 18)]
    /// The exponent e for the significance of hits
    /// Used in the equation 10^{-e} to determine statistical significance
    cutoff_threshold_exp: i32,

    #[arg(short, long, default_value_t = 14)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg(short, long, default_value_t = 150)]
    /// Number of queries to sample
    num_queries: u64,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// Directory to output the file2taxid file
    output_directory: String,

    #[arg(short, long)]
    /// Remove runs with count of ones <= to the input number
    remove_runs: Option<usize>,

    #[arg()]
    /// The ordered file2taxid of the sequences
    ordering_file: String,

    #[arg()]
    /// Directory with fasta files to create reference from
    reference_directory: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_musk_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let cutoff_threshold = 10.0_f64.powi((args.cutoff_threshold_exp).neg());
    let ordering_file_path = Path::new(&args.ordering_file);
    let output_dir_path = Path::new(&args.output_directory);
    let reference_dir_path = Path::new(&args.reference_directory);

    // Load the file2taxid ordering
    let file2taxid_ordering = load_string2taxid(ordering_file_path);

    info!("creating roaring bitmaps for each group...");

    let bitmaps = file2taxid_ordering
        .par_iter()
        .progress()
        .map(|(files, _taxid)| {
            let file_paths = files
                .split("$")
                .map(|file| reference_dir_path.join(file))
                .collect_vec();

            create_bitmap(file_paths, args.kmer_length, args.canonical)
        })
        .collect::<Vec<RoaringBitmap>>();

    info!("raring bitmaps created! constructing database...");

    let mut database = Database::from(
        bitmaps,
        args.canonical,
        cutoff_threshold,
        file2taxid_ordering,
        args.kmer_length,
        args.num_queries,
    );

    match args.remove_runs {
        None => {}
        Some(num_ones) => {
            database.lossy_compression(num_ones);
        }
    }

    dump_data_to_file(
        bincode::serialize(&database).unwrap(),
        &output_dir_path.join("musk.db"),
    )
    .unwrap();

    info!("done!");
}
