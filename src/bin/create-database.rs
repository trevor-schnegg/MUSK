use clap::Parser;
use indicatif::{ParallelProgressIterator, ProgressIterator};
use itertools::Itertools;
use musk::io::{dump_data_to_file, load_string2taxid};
use musk::rle::{NaiveRunLengthEncoding, RunLengthEncoding};
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::{create_bitmap, get_range};
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::path::Path;
use tracing::info;

/// Creates a run length encoding database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 0)]
    /// The index of the block to use
    block_index: usize,

    #[arg(short, long, action)]
    /// Flag that specifies whether or not to use canonical kmers
    canonical: bool,

    #[arg(short, long, default_value_t = 14)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg(short, long, default_value_t = 0)]
    /// 2^{log_blocks} partitions
    log_blocks: u32,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// Directory to output the file2taxid file
    output_directory: String,

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
    let ordering_file_path = Path::new(&args.ordering_file);
    let output_dir_path = Path::new(&args.output_directory);
    let reference_dir_path = Path::new(&args.reference_directory);

    let total_num_kmers = 4_usize.pow(args.kmer_length as u32) as f64;

    let file2taxid_ordering = load_string2taxid(ordering_file_path);

    info!("creating roaring bitmaps for each group...");

    let (lowest_kmer, highest_kmer) =
        get_range(args.kmer_length, args.log_blocks, args.block_index);

    let bitmaps = file2taxid_ordering
        .par_iter()
        .progress()
        .map(|(files, _taxid)| {
            let file_paths = files
                .split("$")
                .map(|file| reference_dir_path.join(file))
                .collect_vec();

            create_bitmap(
                file_paths,
                args.kmer_length,
                lowest_kmer,
                highest_kmer,
                false,
                args.canonical,
            )
        })
        .collect::<Vec<RoaringBitmap>>();

    let p_values = bitmaps.par_iter().map(|bitmap| bitmap.len() as f64 / total_num_kmers).collect::<Vec<f64>>();

    info!("roaring bitmaps computed, creating database...");

    let mut database =
        vec![NaiveRunLengthEncoding::new(); (4 as usize).pow(args.kmer_length as u32)];

    for (index, bitmap) in bitmaps.into_iter().progress().enumerate() {
        for kmer in bitmap {
            database[kmer as usize].push(index);
        }
    }

    let naive_database_runs = database
        .iter()
        .map(|build_rle| build_rle.get_raw_runs().len())
        .sum::<usize>();
    let compressed_database = database
        .into_par_iter()
        .map(|build_rle| build_rle.to_rle())
        .collect::<Vec<RunLengthEncoding>>();
    let compressed_database_runs = compressed_database
        .iter()
        .map(|rle| rle.get_raw_runs().len())
        .sum::<usize>();

    info!(
        "{}\t{}\t{}\t{}",
        args.log_blocks, args.block_index, compressed_database_runs, naive_database_runs
    );

    dump_data_to_file(
        bincode::serialize(&(compressed_database, file2taxid_ordering, p_values)).unwrap(),
        output_dir_path,
    )
    .unwrap();

    info!("done!");
}
