use clap::Parser;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use musk::io::{dump_data_to_file, load_string2taxid};
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::create_bitmap;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::path::Path;
use tracing::info;

/// Computes the lower triangle of a pairwise distance matrix from the input sequences (or sequence groups)
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// Directory to output the file2taxid file
    output_directory: String,

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

    #[arg()]
    /// the file2taxid file
    file2taxid: String,

    #[arg()]
    /// Directory with fasta files to create reference from
    reference_directory: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_musk_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let file2taxid_path = Path::new(&args.file2taxid);
    let output_dir_path = Path::new(&args.output_directory);
    let reference_dir_path = Path::new(&args.reference_directory);

    info!("loading files2taxid at {}", args.file2taxid);

    let file2taxid = load_string2taxid(file2taxid_path);

    info!(
        "{} groups total, creating roaring bitmaps for each group...",
        file2taxid.len()
    );

    let bitmaps = file2taxid
        .par_iter()
        .progress()
        .map(|(files, _taxid)| {
            let file_paths = files
                .split("$")
                .map(|file| reference_dir_path.join(file))
                .collect_vec();

            create_bitmap(file_paths, args.kmer_length, false, args.canonical)
        })
        .collect::<Vec<RoaringBitmap>>();

    info!("roaring bitmaps computed, creating distance matrix...");

    let distances = bitmaps
        .par_iter()
        .progress()
        .enumerate()
        .map(|(index_1, bitmap_1)| {
            bitmaps[..=index_1]
                .iter()
                .enumerate()
                .map(|(index_2, bitmap_2)| {
                    if index_1 == index_2 {
                        0
                    } else {
                        let intersection_size = bitmap_1.intersection_len(bitmap_2);
                        // |A| + |B| - (2 * |A & B|)
                        (bitmap_1.len() + bitmap_2.len() - (2 * intersection_size)) as u32
                    }
                })
                .collect::<Vec<u32>>()
        })
        .collect::<Vec<Vec<u32>>>();

    info!("distance matrix completed! outputting to file...");

    dump_data_to_file(
        bincode::serialize(&(distances, file2taxid)).unwrap(),
        &output_dir_path.join("musk.pairwise_distances"),
    )
    .unwrap();

    info!("done!");
}
