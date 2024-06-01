use clap::Parser;
use indicatif::{ParallelProgressIterator, ProgressIterator};
use itertools::Itertools;
use log::info;
use musk::io::load_string2taxid;
use musk::rle::{BuildRunLengthEncoding, RunLengthEncoding};
use musk::utility::{create_bitmap, get_range};
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::path::Path;

/// Creates a file to tax id mapping where files with the same tax id are grouped
/// together if their k-mer spectra are similar enough.
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 14)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg(short, long, default_value_t = 0)]
    /// 2^{log_blocks} partitions
    log_blocks: u32,

    #[arg(short, long, default_value_t = 0)]
    /// The index of the block to use
    block_index: usize,

    #[arg(short, long)]
    /// The old directory prefix of the fasta files
    old_directory_prefix: Option<String>,

    #[arg(short, long)]
    /// The new directory prefix of the fasta files
    new_directory_prefix: Option<String>,

    #[arg()]
    /// The ordering of the sequences for the full matrix
    ordering_file: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let ordering_file_path = Path::new(&args.ordering_file);

    let mut ordering = load_string2taxid(ordering_file_path);
    if let (Some(old_prefix), Some(new_prefix)) =
        (args.old_directory_prefix, args.new_directory_prefix)
    {
        ordering = ordering
            .into_iter()
            .map(|(files, taxid)| (files.replace(&*old_prefix, &*new_prefix), taxid))
            .collect_vec();
    }

    info!("creating roaring bitmaps for each group...");

    let (lowest_kmer, highest_kmer) =
        get_range(args.kmer_length, args.log_blocks, args.block_index);
    let bitmaps = ordering
        .into_par_iter()
        .progress()
        .map(|(files, _taxid)| {
            (
                files.clone(),
                create_bitmap(
                    &*files,
                    args.kmer_length,
                    lowest_kmer,
                    highest_kmer,
                    false,
                    false,
                ),
            )
        })
        .collect::<Vec<(String, RoaringBitmap)>>();

    info!("roaring bitmaps computed, creating database...");

    let mut database = vec![BuildRunLengthEncoding::new(); 4_usize.pow(args.kmer_length as u32)];
    for (index, (_files, bitmap)) in bitmaps.into_iter().progress().enumerate() {
        for kmer in bitmap {
            database[kmer as usize].push(index);
        }
    }
    let naive_database_runs = database
        .iter()
        .map(|build_rle| build_rle.get_vector().len())
        .sum::<usize>();
    let compressed_database = database
        .into_par_iter()
        .map(|build_rle| build_rle.to_rle())
        .collect::<Vec<RunLengthEncoding>>();
    let compressed_database_runs = compressed_database
        .iter()
        .map(|rle| rle.get_vector().len())
        .sum::<usize>();

    info!(
        "{}\t{}\t{}\t{}",
        args.log_blocks, args.block_index, compressed_database_runs, naive_database_runs
    );

    info!("done!");
}
