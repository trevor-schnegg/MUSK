use clap::Parser;
use indicatif::{ParallelProgressIterator, ProgressIterator};
use itertools::Itertools;
use log::info;
use musk::io::load_data_from_file;
use musk::rle::{BuildRunLengthEncoding, RunLengthEncoding};
use musk::utility::{create_bitmap, get_range};
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::collections::HashMap;
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
    /// The ordering of sequences for the block
    block_ordering_file: String,

    #[arg()]
    /// The ordering of the sequences for the full matrix
    matrix_ordering_file: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let block_ordering_file_path = Path::new(&args.block_ordering_file);
    let matrix_ordering_file_path = Path::new(&args.matrix_ordering_file);

    info!("loading ordering at {}", args.block_ordering_file);

    let mut matrix_ordering = load_data_from_file::<Vec<(String, u32)>>(matrix_ordering_file_path);
    if let (Some(old_prefix), Some(new_prefix)) =
        (args.old_directory_prefix, args.new_directory_prefix)
    {
        matrix_ordering = matrix_ordering
            .into_iter()
            .map(|(files, taxid)| (files.replace(&*old_prefix, &*new_prefix), taxid))
            .collect_vec();
    }
    let block_ordering = load_data_from_file::<Vec<(String, u32)>>(block_ordering_file_path);

    info!("creating roaring bitmaps for each group...");

    let (lowest_kmer, highest_kmer) =
        get_range(args.kmer_length, args.log_blocks, args.block_index);
    let bitmaps: HashMap<String, RoaringBitmap> = HashMap::from_iter(
        matrix_ordering
            .par_iter()
            .progress()
            .map(|(files, _taxid)| {
                (
                    files.clone(),
                    create_bitmap(files.clone(), args.kmer_length, lowest_kmer, highest_kmer),
                )
            })
            .collect::<Vec<(String, RoaringBitmap)>>()
            .into_iter(),
    );

    info!("roaring bitmaps computed, creating block ordering database...");

    let mut block_ordering_database =
        vec![BuildRunLengthEncoding::new(); 4_usize.pow(args.kmer_length as u32)];
    for (index, (files, _taxid)) in block_ordering.into_iter().progress().enumerate() {
        let bitmap = bitmaps.get(&files).unwrap();
        for kmer in bitmap {
            block_ordering_database[kmer as usize].push(index);
        }
    }
    let naive_block_ordering_runs = block_ordering_database
        .iter()
        .map(|build_rle| build_rle.get_vector().len())
        .sum::<usize>();
    let block_ordering_compressed_database = block_ordering_database
        .into_par_iter()
        .map(|build_rle| build_rle.to_rle())
        .collect::<Vec<RunLengthEncoding>>();
    let block_ordering_compressed_runs = block_ordering_compressed_database
        .iter()
        .map(|rle| rle.get_vector().len())
        .sum::<usize>();

    info!("block ordering database completed, creating full matrix ordering database...");

    let mut matrix_ordering_database =
        vec![BuildRunLengthEncoding::new(); 4_usize.pow(args.kmer_length as u32)];
    for (index, (files, _taxid)) in matrix_ordering.into_iter().progress().enumerate() {
        let bitmap = bitmaps.get(&files).unwrap();
        for kmer in bitmap {
            matrix_ordering_database[kmer as usize].push(index);
        }
    }
    let naive_matrix_ordering_runs = matrix_ordering_database
        .iter()
        .map(|build_rle| build_rle.get_vector().len())
        .sum::<usize>();
    let matrix_ordering_compressed_database = matrix_ordering_database
        .into_par_iter()
        .map(|build_rle| build_rle.to_rle())
        .collect::<Vec<RunLengthEncoding>>();
    let matrix_ordering_compressed_runs = matrix_ordering_compressed_database
        .iter()
        .map(|rle| rle.get_vector().len())
        .sum::<usize>();

    println!(
        "{}\t{}\t{}\t{}\t{}\t{}",
        args.log_blocks,
        args.block_index,
        block_ordering_compressed_runs,
        matrix_ordering_compressed_runs,
        naive_block_ordering_runs,
        naive_matrix_ordering_runs
    );

    info!("done!");
}
