use clap::Parser;
use itertools::Itertools;
use log::{debug, info};
use musk::io::{dump_data_to_file, load_data_from_file};
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
    block_i: usize,

    #[arg(short, long)]
    /// The directory prefix of the fasta files
    old_directory_prefix: Option<String>,

    #[arg(short, long)]
    /// The directory prefix of the fasta files
    new_directory_prefix: Option<String>,

    #[arg()]
    /// the file2taxid file
    ordering_file: String,

    #[arg()]
    /// Location to output the serialzed distances
    output_file: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let ordering_file_path = Path::new(&args.ordering_file);
    let output_file_path = Path::new(&args.output_file);

    info!("loading ordering at {}", args.ordering_file);
    let mut ordering = load_data_from_file::<Vec<(String, u32)>>(ordering_file_path);
    if let (Some(old_prefix), Some(new_prefix)) =
        (args.old_directory_prefix, args.new_directory_prefix)
    {
        ordering = ordering
            .into_iter()
            .map(|(files, taxid)| (files.replace(&*old_prefix, &*new_prefix), taxid))
            .collect_vec();
    }
    info!("creating roaring bitmaps for each group...");
    let (lowest_kmer, highest_kmer) = get_range(args.kmer_length, args.log_blocks, args.block_i);
    let bitmaps = ordering
        .into_par_iter()
        .map(|(files, _taxid)| {
            create_bitmap(files, args.kmer_length, lowest_kmer, highest_kmer)
        })
        .collect::<Vec<RoaringBitmap>>();
    info!("kmer vectors computed, creating database...");

    let mut database = vec![BuildRunLengthEncoding::new(); 4_usize.pow(args.kmer_length as u32)];
    for (index, bitmap) in bitmaps.into_iter().enumerate() {
        for kmer in bitmap {
            database[kmer as usize].push(index);
        }
        if index % 1000 == 0 && index != 0 {
            debug!("done inserting {} bitmaps into the database", index);
        }
    }
    let naive_runs = database
        .iter()
        .map(|build_rle| build_rle.get_vector().len())
        .sum::<usize>();
    info!("Total naive runs for the ordering {}", naive_runs);
    let compressed_database = database
        .into_par_iter()
        .map(|build_rle| build_rle.to_rle())
        .collect::<Vec<RunLengthEncoding>>();
    let compressed_runs = compressed_database
        .iter()
        .map(|rle| rle.get_vector().len())
        .sum::<usize>();
    info!("Total compressed runs for the ordering {}", compressed_runs);
    info!("Saving the compressed runs to the output file...");
    dump_data_to_file(
        bincode::serialize(&compressed_database).unwrap(),
        output_file_path,
    )
    .unwrap();
}
