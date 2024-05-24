use clap::Parser;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use log::info;
use musk::io::{dump_data_to_file, load_string2taxid};
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

    #[arg(short, long, default_value_t = 12)]
    /// Length of k-mer to use in the database
    thread_number: usize,

    #[arg(short, long, default_value_t = 0)]
    /// 2^{log_blocks} partitions
    log_blocks: u32,

    #[arg(short, long, default_value_t = 0)]
    /// The index of the block to use
    block_index: usize,

    #[arg(short, long)]
    /// The directory prefix of the fasta files
    old_directory_prefix: Option<String>,

    #[arg(short, long)]
    /// The directory prefix of the fasta files
    new_directory_prefix: Option<String>,

    #[arg()]
    /// Location to output the serialzed distances
    output_file: String,

    #[arg()]
    /// the file2taxid file
    file2taxid: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let file2taxid_path = Path::new(&args.file2taxid);
    let output_file_path = Path::new(&args.output_file);

    info!("loading files2taxid at {}", args.file2taxid);
    let mut file2taxid = load_string2taxid(file2taxid_path);
    if let (Some(old_prefix), Some(new_prefix)) =
        (args.old_directory_prefix, args.new_directory_prefix)
    {
        file2taxid = file2taxid
            .into_iter()
            .map(|(files, taxid)| (files.replace(&*old_prefix, &*new_prefix), taxid))
            .collect_vec();
    }
    info!("{} groups total", file2taxid.len());
    info!("creating roaring bitmaps for each group...");
    let (lowest_kmer, highest_kmer) =
        get_range(args.kmer_length, args.log_blocks, args.block_index);
    let bitmaps = file2taxid
        .into_par_iter()
        .progress()
        .map(|(files, taxid)| {
            (
                create_bitmap(files.clone(), args.kmer_length, lowest_kmer, highest_kmer),
                files,
                taxid,
            )
        })
        .collect::<Vec<(RoaringBitmap, String, u32)>>();
    info!("roaring bitmaps computed, creating distance matrix...");
    let all_distances = bitmaps
        .par_iter()
        .progress()
        .enumerate()
        .map(|(index_1, (bitmap_1, files_1, taxid_1))| {
            let inner_distances = bitmaps
                .par_iter()
                .enumerate()
                .filter_map(|(index_2, (bitmap_2, _files_2, _taxid_2))| {
                    if index_2 <= index_1 {
                        None
                    } else {
                        let intersection_size = bitmap_1.intersection_len(bitmap_2);
                        // |A| + |B| - (2 * |A & B|)
                        let distance =
                            (bitmap_1.len() + bitmap_2.len() - (2 * intersection_size)) as u32;
                        Some(distance)
                    }
                })
                .collect::<Vec<u32>>();
            (inner_distances, files_1.clone(), *taxid_1)
        })
        .collect::<Vec<(Vec<u32>, String, u32)>>();
    info!("distance matrix completed! outputting to file");

    dump_data_to_file(
        bincode::serialize(&all_distances).unwrap(),
        output_file_path,
    )
    .unwrap();
}
