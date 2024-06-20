use clap::Parser;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use musk::io::load_string2taxid;
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::{create_bitmap, get_range};
use rayon::iter::IntoParallelIterator;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::collections::HashMap;
use std::path::Path;
use tracing::info;

/// Prints to stdout a map in the form of <fasta-file-path>\t<tax-id> given a reference location
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 14)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg()]
    /// Location of the files2taxid file
    files2taxid: String,

    #[arg()]
    /// Directory with fasta files to create reference from
    reference_directory: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_musk_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let files2taxid = Path::new(&args.files2taxid);
    let reference_dir_path = Path::new(&args.reference_directory);

    let mut kmer_counts = vec![0_usize; 4_usize.pow(args.kmer_length as u32)];

    info!("getting kmer counts...");

    let bitmaps = load_string2taxid(files2taxid)
        .into_par_iter()
        .progress()
        .map(|(files, _taxid)| {
            let file_paths = files
                .split("$")
                .map(|file| reference_dir_path.join(file))
                .collect_vec();

            create_bitmap(
                file_paths,
                args.kmer_length,
                false,
                false,
            )
        })
        .collect::<Vec<RoaringBitmap>>();

    for bitmap in bitmaps {
        for kmer in bitmap {
            kmer_counts[kmer as usize] += 1;
        }
    }

    info!("kmer counts collected! finding the number of times each count occurs...");

    for (kmer, number_of_ones) in kmer_counts.iter().enumerate() {
        println!("{}\t{}", kmer, number_of_ones);
    }

    println!("========== END INDIVIDUAL KMERS ==========");

    let mut num_ones_to_count = HashMap::new();

    for (_kmer, num_ones) in kmer_counts.into_iter().enumerate() {
        match num_ones_to_count.get_mut(&num_ones) {
            None => {
                num_ones_to_count.insert(num_ones, 1_usize);
            }
            Some(count) => *count += 1,
        }
    }

    for (num_ones, count) in num_ones_to_count {
        println!("{}\t{}", num_ones, count);
    }
}
