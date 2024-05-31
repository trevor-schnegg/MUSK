use clap::Parser;
use indicatif::ParallelProgressIterator;
use log::info;
use musk::io::load_string2taxid;
use musk::utility::{create_bitmap, get_range};
use rayon::iter::IntoParallelIterator;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::collections::HashMap;
use std::path::Path;

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
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let files2taxid = Path::new(&args.files2taxid);

    let mut kmer_counts = vec![0_usize; 4_usize.pow(args.kmer_length as u32)];

    info!("getting kmer counts...");

    let (lowest_kmer, highest_kmer) = get_range(args.kmer_length, 12, 667);

    let bitmaps = load_string2taxid(files2taxid)
        .into_par_iter()
        .progress()
        .map(|(files, _taxid)| create_bitmap(files, args.kmer_length, lowest_kmer, highest_kmer, false, false))
        .collect::<Vec<RoaringBitmap>>();

    for bitmap in bitmaps {
        for kmer in bitmap {
            kmer_counts[kmer as usize] += 1;
        }
    }

    info!("kmer counts collected! finding the number of times each count occurs...");

    for (kmer, number_of_ones) in kmer_counts.iter().enumerate() {
        if lowest_kmer <= kmer && kmer < highest_kmer {
            println!("{}\t{}", kmer, number_of_ones);
        }
    }

    println!("========== END INDIVIDUAL KMERS ==========");

    let mut num_ones_to_count = HashMap::new();

    for (kmer, num_ones) in kmer_counts.into_iter().enumerate() {
        if lowest_kmer <= kmer && kmer < highest_kmer {
            match num_ones_to_count.get_mut(&num_ones) {
                None => {
                    num_ones_to_count.insert(num_ones, 1_usize);
                }
                Some(count) => *count += 1,
            }
        }
    }

    for (num_ones, count) in num_ones_to_count {
        println!("{}\t{}", num_ones, count);
    }
}
