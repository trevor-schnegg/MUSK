use clap::Parser;
use indicatif::ParallelProgressIterator;
use log::info;
use musk::io::load_string2taxid;
use musk::utility::{create_bitmap, get_range};
use rayon::iter::IntoParallelIterator;
use std::collections::HashMap;
use std::path::Path;
use std::sync::mpsc;
use rayon::prelude::*;

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

    let (lowest_kmer, highest_kmer) =
        get_range(args.kmer_length, 0, 0);
    
    info!("getting kmer counts...");
    
    let (sender, reciever) = mpsc::channel();
    let _done_with_string2taxid = load_string2taxid(files2taxid).into_par_iter().progress().map_with(sender, |s, x| {
        s.send(create_bitmap(x.0.clone(), args.kmer_length, lowest_kmer, highest_kmer)).unwrap();
        x
    }).collect::<Vec<(String, u32)>>();
    for bitmap in reciever.into_iter() {
        for kmer in bitmap {
            kmer_counts[kmer as usize] += 1;
        }
    }

    info!("kmer counts collected! finding the number of times each count occurs...");

    let mut num_ones_to_count = HashMap::new();

    for count in kmer_counts {
        match num_ones_to_count.get_mut(&count) {
            None => {num_ones_to_count.insert(count, 1);},
            Some(count) => *count += 1,
        }
    }

    for (num_ones, count) in num_ones_to_count {
        println!("{}\t{}", num_ones, count);
    }

}