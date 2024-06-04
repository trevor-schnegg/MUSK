use clap::Parser;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use log::info;
use musk::io::load_string2taxid;
use musk::kmer_iter::KmerIter;
use musk::utility::get_fasta_iterator_of_file;
use rayon::iter::IntoParallelIterator;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use rayon::prelude::*;


fn create_kmer_vec(
    files: String,
    kmer_length: usize,
) -> Vec<usize> {
    let mut set = HashSet::new();
    for file in files.split(",") {
        let mut record_iter = get_fasta_iterator_of_file(Path::new(&file));
        while let Some(Ok(record)) = record_iter.next() {
            if record.seq().len() < kmer_length {
                continue;
            }
            for kmer in KmerIter::from(record.seq(), kmer_length) {
                set.insert(kmer);
            }
        }
    }
    set.into_iter().collect_vec()
}

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
    
    let kmer_vecs = load_string2taxid(files2taxid).into_par_iter().progress().map(|(files, _taxid)| {
        create_kmer_vec(files, args.kmer_length)
    }).collect::<Vec<Vec<usize>>>();

    for kmer_vec in kmer_vecs {
        for kmer in kmer_vec {
            kmer_counts[kmer] += 1;
        }
    }

    info!("kmer counts collected! finding the number of times each count occurs...");

    let mut num_ones_to_count = HashMap::new();

    for (kmer, number_of_ones) in kmer_counts.iter().enumerate() {
        println!("{}\t{}", kmer, number_of_ones);
    }

    println!("========== END INDIVIDUAL KMERS ==========");

    for num_ones in kmer_counts {
        match num_ones_to_count.get_mut(&num_ones) {
            None => {num_ones_to_count.insert(num_ones, 1_usize);},
            Some(count) => *count += 1,
        }
    }

    for (num_ones, count) in num_ones_to_count {
        println!("{}\t{}", num_ones, count);
    }

}
