use clap::Parser;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use log::info;
use musk::{
    io::{dump_data_to_file, load_data_from_file},
    kmer_iter::KmerIter,
    utility::get_fasta_iterator_of_file,
};
use rand::{
    distributions::{Distribution, Uniform},
    thread_rng,
};
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::{collections::HashSet, path::Path};

fn create_bitmap(files: &str, subset: &HashSet<u32>, kmer_length: usize) -> RoaringBitmap {
    let mut bitset = RoaringBitmap::new();
    for file in files.split(",") {
        let mut record_iter = get_fasta_iterator_of_file(Path::new(&file));
        while let Some(Ok(record)) = record_iter.next() {
            if record.seq().len() < kmer_length {
                continue;
            }
            for kmer in KmerIter::from(record.seq(), kmer_length).map(|kmer| kmer as u32) {
                if subset.contains(&kmer) {
                    bitset.insert(kmer);
                }
            }
        }
    }
    bitset
}

/// Creates a sample of k-mers from the matrix
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 14)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg()]
    /// The old directory prefix of the fasta files
    old_directory_prefix: String,

    #[arg()]
    /// The old directory prefix of the fasta files
    new_directory_prefix: String,

    #[arg()]
    /// the ordering file
    ordering: String,

    #[arg()]
    /// Location to output the serialzed bitmaps
    output_file: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let ordering_path = Path::new(&args.ordering);
    let output_path = Path::new(&args.output_file);
    let total_kmers = 4_usize.pow(args.kmer_length as u32);
    let subset_size = 4_usize.pow(8);

    let ordering = load_data_from_file::<Vec<(String, u32)>>(ordering_path)
        .into_iter()
        .map(|(files, taxid)| {
            (
                files.replace(&*args.old_directory_prefix, &*args.new_directory_prefix),
                taxid,
            )
        })
        .collect_vec();

    let mut kmer_subset = HashSet::new();
    let mut rng = thread_rng();
    let distribution = Uniform::new(0_u32, total_kmers as u32);
    while kmer_subset.len() < subset_size {
        let sample = distribution.sample(&mut rng);
        kmer_subset.insert(sample);
    }

    info!("creating roaring bitmaps for each group...");
    let outputs = ordering
        .into_par_iter()
        .progress()
        .map(|(files, taxid)| {
            let bitmap = create_bitmap(&*files, &kmer_subset, args.kmer_length);
            (files, bitmap, taxid)
        })
        .collect::<Vec<(String, RoaringBitmap, u32)>>();
    info!("bitmaps created! outputting to file...");
    dump_data_to_file(
        bincode::serialize(&(kmer_subset, outputs)).unwrap(),
        output_path,
    )
    .unwrap();
}
