use clap::Parser;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use musk::{
    io::{create_output_file, dump_data_to_file, load_string2taxid},
    kmer_iter::KmerIter,
    tracing::start_musk_tracing_subscriber,
    utility::get_fasta_iter_of_file,
};
use rand::{
    distributions::{Distribution, Uniform},
    thread_rng,
};
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::{
    collections::HashSet,
    path::{Path, PathBuf},
};
use tracing::info;

fn create_bitmap(
    files: Vec<PathBuf>,
    subset: &HashSet<u32>,
    kmer_length: usize,
    canonical: bool,
) -> RoaringBitmap {
    let mut bitset = RoaringBitmap::new();
    for file in files {
        let mut record_iter = get_fasta_iter_of_file(&file);
        while let Some(Ok(record)) = record_iter.next() {
            if record.seq().len() < kmer_length {
                continue;
            }
            for kmer in KmerIter::from(record.seq(), kmer_length, canonical).map(|kmer| kmer as u32)
            {
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
    #[arg(short, long, action)]
    /// Flag that specifies whether or not to use canonical kmers
    canonical: bool,

    #[arg(short, long, default_value_t = 14)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg()]
    /// the ordering file
    ordering: String,

    #[arg()]
    /// Location to output the serialzed bitmaps
    output_location: String,

    #[arg()]
    /// Directory with fasta files to create reference from
    reference_location: String,
}

const SUBSET_SIZE: usize = 4_usize.pow(9);

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_musk_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let output_loc_path = Path::new(&args.output_location);
    let ordering_path = Path::new(&args.ordering);
    let reference_path = Path::new(&args.reference_location);

    let total_kmers = 4_usize.pow(args.kmer_length as u32);

    let mut output_file = create_output_file(output_loc_path, "musk.subset.bitmaps");

    info!("loading ordering at {:?}", ordering_path);

    let ordering = load_string2taxid(ordering_path);

    info!(
        "ordering loaded! creating a random sample of k-mers of size {} out of {} total k-mers",
        SUBSET_SIZE, total_kmers
    );

    let mut kmer_subset = HashSet::new();
    let mut rng = thread_rng();
    let distribution = Uniform::new(0_u32, total_kmers as u32);
    while kmer_subset.len() < SUBSET_SIZE {
        let sample = distribution.sample(&mut rng);
        kmer_subset.insert(sample);
    }

    info!("sample created! creating roaring bitmaps for each group...");

    let outputs = ordering
        .par_iter()
        .progress()
        .map(|(files, _taxid)| {
            let file_paths = files
                .split("$")
                .map(|file| reference_path.join(file))
                .collect_vec();

            create_bitmap(file_paths, &kmer_subset, args.kmer_length, args.canonical)
        })
        .collect::<Vec<RoaringBitmap>>();

    info!("bitmaps created! outputting to file...");

    dump_data_to_file(
        bincode::serialize(&(kmer_subset, outputs)).unwrap(),
        &mut output_file,
    )
    .unwrap();
}
