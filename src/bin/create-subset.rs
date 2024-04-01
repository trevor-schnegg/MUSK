use clap::Parser;
use rand::distributions::{Distribution, Uniform};
use std::{collections::HashMap, path::Path};
use musk::{io::{load_data_from_file, load_string2taxid}, kmer_iter::KmerIter, utility::get_fasta_iterator_of_file};
use log::info;

pub fn push_index(bitset: &mut Vec<u8>, bit_to_set: usize) -> () {
    let (byte, bit) = (bit_to_set / 8, bit_to_set % 8);
    bitset[byte] |= 1 << bit;
}

/// Creates a sample of k-mers from the matrix
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = false)]
    /// If the file2taxid file is actually an ordering file, use this flag
    is_ordering: bool,

    #[arg(short, long, default_value_t = 14)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg()]
    /// the file2taxid file
    file2taxid: String,

    #[arg()]
    /// Location to output the sample
    output_file: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let ordering_path = Path::new(&args.file2taxid);
    let total_kmers = 4_usize.pow(args.kmer_length as u32);

    let ordering = {if args.is_ordering {
        info!("deserializing ordering from {}", args.file2taxid);
        let ordering = load_data_from_file::<Vec<(String, u32)>>(ordering_path);
        info!("ordering deserialzed! creating subset...");
        ordering
    } else {
        info!("loading file2taxid from {}", args.file2taxid);
        let ordering = load_string2taxid(ordering_path).into_iter().collect::<Vec<(String, u32)>>();
        info!("file loaded! creating subset...");
        ordering
    }};

    let mut subset = HashMap::new();
    let subset_size = 4_usize.pow(9);
    let uniform_distribution = Uniform::new(0, total_kmers);
    let mut rng = rand::thread_rng();
    while subset.len() < subset_size {
        let kmer = uniform_distribution.sample(&mut rng);
        match subset.get(&kmer) {
            None => {subset.insert(kmer, vec![0_u8; (ordering.len() / 8) + 1]);},
            Some(_) => {continue;},
        }
    }

    for (index, (file, _)) in ordering.into_iter().enumerate() {
        let mut fasta_iterator = get_fasta_iterator_of_file(&Path::new(&file));
        while let Some(Ok(record)) = fasta_iterator.next() {
            let kmer_iterator = KmerIter::from(record.seq(), args.kmer_length);
            for kmer in kmer_iterator {
                match subset.get_mut(&kmer) {
                    None => {continue;},
                    Some(bit_vector) => {push_index(bit_vector, index)},
                }
            }
        }
    }
}
