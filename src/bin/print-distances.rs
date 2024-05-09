use clap::Parser;
use itertools::Itertools;
use log::{debug, info};
use musk::io::load_data_from_file;
use roaring::RoaringBitmap;
use std::{collections::{HashMap, HashSet}, path::Path};

/// Creates an ordering of files based on distances between bitmaps
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg()]
    /// the distances file
    subset_kmers: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let subset_bitmaps_path= Path::new(&args.subset_kmers);

    info!("loading kmer subset at {}", args.subset_kmers);
    let (subset_kmers, subset_bitmaps) = load_data_from_file::<(
        HashSet<u32>,
        Vec<(String, RoaringBitmap, u32)>,
    )>(subset_bitmaps_path);

    debug!("{}", subset_bitmaps.len());

    let mut subset_kmers_map: HashMap<u32, Vec<u32>> = HashMap::from_iter(subset_kmers.into_iter().map(|kmer| (kmer, vec![])));

    for (index, (_files, kmers, _taxid)) in subset_bitmaps.into_iter().enumerate() {
        for kmer in kmers {
            subset_kmers_map.get_mut(&kmer).unwrap().push(index as u32);
        }
    }
    for (kmer, bits_set) in subset_kmers_map.into_iter().sorted_by_key(|(kmer, _set_bits)| *kmer) {
        let mut print_string = kmer.to_string();
        print_string += "\t";
        let mut bits_set_iter = bits_set.into_iter();
        if let Some(bit_set) = bits_set_iter.next() {
            print_string += &*bit_set.to_string();
        }
        while let Some(bit_set) = bits_set_iter.next() {
            print_string += " ";
            print_string += &*bit_set.to_string();
        }
        println!("{}", print_string);
    }
}
