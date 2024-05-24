use clap::Parser;
use indicatif::ProgressIterator;
use log::info;
use musk::{
    io::{dump_data_to_file, load_data_from_file},
    rle::{BuildRunLengthEncoding, RunLengthEncoding},
};
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::{
    collections::{HashMap, HashSet},
    path::Path,
};

/// Creates a sample of k-mers from the matrix
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg()]
    /// Location to output the serialzed bitmaps
    subset_bitmaps: String,

    #[arg()]
    /// Location to output the serialzed bitmaps
    output_file: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let output_path = Path::new(&args.output_file);
    let subset_bitmaps_path = Path::new(&args.subset_bitmaps);

    let (subset_kmers, subset_bitmaps) = load_data_from_file::<(
        HashSet<u32>,
        Vec<(String, RoaringBitmap, u32)>,
    )>(subset_bitmaps_path);

    info!("creating rle...");
    let mut rle_map: HashMap<u32, BuildRunLengthEncoding> = HashMap::from_iter(
        subset_kmers
            .iter()
            .map(|kmer| (*kmer, BuildRunLengthEncoding::new())),
    );
    for (index, (_files, bitmap, _taxid)) in subset_bitmaps.iter().progress().enumerate() {
        for kmer in bitmap {
            rle_map.get_mut(&kmer).unwrap().push(index);
        }
    }
    let mut compressed_database = rle_map
        .into_par_iter()
        .map(|(kmer, build_rle)| (kmer, build_rle.to_rle()))
        .collect::<Vec<(u32, RunLengthEncoding)>>();
    compressed_database.sort_by_key(|(kmer, _rle)| *kmer);
    dump_data_to_file(
        bincode::serialize(&compressed_database).unwrap(),
        output_path,
    )
    .unwrap();
}
