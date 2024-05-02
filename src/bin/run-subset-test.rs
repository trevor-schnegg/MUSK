use clap::Parser;
use log::{debug, info};
use musk::{
    io::load_data_from_file,
    rle::{BuildRunLengthEncoding, RunLengthEncoding},
};
use rand::{seq::SliceRandom, thread_rng};
use rayon::prelude::*;
use roaring::RoaringBitmap;
use statrs::statistics::Statistics;
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
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let subset_bitmaps_path = Path::new(&args.subset_bitmaps);

    let (subset_kmers, subset_bitmaps) = load_data_from_file::<(
        HashSet<u32>,
        Vec<(String, RoaringBitmap, u32)>,
    )>(subset_bitmaps_path);

    info!("testing rle...");
    let mut rle_map: HashMap<u32, BuildRunLengthEncoding> = HashMap::from_iter(
        subset_kmers
            .iter()
            .map(|kmer| (*kmer, BuildRunLengthEncoding::new())),
    );
    for (index, (_files, bitmap, _taxid)) in subset_bitmaps.iter().enumerate() {
        for kmer in bitmap {
            rle_map.get_mut(&kmer).unwrap().push(index);
        }
    }
    let compressed_database = rle_map
        .into_par_iter()
        .map(|(_kmer, build_rle)| build_rle.to_rle())
        .collect::<Vec<RunLengthEncoding>>();
    let compressed_rle_length = bincode::serialize(&compressed_database).unwrap().len();
    info!("compressed rle length {}", compressed_rle_length);

    info!("testing permutations of rle...");
    let mut rng = thread_rng();
    let mut rle_lengths = vec![];
    for i in 0..100 {
        let mut new_permutation = subset_bitmaps.clone();
        new_permutation.shuffle(&mut rng);
        let mut rle_map: HashMap<u32, BuildRunLengthEncoding> = HashMap::from_iter(
            subset_kmers
                .iter()
                .map(|kmer| (*kmer, BuildRunLengthEncoding::new())),
        );
        for (index, (_files, bitmap, _taxid)) in new_permutation.iter().enumerate() {
            for kmer in bitmap {
                rle_map.get_mut(&kmer).unwrap().push(index);
            }
        }
        let compressed_database = rle_map
            .into_par_iter()
            .map(|(_kmer, build_rle)| build_rle.to_rle())
            .collect::<Vec<RunLengthEncoding>>();
        let compressed_rle_length = bincode::serialize(&compressed_database).unwrap().len();
        rle_lengths.push(compressed_rle_length);
        if i % 10 == 0 && i != 0 {
            debug!("completed {} permutations", i);
        }
    }
    info!(
        "average rle length {}, std dev rle length {}",
        rle_lengths.iter().map(|x| *x as f64).mean(),
        rle_lengths.iter().map(|x| *x as f64).std_dev()
    );
}
