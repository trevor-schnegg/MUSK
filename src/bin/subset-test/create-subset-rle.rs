use clap::Parser;
use indicatif::ProgressIterator;
use musk::{
    io::{create_output_file, dump_data_to_file, load_data_from_file},
    rle::{NaiveRunLengthEncoding, RunLengthEncoding},
    tracing::start_musk_tracing_subscriber,
};
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::{
    collections::{HashMap, HashSet},
    path::Path,
};
use tracing::info;

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
    output_location: String,
}

fn main() {
    start_musk_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let output_loc_path = Path::new(&args.output_location);
    let subset_bitmaps_path = Path::new(&args.subset_bitmaps);

    let mut output_file = create_output_file(output_loc_path, "musk.subset.rle");

    let (subset_kmers, subset_bitmaps) = load_data_from_file::<(
        HashSet<u32>,
        Vec<(String, RoaringBitmap, u32)>,
    )>(subset_bitmaps_path);

    info!("creating rle...");
    let mut rle_map: HashMap<u32, NaiveRunLengthEncoding> = HashMap::from_iter(
        subset_kmers
            .iter()
            .map(|kmer| (*kmer, NaiveRunLengthEncoding::new())),
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
        &mut output_file,
    )
    .unwrap();
}
