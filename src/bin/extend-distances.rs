use clap::Parser;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use musk::io::{create_output_file, dump_data_to_file, load_data_from_file, load_string2taxid};
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::create_bitmap;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::path::Path;
use tracing::info;

/// Computes the lower triangle of a pairwise distance matrix from the input sequences (or sequence groups)
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

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// The location of the output
    /// If a file, an extension is added
    /// If a directory, the normal extension is the file name
    output_location: String,

    #[arg()]
    /// the distances file
    distances: String,

    #[arg()]
    /// the file2taxid file of the new fasta files
    file2taxid: String,

    #[arg()]
    /// Directory with fasta files to extend the distances with
    new_reference_directory: String,

    #[arg()]
    /// Directory with fasta files that the distances were created with
    old_reference_directory: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_musk_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let distances_path = Path::new(&args.distances);
    let file2taxid_path = Path::new(&args.file2taxid);
    let kmer_len = args.kmer_length;
    let new_ref_dir_path = Path::new(&args.new_reference_directory);
    let old_ref_dir_path = Path::new(&args.old_reference_directory);
    let output_loc_path = Path::new(&args.output_location);

    // If the file2taxid is grouped and canonical was used, override command line to use canonical
    // Otherwise, use the argument from the command line
    let canonical = if file2taxid_path
        .file_name()
        .expect("provided file2taxid is not a file")
        .to_str()
        .unwrap()
        .contains(".c.")
    {
        true
    } else {
        args.canonical
    };

    info!("use canonical k-mers: {}", canonical);

    // Create the output file
    let output_file = if canonical {
        create_output_file(output_loc_path, "musk.c.pd")
    } else {
        create_output_file(output_loc_path, "musk.pd")
    };

    info!("loading pairwise distances at {}", args.distances);

    let (old_distances, old_file2taxid) =
        load_data_from_file::<(Vec<Vec<u32>>, Vec<(String, usize)>)>(distances_path);
    let old_file2taxid_len = old_file2taxid.len();

    info!("loading files2taxid at {}", args.file2taxid);

    let new_file2taxid = load_string2taxid(file2taxid_path);

    info!("creating bitmaps for the old file2taxid...");

    let old_bitmaps = old_file2taxid
        .par_iter()
        .progress()
        .map(|(files, _taxid)| {
            let file_paths = files
                .split("$")
                .map(|file| old_ref_dir_path.join(file))
                .collect_vec();

            create_bitmap(file_paths, kmer_len, canonical)
        })
        .collect::<Vec<RoaringBitmap>>();

    info!(
        "{} groups need to be added, creating roaring bitmaps for new file2taxid...",
        new_file2taxid.len()
    );

    let new_bitmaps = new_file2taxid
        .par_iter()
        .progress()
        .map(|(files, _taxid)| {
            let file_paths = files
                .split("$")
                .map(|file| new_ref_dir_path.join(file))
                .collect_vec();

            create_bitmap(file_paths, kmer_len, canonical)
        })
        .collect::<Vec<RoaringBitmap>>();

    info!("roaring bitmaps computed, filling out distance matrix...");

    let all_bitmaps = old_bitmaps
        .into_iter()
        .chain(new_bitmaps.into_iter())
        .collect_vec();

    let new_distances = all_bitmaps
        .par_iter()
        .progress()
        .enumerate()
        .filter_map(|(index_1, bitmap_1)| {
            if index_1 < old_file2taxid_len {
                None
            } else {
                Some(
                    all_bitmaps[..=index_1]
                        .iter()
                        .enumerate()
                        .map(|(index_2, bitmap_2)| {
                            if index_1 == index_2 {
                                0
                            } else {
                                let intersection_size = bitmap_1.intersection_len(bitmap_2);
                                // |A| + |B| - (2 * |A & B|)
                                (bitmap_1.len() + bitmap_2.len() - (2 * intersection_size)) as u32
                            }
                        })
                        .collect::<Vec<u32>>(),
                )
            }
        })
        .collect::<Vec<Vec<u32>>>();

    info!("new distances computed! combining and outputting to file...");

    let all_file2taxid = old_file2taxid
        .into_iter()
        .chain(new_file2taxid.into_iter())
        .collect_vec();

    let all_distances = old_distances
        .into_iter()
        .chain(new_distances.into_iter())
        .collect_vec();

    dump_data_to_file(&(all_distances, all_file2taxid), output_file).unwrap();

    info!("done!");
}
