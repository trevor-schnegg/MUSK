use clap::Parser;
use indicatif::ParallelProgressIterator;
use musk::explore::connected_components;
use musk::io::load_taxid2files;
use musk::utility::{create_bitmap, get_range};
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::path::Path;
use tracing::{debug, info};

/// Creates a matrix of (hamming) distances between bitmaps
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
    /// the file2taxid file
    file2taxid: String,
}

fn main() {
    tracing_subscriber::fmt::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let file2taxid_path = Path::new(&args.file2taxid);

    info!("loading file2taxid at {}", args.file2taxid);

    let taxid2files = load_taxid2files(file2taxid_path);

    info!("file2taxid loaded! exploring files with the same tax id");

    let (lowest_kmer, highest_kmer) = get_range(args.kmer_length, 0, 0);
    for (taxid, files) in taxid2files {
        if files.len() == 1 {
            println!("{}\t{}", files[0], taxid);
            continue;
        }
        debug!(
            "creating sets for taxid '{}' with {} files...",
            taxid,
            files.len()
        );
        let bitmaps = files
            .par_iter()
            .progress()
            .map(|file| {
                create_bitmap(
                    &*file,
                    args.kmer_length,
                    lowest_kmer,
                    highest_kmer,
                    false,
                    args.canonical,
                )
            })
            .collect::<Vec<RoaringBitmap>>();
        debug!("bitmaps created! performing comparisons...");
        let connected_components = connected_components(bitmaps, 0.9);
        for component in connected_components {
            let mut files_string = String::new();
            for file_index in component {
                if files_string.is_empty() {
                    files_string += &*files[file_index];
                } else {
                    files_string += &*(String::from(",") + &*files[file_index])
                }
            }
            println!("{}\t{}", files_string, taxid);
        }
    }
}
