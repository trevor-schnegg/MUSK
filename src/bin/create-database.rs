use clap::Parser;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use musk::database::Database;
use musk::io::{create_output_file, dump_data_to_file, load_string2taxid};
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::create_bitmap;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::cmp::min;
use std::path::Path;
use tracing::info;

const CANONICAL: bool = true;

/// Creates a run length encoding database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long)]
    /// Level of compression. If not supplied no lossy compression is used.
    /// Otherwise, 1 for minimal, 2 for medium, and 3 for heavy compression
    compression_level: Option<usize>,

    #[arg(short, long, default_value_t = 14)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// Where to write the output
    /// If a file, '.musk.db' is added
    /// If a directory, 'musk.db' will be the file name
    /// Name means: musk, (d)ata(b)ase
    output_location: String,

    #[arg()]
    /// The (preferrably ordered) file2taxid map
    file2taxid: String,

    #[arg()]
    /// Directory with fasta files to create reference from
    reference_directory: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_musk_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let kmer_len = args.kmer_length;
    let file2taxid_path = Path::new(&args.file2taxid);
    let output_loc_path = Path::new(&args.output_location);
    let ref_dir_path = Path::new(&args.reference_directory);
    let compresssion_level = args.compression_level;

    // Create the output file so it errors if an incorrect output file is provided before computation
    let output_file = create_output_file(output_loc_path, "musk.db");

    // Load the file2taxid ordering
    info!("loading file2taxid at {}", args.file2taxid);
    let file2taxid_ordering = load_string2taxid(file2taxid_path);

    info!("creating roaring bitmaps for each group...");
    let bitmaps = file2taxid_ordering
        .par_iter()
        .progress()
        .map(|(files, _taxid)| {
            // Split the files up if they are grouped
            let file_paths = files
                .split("$")
                .map(|file| ref_dir_path.join(file))
                .collect_vec();

            create_bitmap(file_paths, kmer_len, CANONICAL)
        })
        .collect::<Vec<RoaringBitmap>>();

    info!("roaring bitmaps created! constructing database...");
    let mut database = Database::from(bitmaps, CANONICAL, file2taxid_ordering, kmer_len);

    if let Some(comp_level) = compresssion_level {
        if comp_level >= 1 {
            database.lossy_compression(min(comp_level, 3));
        }
    };

    dump_data_to_file(&database, output_file).expect("could not output database to file");

    info!("done!");
}
