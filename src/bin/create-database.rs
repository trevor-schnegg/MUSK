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

/// Creates a run length encoding database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, action)]
    /// Flag that specifies whether or not to use canonical kmers
    canonical: bool,

    #[arg(short, long)]
    /// Level of compression. If not supplied no lossy compression is used.
    /// Otherwise, 1 for minimal, 2 for medium, and 3 for heavy compression
    compression_level: Option<usize>,

    #[arg(short, long, default_value_t = 14)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// The location of the output
    /// If a file, an extension is added
    /// If a directory, the normal extension is the file name
    output_location: String,

    #[arg(short, long, default_value_t = 150)]
    /// Number of queries to sample
    query_number: u64,

    #[arg()]
    /// The ordered file2taxid of the sequences
    ordering_file: String,

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
    let ordering_file_path = Path::new(&args.ordering_file);
    let output_loc_path = Path::new(&args.output_location);
    let ref_dir_path = Path::new(&args.reference_directory);
    let compresssion_level = args.compression_level;

    // If canonical was used for the file2taxid, override command line to use canonical
    // Otherwise, use the argument from the command line
    let canonical = if ordering_file_path
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

    // Create the output file so it errors if an incorrect output file is provided before computation
    let output_file = create_output_file(output_loc_path, "musk.db");

    // Load the file2taxid ordering
    let file2taxid_ordering = load_string2taxid(ordering_file_path);

    info!("creating roaring bitmaps for each group...");

    let bitmaps = file2taxid_ordering
        .par_iter()
        .progress()
        .map(|(files, _taxid)| {
            let file_paths = files
                .split("$")
                .map(|file| ref_dir_path.join(file))
                .collect_vec();

            create_bitmap(file_paths, kmer_len, canonical)
        })
        .collect::<Vec<RoaringBitmap>>();

    info!("roaring bitmaps created! constructing database...");

    let mut database = Database::from(bitmaps, canonical, file2taxid_ordering, kmer_len);

    match compresssion_level {
        None => {}
        Some(compression_level) => {
            if compression_level >= 1 {
                database.lossy_compression(min(compression_level, 3));
            }
        }
    }

    dump_data_to_file(&database, output_file).expect("could not output database to file");

    info!("done!");
}
