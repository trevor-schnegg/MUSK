use clap::Parser;
use indicatif::ParallelProgressIterator;
use musk::explore::connected_components;
use musk::io::{create_output_file, load_string2taxid};
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::create_bitmap;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use tracing::{debug, info};

/// Groups an input file2taxid
/// Files with the same taxid are compared and if they are similar enough, they are combined
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

    #[arg(short, long, default_value_t = 0.95)]
    /// The Jaccard similarity required to combine reference sequences
    minimum_similarity: f64,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// The location of the output
    /// If a file, an extension is added
    /// If a directory, the normal extension is the file name
    output_location: String,

    #[arg()]
    /// the file2taxid file
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
    let canonical = args.canonical;
    let file2taxid_path = Path::new(&args.file2taxid);
    let kmer_len = args.kmer_length;
    let output_loc_path = Path::new(&args.output_location);
    let ref_dir_path = Path::new(&args.reference_directory);

    info!("use canonical k-mers: {}", canonical);

    // Create the output file
    let mut output_file = if args.canonical {
        create_output_file(output_loc_path, "musk.g.c.f2t")
    } else {
        create_output_file(output_loc_path, "musk.g.f2t")
    };

    info!("loading file2taxid at {}", args.file2taxid);

    // Create a taxid to files hashmap where each taxid has a list of files with that taxid
    let mut taxid2files = HashMap::new();
    for (file, taxid) in load_string2taxid(file2taxid_path) {
        match taxid2files.get_mut(&taxid) {
            None => {
                taxid2files.insert(taxid, vec![file]);
            }
            Some(files_vec) => files_vec.push(file),
        }
    }

    info!("file2taxid loaded! exploring files with the same tax id");

    for (taxid, files) in taxid2files {
        // If there is only 1 file, no comparisons are needed
        if files.len() == 1 {
            output_file
                .write(format!("{}\t{}\n", files[0], taxid).as_bytes())
                .expect("could not write to output file");
            continue;
        }

        debug!(
            "creating sets for taxid '{}' with {} files...",
            taxid,
            files.len()
        );

        let file_paths = files
            .par_iter()
            .map(|file| ref_dir_path.join(file))
            .collect::<Vec<PathBuf>>();

        // Create a bitmap for each file
        let bitmaps = file_paths
            .into_par_iter()
            .progress()
            .map(|file| create_bitmap(vec![file], kmer_len, canonical))
            .collect::<Vec<RoaringBitmap>>();

        debug!("bitmaps created! performing comparisons...");

        let connected_components = connected_components(bitmaps, args.minimum_similarity);

        for component in connected_components {
            let mut files_string = String::new();

            for file_index in component {
                if files_string.is_empty() {
                    files_string += &*files[file_index];
                } else {
                    files_string += &*(String::from("$") + &*files[file_index])
                }
            }

            output_file
                .write(format!("{}\t{}\n", files_string, taxid).as_bytes())
                .expect("could not write to output file");
        }
    }
}
