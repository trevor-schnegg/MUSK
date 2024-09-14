use clap::Parser;
use indicatif::ParallelProgressIterator;
use musk::io::{create_output_file, load_string2taxid};
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::{get_fasta_files, get_fasta_iter_of_file};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use taxonomy::{ncbi, TaxRank, Taxonomy};
use tracing::{error, info, warn};

/// Creates a file2taxid file of the form <fasta-file>\t<taxid> for a given a reference location
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// Where to write the output
    /// If a file, '.musk.f2t' is added
    /// If a directory, 'musk.f2t' will be the file name
    /// Name means: musk, (f)ile(2)(t)axid
    output_location: String,

    #[arg()]
    /// The accession2taxid map file
    accession2taxid: String,

    #[arg()]
    /// Directory with fasta files to create reference from
    reference_directory: String,

    #[arg()]
    /// Directory containing names.dmp and nodes.dmp
    taxonomy_directory: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_musk_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let accession2taxid_path = Path::new(&args.accession2taxid);
    let output_loc_path = Path::new(&args.output_location);
    let reference_dir_path = Path::new(&args.reference_directory);
    let taxonomy_dir = Path::new(&args.taxonomy_directory);

    // Create the output file
    let mut output_file = BufWriter::new(create_output_file(output_loc_path, "musk.f2t"));

    info!("reading accession2taxid at {}", args.accession2taxid);
    let accession2taxid: HashMap<String, usize> =
        HashMap::from_iter(load_string2taxid(accession2taxid_path).into_iter());

    info!("reading taxonomy at {}", args.taxonomy_directory);
    let taxonomy = ncbi::load(taxonomy_dir).unwrap();

    info!("looping through files at {}", args.reference_directory);
    let file2taxid = get_fasta_files(reference_dir_path)
        .into_par_iter()
        .progress()
        .filter_map(|file| {
            // Get the first record from the fasta file
            match get_fasta_iter_of_file(&file).next() {
                None => {
                    warn!(
                        "no first record found in fasta file at {:?}. skipping...",
                        file
                    );
                    None
                }
                Some(record_result) => match record_result {
                    Ok(record) => {
                        let mut taxid = accession2taxid[record.id()];

                        // Try to move the taxid up to the species level, if possible
                        if let Ok(Some((species_taxid, _distance))) =
                            taxonomy.parent_at_rank(taxid as usize, TaxRank::Species)
                        {
                            taxid = species_taxid;
                        };

                        Some((file, taxid))
                    }
                    Err(e) => {
                        error!("error encountered while parsing fasta file {:?}", file);
                        error!("{:?}", e);
                        warn!("skipping...");
                        None
                    }
                },
            }
        })
        .collect::<Vec<(PathBuf, usize)>>();

    for (file, taxid) in file2taxid {
        // Write the result to the output file
        output_file
            .write(
                format!(
                    "{}\t{}\n",
                    file.file_name().unwrap().to_str().unwrap(),
                    taxid
                )
                .as_bytes(),
            )
            .expect("could not write to output file");
    }

    output_file.flush().unwrap();
}
