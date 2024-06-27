use clap::Parser;
use indicatif::ParallelProgressIterator;
use musk::io::{create_output_file, load_string2taxid};
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::{get_fasta_files, get_fasta_iter_of_file};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use taxonomy::{ncbi, TaxRank, Taxonomy};
use tracing::{error, info, warn};

/// Creates a file2taxid file of the form <fasta-file>\t<tax-id> for a given a reference location
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// The location of the output
    /// If a file, an extension is added
    /// If a directory, the normal extension is the file name
    output_location: String,

    #[arg()]
    /// Accession2taxid file
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
    let mut output_file = create_output_file(output_loc_path, "musk.f2t");

    info!("reading accession2taxid at {}", args.accession2taxid);

    let accession2taxid: HashMap<String, usize> =
        HashMap::from_iter(load_string2taxid(accession2taxid_path).into_iter());

    info!(
        "accession2taxid loaded! reading taxonomy at {}",
        args.taxonomy_directory
    );

    let taxonomy = ncbi::load(taxonomy_dir).unwrap();

    info!(
        "taxonomy read! looping through sequences at {}",
        args.reference_directory
    );

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
                        match taxonomy.parent_at_rank(taxid as usize, TaxRank::Species) {
                            Ok(Some((species_taxid, _distance))) => {
                                taxid = species_taxid;
                            }
                            _ => {}
                        };

                        Some((file, taxid))
                    }
                    Err(e) => {
                        error!("{:?}", e);
                        warn!(
                            "previous error found while parsing fasta file at {:?}. skipping...",
                            file
                        );
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
}
