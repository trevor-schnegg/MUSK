use clap::Parser;
use indicatif::ProgressIterator;
use musk::io::load_string2taxid;
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::{get_fasta_files, get_fasta_iter_of_file};
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use taxonomy::{ncbi, TaxRank, Taxonomy};
use tracing::{error, info, warn};

/// Creates a file2taxid file of the form <fasta-file>\t<tax-id> for a given a reference location
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// Name of the output file
    output_file: String,

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
    let output_file_path = Path::new(&args.output_file);
    let reference_dir_path = Path::new(&args.reference_directory);

    // Add extension to the output file
    let mut output_file =
        File::create(output_file_path.join(".musk.f2t")).expect("could not create output file");

    info!("reading accession2taxid at {}", args.accession2taxid);

    let accession2taxid: HashMap<String, usize> =
        HashMap::from_iter(load_string2taxid(accession2taxid_path).into_iter());

    info!(
        "accession2taxid loaded! reading taxonomy at {}",
        args.taxonomy_directory
    );

    let taxonomy = ncbi::load(Path::new(&args.taxonomy_directory)).unwrap();

    info!(
        "taxonomy read! looping through sequences at {}",
        args.reference_directory
    );

    for file in get_fasta_files(reference_dir_path).into_iter().progress() {
        // Get the first record from the fasta file
        let first_record = match get_fasta_iter_of_file(&file).next() {
            None => {
                warn!(
                    "no first record found in fasta file at {:?}. skipping...",
                    file
                );
                continue;
            }
            Some(record_result) => match record_result {
                Ok(record) => record,
                Err(e) => {
                    error!("{:?}", e);
                    warn!(
                        "previous error found while parsing fasta file at {:?}. skipping...",
                        file
                    );
                    continue;
                }
            },
        };

        let mut taxid = accession2taxid[first_record.id()];

        // Try to move the taxid up to the species level, if possible
        match taxonomy.parent_at_rank(taxid as usize, TaxRank::Species) {
            Ok(Some((species_taxid, _distance))) => {
                taxid = species_taxid;
            }
            _ => {}
        }

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
