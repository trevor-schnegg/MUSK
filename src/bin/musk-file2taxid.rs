use clap::Parser;
use indicatif::ParallelProgressIterator;
use musk::io::{create_output_file, load_string2taxid};
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::{get_fasta_files, get_fasta_iter_of_file};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use tracing::{error, info, warn};

/// Creates a file2taxid (.f2t) file of the form <fasta-file>\t<taxid> for a given a reference location
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, verbatim_doc_comment)]
    /// The accession2taxid/seqid2taxid file.
    /// If not provided, all tax ids will be set to 0.
    accession2taxid: Option<String>,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string(), verbatim_doc_comment)]
    /// Where to write the file2taxid (.f2t) file.
    /// If a file is provided, the extention '.musk.f2t' is added.
    /// If a directory is provided, 'musk.f2t' will be the file name.
    output_location: String,

    #[arg()]
    /// Directory with fasta files to create reference from
    reference_directory: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_musk_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let output_loc_path = Path::new(&args.output_location);
    let reference_dir_path = Path::new(&args.reference_directory);

    // Create the output file
    let mut output_file = BufWriter::new(create_output_file(output_loc_path, "musk.f2t"));

    // Get the accession2taxid, if one was provided
    let accession2taxid: Option<HashMap<String, usize>> = match args.accession2taxid {
        None => {
            warn!("no accession2taxid was provided - setting all tax ids to 0");
            warn!("please be sure this is intentional");
            None
        }
        Some(accession2taxid) => {
            let accession2taxid_path = Path::new(&accession2taxid);
            info!("reading accession2taxid at {}", accession2taxid);
            Some(HashMap::from_iter(
                load_string2taxid(accession2taxid_path).into_iter(),
            ))
        }
    };

    info!("searching through files in {}", args.reference_directory);
    let file2taxid = get_fasta_files(reference_dir_path)
        .into_par_iter()
        .progress()
        .filter_map(|file| {
            match &accession2taxid {
                None => Some((file, 0)),
                Some(accession2taxid) => {
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
                                let taxid = *accession2taxid.get(record.id()).expect(&*format!(
                                    "record id {} not in the provided accession2taxid",
                                    record.id(),
                                ));
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
                }
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
