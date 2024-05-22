use clap::Parser;
use log::{debug, info};
use musk::io::load_string2taxid;
use musk::utility::{get_fasta_files, get_fasta_iterator_of_file};
use std::collections::HashMap;
use std::path::Path;
use taxonomy::{ncbi, TaxRank, Taxonomy};

/// Prints to stdout a map in the form of <fasta-file-path>\t<tax-id> given a reference location
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg()]
    /// Accession2taxid file
    accession2taxid: String,

    #[arg()]
    /// Directory with fasta files to create reference from
    reference_location: String,

    #[arg()]
    /// Directory containing names.dmp and nodes.dmp
    taxonomy_directory: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let reference_loc = Path::new(&args.reference_location);

    info!("reading accession2taxid at {}", args.accession2taxid);
    let accession2taxid: HashMap<String, u32> = HashMap::from_iter(load_string2taxid(Path::new(&args.accession2taxid)).into_iter());
    info!(
        "accession2taxid loaded! reading taxonomy at {}",
        args.taxonomy_directory
    );
    let taxonomy = ncbi::load(Path::new(&args.taxonomy_directory)).unwrap();
    info!(
        "taxonomy read! looping through sequences at {}",
        args.reference_location
    );
    let mut fasta_files = get_fasta_files(reference_loc).into_iter().enumerate();
    while let Some((i, file)) = fasta_files.next() {
        if i % 2500 == 0 && i != 0 {
            debug!("done with {} files", i);
        }
        let mut record_iter = get_fasta_iterator_of_file(Path::new(&file));
        let first_record = record_iter.next().unwrap().unwrap();
        let mut taxid = accession2taxid[first_record.id()];
        match taxonomy.parent_at_rank(&*taxid.to_string(), TaxRank::Species) {
            Ok(Some(x)) => {
                taxid = x.0.parse().unwrap();
            }
            _ => {}
        }
        println!("{}\t{}", file, taxid);
    }
}
