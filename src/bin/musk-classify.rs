use clap::Parser;
use log::{debug, info};
use musk::io::{load_accession2taxid, load_database};
use musk::utility::{convert_to_uppercase, create_fasta_iterator_from_file};
use std::path::Path;

/// Converts a fasta file to a database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg()]
    /// Accession2taxid
    accession2taxid: String,

    #[arg()]
    /// Name of the file to create a reference from
    index: String,

    #[arg()]
    /// The file containing fasta reads to classify
    reads_file: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let accession2taxid = Path::new(&args.accession2taxid);
    let reads_file = Path::new(&args.reads_file);

    // Get accession2taxid
    info!("Loading accession2taxid from: {}", args.accession2taxid);
    let accession2taxid = load_accession2taxid(accession2taxid);
    info!("accession2taxid loaded!");

    info!("Loading database");
    let database = load_database(Path::new(&args.index));
    info!("Database loaded!");

    info!("Beginning classification");
    let mut read_iter = create_fasta_iterator_from_file(reads_file);
    let mut read_query_count = 0_usize;
    while let Some(Ok(read)) = read_iter.next() {
        let accession = database.query_read(convert_to_uppercase(read.seq()));
        match accession {
            None => {
                println!("{}\t0", read.id())
            }
            Some(accession) => {
                println!(
                    "{}\t{}",
                    read.id(),
                    accession2taxid.get(&accession).unwrap()
                )
            }
        }
        read_query_count += 1;
        if read_query_count % 100000 == 0 {
            debug!("{} reads processed", read_query_count);
        }
    } // end read iterator
    info!("Done!")
} // end main
