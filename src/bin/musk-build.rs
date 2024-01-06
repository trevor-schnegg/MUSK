use bincode::serialize;
use clap::Parser;
use log::{debug, info};
use musk::database::Database;
use musk::io::dump_data_to_file;
use musk::utility::{get_fasta_files, get_fasta_iterator_of_file};
use std::path::Path;

/// Converts a fasta directory to a database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 15)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg()]
    /// Location to store the resulting index
    index_out: String,

    #[arg()]
    /// Directory with fasta files to create reference from
    reference_loc: String,

    #[arg()]
    /// Directory containing names.dmp and nodes.dmp
    taxonomy_dir: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let index_out = Path::new(&args.index_out);
    let reference_loc = Path::new(&args.reference_loc);

    // Create database variable
    let mut database = Database::new(args.kmer_length);

    // Create database
    info!("Creating database");
    let mut fasta_files = get_fasta_files(reference_loc).into_iter();
    while let Some(file) = fasta_files.next() {
        debug!("reading file: {}", file);
        let mut record_iter = get_fasta_iterator_of_file(Path::new(&file));
        while let Some(Ok(record)) = record_iter.next() {
            if record.seq().len() < args.kmer_length {
                continue;
            }
            database.insert_record(record);
        }
    }
    dump_data_to_file(
        serialize(&database).expect("could not serialize database"),
        index_out,
    )
    .expect(&*format!("could not write database to {:?}", index_out));

    info!("Database creation complete!");
} // end main
