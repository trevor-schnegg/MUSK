use clap::Parser;
use log::{info};
use probabilitic_classifier::database::Database;
use probabilitic_classifier::utility::{
    convert_to_uppercase, create_fasta_iterator_from_file, reverse_complement,
};
use std::path::Path;
use bincode::serialize;
use probabilitic_classifier::io::{dump_data_to_file};

/// Converts a fasta directory to a database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
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
    let mut database = Database::new(16);

    info!("Creating database");
    // Create database string and get probabilities
    let mut record_iter = create_fasta_iterator_from_file(reference_loc);
    while let Some(Ok(record)) = record_iter.next() {
        if record.seq().len() < 16 {
            continue;
        }
        let uppercase_record_seq = convert_to_uppercase(record.seq());
        let reverse_complement_seq = reverse_complement(&*uppercase_record_seq);
        database.insert_record(
            uppercase_record_seq,
            reverse_complement_seq,
            record.id().to_string(),
        );
    }
    dump_data_to_file(serialize(&database).expect("could not serialize database"), index_out).expect("could not write database to file");

    info!("Database created!");
} // end main