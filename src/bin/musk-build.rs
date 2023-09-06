use bincode::serialize;
use clap::Parser;
use log::{debug, info};
use musk::database::Database;
use musk::io::dump_data_to_file;
use musk::utility::{
    convert_to_uppercase, create_fasta_iterator_from_file, get_fasta_files, reverse_complement,
};
use std::path::Path;

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

    // Create database
    info!("Creating database");
    let mut fasta_files = get_fasta_files(reference_loc).into_iter();
    while let Some(file) = fasta_files.next() {
        debug!("reading file: {}", file);
        let mut record_iter = create_fasta_iterator_from_file(Path::new(&file));
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
    }
    dump_data_to_file(
        serialize(&database).expect("could not serialize database"),
        index_out,
    )
    .expect("could not write database to file");

    info!("Database created!");
} // end main
