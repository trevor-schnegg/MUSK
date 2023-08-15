use clap::Parser;
use log::{debug, info};
use probabilitic_classifier::database::Database;
use probabilitic_classifier::taxonomy::get_accession_to_tax_id;
use probabilitic_classifier::utility::{
    convert_to_uppercase, create_fasta_iterator_from_file, reverse_complement,
};
use std::path::Path;
use threadpool::ThreadPool;
use std::sync::{Arc, mpsc};

/// Converts a fasta file to a database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg()]
    /// Name of the file to create a reference from
    reference_file: String,

    #[arg()]
    /// Name of the file to create a reference from
    reads_file: String,

    #[arg()]
    /// Directory containing names.dmp, nodes.dmp, nucl_extra.accession2taxid,
    /// nucl_gb.accession2taxid, and nucl_wgs.accession2taxid
    taxonomy_dir: String,

    #[arg()]
    /// k-mers to test
    kmer_len: usize,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let taxonomy_dir = Path::new(&args.taxonomy_dir);
    let reference_file = Path::new(&args.reference_file);
    let reads_file = Path::new(&args.reads_file);

    // Get accession2taxid
    info!(
        "Loading accession2taxid from taxonomy directory: {}",
        args.taxonomy_dir
    );
    let accession2taxid = get_accession_to_tax_id(taxonomy_dir, reference_file);
    info!("accession2taxid loaded!");
    let mut database = Database::new(args.kmer_len);

    info!("Creating database");
    // Create database string and get probabilities
    let mut record_iter = create_fasta_iterator_from_file(reference_file);
    while let Some(Ok(record)) = record_iter.next() {
        if record.seq().len() < args.kmer_len {
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
    let database = Arc::new(database);
    info!("Database created!");

    info!("Beginning classification");
    let mut read_iter = create_fasta_iterator_from_file(reads_file);
    let mut read_query_count = 0_usize;
    let pool = ThreadPool::new(15);
    let (tx, rx) = mpsc::channel();
    while let Some(Ok(read)) = read_iter.next() {
        let tx = tx.clone();
        let database = Arc::clone(&database);
        pool.execute(move || tx.send((read.id().to_string(), database.query_read(convert_to_uppercase(read.seq())))).expect("failed to pass a message to receiver"));
    } // end read iterator
    drop(tx); // drop the original transmitter so that the stream closes when finished
    for (read_id, accession) in rx {
        match accession {
            None => {println!("{}\t0", read_id)}
            Some(accession) => {println!("{}\t{}", read_id, accession2taxid.get(&accession).unwrap())}
        }
        read_query_count += 1;
        if read_query_count % 100000 == 0 {
            debug!("{} reads processed", read_query_count);
        }
    }

    info!("Done!")
} // end main
