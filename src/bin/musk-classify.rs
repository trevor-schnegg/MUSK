use clap::Parser;
use log::{debug, info};
use musk::io::{load_accession2taxid, load_database};
use musk::utility::{convert_to_uppercase, get_fasta_iterator_of_file};
use std::fs::File;
use std::io::Write;
use std::ops::Neg;
use std::path::Path;

/// Converts a fasta file to a database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 400)]
    /// Number of queries to perform for each read
    num_queries: usize,

    #[arg(short, long, default_value_t = 25)]
    /// The exponent of the required probability threshold
    exp: usize,

    #[arg()]
    /// Accession2taxid
    accession2taxid: String,

    #[arg()]
    /// Name of the file to create a reference from
    index: String,

    #[arg()]
    /// Name of the output file
    output_file: String,

    #[arg()]
    /// The file containing fasta reads to classify
    reads_file: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let accession2taxid = Path::new(&args.accession2taxid);
    let true_exponent = (args.exp as i32).neg();
    let output_file = Path::new(&args.output_file);
    let reads_file = Path::new(&args.reads_file);

    // Get accession2taxid
    info!("Loading accession2taxid from: {}", args.accession2taxid);
    let accession2taxid = load_accession2taxid(accession2taxid);
    info!("accession2taxid loaded!");

    info!("Loading database");
    let database = load_database(Path::new(&args.index));
    info!("Database loaded!");

    let mut output_file = File::create(output_file)
        .expect(&*format!("Could not create output file {:?}", output_file));
    info!("Beginning classification");
    let mut read_iter = get_fasta_iterator_of_file(reads_file);
    let mut read_query_count = 0_usize;
    while let Some(Ok(read)) = read_iter.next() {
        let accession = database.query_read(
            convert_to_uppercase(read.seq()),
            args.num_queries,
            true_exponent,
        );
        match accession {
            None => {
                let string = format!("{}\t0\n", read.id());
                output_file
                    .write(string.as_bytes())
                    .expect("failed to write to output file");
            }
            Some(accession) => {
                println!("{}", read.id());
                let string = format!(
                    "{}\t{}\n",
                    read.id(),
                    accession2taxid.get(accession).unwrap()
                );
                output_file
                    .write(string.as_bytes())
                    .expect("failed to write to output file");
            }
        }
        read_query_count += 1;
        if read_query_count % 100000 == 0 {
            debug!("{} reads processed", read_query_count);
        }
    } // end read iterator
    info!("Done!")
} // end main
