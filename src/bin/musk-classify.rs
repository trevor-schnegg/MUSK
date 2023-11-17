use clap::Parser;
use log::{debug, info};
use musk::database::Database;
use musk::io::load_accession2taxid;
use musk::utility::get_fasta_iterator_of_file;
use std::ops::Neg;
use std::path::Path;

/// Converts a fasta file to a database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long)]
    /// Number of queries to perform for each read
    num_queries: Option<usize>,

    #[arg(short, long)]
    /// The exponent of the required probability threshold
    exp: Option<usize>,

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
    let true_exponent = match args.exp {
        None => None,
        Some(exp) => Some((exp as i32).neg()),
    };
    let reads_file = Path::new(&args.reads_file);

    // Get accession2taxid
    info!("Loading accession2taxid from: {}", args.accession2taxid);
    let accession2taxid = load_accession2taxid(accession2taxid);
    info!("accession2taxid loaded!");

    info!("Loading database");
    let database = Database::load(Path::new(&args.index));
    info!("Database loaded!");

    info!("Beginning classification");
    let mut read_iter = get_fasta_iterator_of_file(reads_file);
    let mut read_query_count = 0_usize;
    let mut prob_sum = 0.0;
    let mut num_prob_sum = 0_usize;
    let mut lowest_prob = 0.0;
    while let Some(Ok(read)) = read_iter.next() {
        let read_id = read.id().to_string();
        let (accession, prob) = database.classify_read(read, args.num_queries, true_exponent);
        if prob < lowest_prob {
            lowest_prob = prob;
        }
        prob_sum += prob;
        num_prob_sum += 1;
        match accession {
            None => {
                println!("{}\t0", read_id);
            }
            Some(accession) => {
                println!("{}\t{}", read_id, accession2taxid.get(accession).unwrap());
            }
        }
        read_query_count += 1;
        if read_query_count % 100000 == 0 {
            debug!("{} reads processed", read_query_count);
        }
    } // end read iterator
    debug!("average prob was {}", prob_sum / num_prob_sum as f64);
    debug!("lowest observed was {}", lowest_prob);
    info!("Done!")
} // end main
