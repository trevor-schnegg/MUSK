use clap::Parser;
use log::{debug, info};
use musk::big_exp_float::BigExpFloat;
use musk::database::Database;
use musk::io::load_string2taxid;
use musk::utility::get_fasta_iterator_of_file;
use num_traits::{One, Zero};
use std::ops::Neg;
use std::path::Path;

/// Converts a fasta file to a database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 180)]
    /// Number of k-mers to query for each read
    num_queries: usize,

    #[arg(short, long, default_value_t = 6)]
    /// The positive integer e in the equation: prob = 10^-e. Probabilities < prob are considered matches.
    exp: usize,

    #[arg()]
    /// Accession2taxid file
    accession2taxid: String,

    #[arg()]
    /// The database file
    database: String,

    #[arg()]
    /// The file containing fasta reads to classify
    reads: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let accession2taxid = Path::new(&args.accession2taxid);
    let prob_threshold = BigExpFloat::from_f64(10.0_f64.powi((args.exp as i32).neg()));
    let reads_file = Path::new(&args.reads);

    // Get accession2taxid
    info!("Loading accession2taxid from: {}", args.accession2taxid);
    let accession2taxid = load_string2taxid(accession2taxid);
    info!("accession2taxid loaded!");

    info!("Loading database");
    let database = Database::load(Path::new(&args.database));
    info!("Database loaded!");

    info!("Beginning classification");
    let mut read_iter = get_fasta_iterator_of_file(reads_file);
    let mut read_query_count = 0_usize;
    let mut prob_sum = BigExpFloat::zero();
    let mut num_prob_sum = 0_usize;
    let mut lowest_prob = BigExpFloat::one();
    while let Some(Ok(read)) = read_iter.next() {
        let read_id = read.id().to_string();
        let (prob, accession) = match database.classify_read(read, args.num_queries) {
            None => {
                println!("{}\t0", read_id);
                continue;
            }
            Some(t) => t,
        };
        if prob < lowest_prob {
            lowest_prob = prob;
        }
        prob_sum = prob_sum + prob;
        num_prob_sum += 1;
        if prob < prob_threshold {
            println!("{}\t{}", read_id, accession2taxid.get(accession).unwrap());
        } else {
            println!("{}\t0", read_id);
        }
        read_query_count += 1;
        if read_query_count % 100000 == 0 {
            debug!("{} reads processed", read_query_count);
        }
    } // end read iterator
    debug!(
        "average prob was {:?}",
        prob_sum / BigExpFloat::from_f32(num_prob_sum as f32)
    );
    debug!("lowest observed was {:?}", lowest_prob);
    info!("Done!")
} // end main
