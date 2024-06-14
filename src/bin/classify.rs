use clap::Parser;
use musk::consts::Consts;
use musk::io::load_data_from_file;
use musk::kmer_iter::KmerIter;
use musk::rle::RunLengthEncoding;
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::get_fasta_iter_of_file;
use musk::{big_exp_float::BigExpFloat, binomial_sf::sf};
use num_traits::One;
use statrs::distribution::{Binomial, DiscreteCDF};
use std::fs::File;
use std::io::Write;
use std::ops::Neg;
use std::path::Path;
use tracing::{info, warn};

/// Creates a run length encoding database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, action)]
    /// Flag that specifies whether or not to use canonical kmers
    /// If canonical kmers are used, each read will only be queried once
    /// Otherwise, the forward and reverse complement will be queried
    canonical: bool,

    #[arg(short, long, default_value_t = 12)]
    cutoff_threshold_exp: i32,

    #[arg(short, long, default_value_t = 14)]
    /// Length of k-mer in the database
    kmer_length: usize,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// Directory to output the read2taxid
    output_directory: String,

    #[arg()]
    /// The database file
    database: String,

    #[arg()]
    /// Directory with fasta reads to query
    reads: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_musk_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let database_path = Path::new(&args.database);
    let output_dir_path = Path::new(&args.output_directory);
    let reads_path = Path::new(&args.reads);

    let mut output_file =
        File::create(output_dir_path.join("readid2taxid")).expect("could not create output file");

    let cutoff_threshold = BigExpFloat::from_f64(10.0_f64.powi((args.cutoff_threshold_exp).neg()));
    let consts = Consts::new();

    let (database, file2taxid, p_values) =
        load_data_from_file::<(Vec<RunLengthEncoding>, Vec<(String, usize)>, Vec<f64>)>(
            database_path,
        );

    let mut read_iter = get_fasta_iter_of_file(reads_path);

    while let Some(Ok(read)) = read_iter.next() {
        let mut collected_hits = vec![0_u64; file2taxid.len()];
        if args.canonical {
            // Collect the hits from the read
            let mut query_count = 0;
            for kmer in KmerIter::from(read.seq(), args.kmer_length, true) {
                query_count += 1;
                for sequence in database[kmer].iter() {
                    collected_hits[sequence] += 1;
                }
            }

            // Classify the hits
            // Would do this using min_by_key but the Ord trait is difficult to implement for float types
            let (mut lowest_prob_index, mut lowest_prob) = (0, BigExpFloat::one());
            for (index, probability) in
                collected_hits
                    .into_iter()
                    .enumerate()
                    .filter_map(|(index, hit_count)| {
                        // Only compute if there was at least 1 hit
                        if hit_count > 0 {
                            let p = p_values[index];
                            let n = query_count;
                            let x = hit_count;
                            // Perform the computation using f64
                            let prob = Binomial::new(p, n).unwrap().sf(x);
                            // If the probability is greater than 0.0, use it
                            let big_exp_float_prob = if prob > 0.0 {
                                BigExpFloat::from_f64(prob)
                            } else {
                                // Otherwise, compute the probability using higher precision
                                sf(p, n, x, &consts)
                            };
                            Some((index, big_exp_float_prob))
                        } else {
                            // If there were 0 hits, don't compute
                            None
                        }
                    })
            {
                // For each index that we computed, compare to find the lowest probability
                // If, for whatever reason, the probabilities are the same, this will use the first one
                if probability < lowest_prob {
                    (lowest_prob_index, lowest_prob) = (index, probability);
                }
            }

            // Print the classification to a file
            let taxid = if lowest_prob < cutoff_threshold {
                file2taxid[lowest_prob_index].1
            } else {
                0
            };
            output_file
                .write(format!("{}\t{}\n", read.id(), taxid).as_bytes())
                .expect("could not write to output file");
        } else {
            warn!("Haven't implemented non-canonical kmers yet");
        }
    }

    info!("done!");
}
