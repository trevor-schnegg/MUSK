use clap::Parser;
use musk::big_exp_float::BigExpFloat;
use musk::database::Database;
use musk::io::{create_output_file, load_data_from_file};
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::get_fastq_iter_of_file;
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::ops::Neg;
use std::path::Path;
use std::sync::Mutex;
use tracing::{info, warn};

/// Creates a run length encoding database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 12)]
    /// The exponent e for the significance of hits
    /// Used in the equation 10^{-e} to determine statistical significance
    /// MUST be lower than the cutoff provided for database construction
    exp_cutoff: i32,

    #[arg(short, long, default_value_t = 100)]
    // The maximum number of queries to use in the binomial function
    max_queries: u64,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// Where to write the output.
    /// If a file, '.musk.r2t' is added.
    /// If a directory, 'musk.r2t' will be the file name.
    /// Name means: musk, (r)ead ID (2) (t)ax ID
    output_location: String,

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
    let cutoff_threshold = BigExpFloat::from_f64(10.0_f64.powi((args.exp_cutoff).neg()));
    let database_path = Path::new(&args.database);
    let output_loc_path = Path::new(&args.output_location);
    let reads_path = Path::new(&args.reads);

    let output_file = create_output_file(output_loc_path, "musk.r2t");
    let writer = Mutex::new(BufWriter::new(output_file));

    info!("loading database at {:?}", database_path);
    let database = load_data_from_file::<Database>(database_path);

    info!("classifying reads...");
    let read_iter = get_fastq_iter_of_file(reads_path);
    read_iter
        .par_bridge()
        .into_par_iter()
        .for_each(|record_result| match record_result {
            Err(_) => {
                warn!("error encountered while reading fastq file");
                warn!("skipping the read that caused the error")
            }
            Ok(record) => {
                let classification =
                    database.classify(record.seq(), cutoff_threshold, args.max_queries);
                let mut writer = writer.lock().unwrap();
                match classification {
                    Some((file, taxid)) => {
                        writer
                            .write(format!("{}\t{}\t{}\n", record.id(), file, taxid).as_bytes())
                            .expect("could not write to output file");
                    }
                    None => {
                        writer
                            .write(format!("{}\tNA\t0\n", record.id()).as_bytes())
                            .expect("could not write to output file");
                    }
                };
            }
        });

    info!("done!");
}
