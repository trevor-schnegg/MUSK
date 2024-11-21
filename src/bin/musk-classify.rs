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

/// Classifies the input reads using a musk database (.db/.cdb) file.
/// Output is a readid2file (.r2f) mapping, including the taxid for the file if it was provided during database construction.
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 6, verbatim_doc_comment)]
    /// The exponent 'e' used in the equation 10^{-e}.
    /// Any calculated p-value below 10^{-e} will result in a classification.
    exp_cutoff: i32,

    #[arg(short, long, default_value_t = 100)]
    // The maximum number of queries to use in the binomial function
    max_queries: u64,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string(), verbatim_doc_comment)]
    /// Where to write the readid2file (.r2f) file.
    /// If a file is provided, the extension '.musk.r2f' is added.
    /// If a directory is provided, 'musk.r2f' will be the file name.
    output_location: String,

    #[arg()]
    /// The database (.db/.cdb) file
    database: String,

    #[arg()]
    /// FASTQ reads file to query
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

    // Create the output file so it errors if an incorrect output file is provided before computation
    let output_file = create_output_file(output_loc_path, "musk.r2f");

    // Create a mutex over a writer to allow multiple threads to write to the output file
    let output_writer = Mutex::new(BufWriter::new(output_file));

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
                let mut writer = output_writer.lock().unwrap();
                match classification {
                    Some((file, taxid)) => {
                        writer
                            .write(format!("{}\t{}\t{}\n", record.id(), file, taxid).as_bytes())
                            .expect("could not write to output file");
                    }
                    None => {
                        writer
                            .write(format!("{}\tU\t0\n", record.id()).as_bytes())
                            .expect("could not write to output file");
                    }
                };
            }
        });

    info!("done!");
}
