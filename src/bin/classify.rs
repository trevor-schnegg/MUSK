use clap::Parser;
use musk::big_exp_float::BigExpFloat;
use musk::database::Database;
use musk::io::{create_output_file, load_data_from_file};
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::get_fastq_iter_of_file;
use std::io::Write;
use std::ops::Neg;
use std::path::Path;
use std::sync::{mpsc, Arc};
use threadpool::ThreadPool;
use tracing::info;

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

    #[arg(short, long, default_value_t = 18)]
    /// The exponent e for the significance of hits
    /// Used in the equation 10^{-e} to determine statistical significance
    /// MUST be lower than the cutoff provided for database construction
    exp_cutoff: i32,

    #[arg(short, long, default_value_t = 14)]
    /// Length of k-mer in the database
    kmer_length: usize,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// The location of the output
    /// If a file, an extension is added
    /// If a directory, the normal extension is the file name
    output_location: String,

    #[arg(short, long, default_value_t = 12)]
    /// Number of threads to use in classification
    thread_num: usize,

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

    let mut output_file = create_output_file(output_loc_path, "musk.r2t");

    info!("loading database at {:?}", database_path);

    let database = load_data_from_file::<Database>(database_path);

    info!("database loaded! classifying reads...");

    let mut read_iter = get_fastq_iter_of_file(reads_path);

    let (sender, receiver) = mpsc::channel();
    let pool = ThreadPool::new(args.thread_num);
    let database_arc = Arc::new(database);

    while let Some(Ok(read)) = read_iter.next() {
        let sender_clone = sender.clone();
        let database_arc_clone = database_arc.clone();

        pool.execute(move || {
            sender_clone
                .send((
                    read.id().to_string(),
                    database_arc_clone.classify(read.seq(), cutoff_threshold),
                ))
                .unwrap();
        })
    }

    drop(sender);

    for (read, taxid) in receiver {
        // Print the classification to a file
        output_file
            .write(format!("{}\t{}\n", read, taxid).as_bytes())
            .expect("could not write to output file");
    }

    info!("done!");
}
