use clap::Parser;
use musk::database::Database;
use musk::io::load_data_from_file;
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::get_fasta_iter_of_file;
use std::fs::File;
use std::io::Write;
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

    #[arg(short, long, default_value_t = 14)]
    /// Length of k-mer in the database
    kmer_length: usize,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// Directory to output the read2taxid
    output_directory: String,

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
    let database_path = Path::new(&args.database);
    let output_dir_path = Path::new(&args.output_directory);
    let reads_path = Path::new(&args.reads);

    let mut output_file =
        File::create(output_dir_path.join("readid2taxid")).expect("could not create output file");

    info!("loading database at {:?}", database_path);

    let database = load_data_from_file::<Database>(database_path);

    info!("database loaded! classifying reads...");

    let mut read_iter = get_fasta_iter_of_file(reads_path);

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
                    database_arc_clone.classify(read.seq()),
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
