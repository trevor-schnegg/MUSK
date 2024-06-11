use clap::Parser;
use musk::io::load_data_from_file;
use musk::rle::RunLengthEncoding;
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::get_fasta_iterator_of_file;
use std::path::Path;
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

    let (database, file2taxid, p_values) =
        load_data_from_file::<(Vec<RunLengthEncoding>, Vec<(String, usize)>, Vec<f64>)>(
            database_path,
        );

    let mut read_iter = get_fasta_iterator_of_file(reads_path);

    while let Some(Ok(read)) = read_iter.next() {
        if args.canonical {
        } else {
        }
    }

    info!("done!");
}
