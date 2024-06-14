use bio::io::fasta;
use clap::Parser;
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::{get_fasta_iter_of_file, get_fastq_iter_of_file};
use std::fs::File;
use std::path::Path;
use tracing::info;

/// Creates a run length encoding database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, action)]
    /// Specify this flag if the file is a fasta file
    fasta: bool,

    #[arg(short, long, default_value_t = 180)]
    /// Length to chop the read into
    length: usize,

    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// Directory to output the chopped reads
    output_directory: String,

    #[arg()]
    /// Fasta reads file to chop
    reads: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_musk_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let output_dir_path = Path::new(&args.output_directory);
    let reads_path = Path::new(&args.reads);

    let file = File::create(output_dir_path.join("chopped_reads.fasta")).unwrap();
    let mut writer = fasta::Writer::new(file);

    if args.fasta {
        let mut fasta_reads_iter = get_fasta_iter_of_file(reads_path);

        while let Some(Ok(read)) = fasta_reads_iter.next() {
            let seq = if read.seq().len() < 180 {
                read.seq()
            } else {
                &read.seq()[..args.length]
            };
            writer.write(read.id(), None, seq).unwrap();
        }
    } else {
        let mut fastq_reads_iter = get_fastq_iter_of_file(reads_path);

        while let Some(Ok(read)) = fastq_reads_iter.next() {
            let seq = if read.seq().len() < 180 {
                read.seq()
            } else {
                &read.seq()[..args.length]
            };
            writer.write(read.id(), None, seq).unwrap();
        }
    }

    info!("done!");
}
