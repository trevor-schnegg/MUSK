use bio::io::{fasta, fastq};
use clap::Parser;
use musk::io::create_output_file;
use musk::tracing::start_musk_tracing_subscriber;
use musk::utility::{get_fasta_iter_of_file, get_fastq_iter_of_file};
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
    /// The location of the output
    /// If a file, an extension is added
    /// If a directory, the normal extension is the file name
    output_location: String,

    #[arg()]
    /// Fasta reads file to chop
    reads: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_musk_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let output_loc_path = Path::new(&args.output_location);
    let reads_path = Path::new(&args.reads);

    if args.fasta {
        let output_file = create_output_file(output_loc_path, "chopped_reads.fasta");
        let mut writer = fasta::Writer::new(output_file);

        let mut fasta_reads_iter = get_fasta_iter_of_file(reads_path);

        while let Some(Ok(read)) = fasta_reads_iter.next() {
            let seq = if read.seq().len() < 180 {
                read.seq()
            } else {
                &read.seq()[..args.length]
            };
            writer.write(read.id(), read.desc(), seq).unwrap();
        }
    } else {
        let output_file = create_output_file(output_loc_path, "chopped_reads.fastq");
        let mut writer = fastq::Writer::new(output_file);

        let mut fastq_reads_iter = get_fastq_iter_of_file(reads_path);

        while let Some(Ok(read)) = fastq_reads_iter.next() {
            let (seq, qual) = if read.seq().len() < 180 {
                (read.seq(), read.qual())
            } else {
                (&read.seq()[..args.length], &read.qual()[..args.length])
            };
            writer.write(read.id(), read.desc(), seq, qual).unwrap();
        }
    }

    info!("done!");
}
