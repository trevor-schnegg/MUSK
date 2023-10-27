use bio::io::fasta;
use clap::Parser;
use musk::utility::get_fasta_iterator_of_file;
use std::io;
use std::path::Path;

/// Converts a fasta directory to a database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg()]
    /// File path of the human genome assembly
    human_genome_file: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let human_genome = Path::new(&args.human_genome_file);

    let mut writer = fasta::Writer::new(io::stdout());

    let mut record_iter = get_fasta_iterator_of_file(human_genome);
    while let Some(Ok(record)) = record_iter.next() {
        if record.id().starts_with("NC_") {
            writer
                .write_record(&record)
                .expect(&*format!("error printing record with id {}", record.id()))
        }
    }
} // end main
