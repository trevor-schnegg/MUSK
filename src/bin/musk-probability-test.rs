use bincode::serialize;
use bio::io::fasta::Record;
use clap::Parser;
use log::info;
use musk::database::Database;
use musk::generator::create_random_read;
use musk::io::dump_data_to_file;
use std::path::Path;

/// Converts a fasta file to a database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 1000)]
    /// The length of random reads to generate
    num_reads: usize,

    #[arg(short, long, default_value_t = 7000)]
    /// The length of random reads to generate
    read_length: usize,

    #[arg()]
    /// Name of the file to create a reference from
    index: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let index = Path::new(&args.index);

    info!("Loading database");
    let database = Database::load(index);
    info!("Database loaded!");

    info!("Beginning test");

    for record in (0..args.num_reads).map(|i| {
        Record::with_attrs(
            &*i.to_string(),
            None,
            create_random_read(args.read_length).as_bytes(),
        )
    }) {
        database.classify_read(record, Some(args.read_length), None);
    }

    dump_data_to_file(
        serialize(&database).expect("could not serialize database"),
        index,
    )
    .expect(&*format!("could not write database to {:?}", index));

    info!("Done!")
} // end main
