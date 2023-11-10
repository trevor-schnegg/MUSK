use bincode::serialize;
use clap::Parser;
use log::info;
use musk::database::Database;
use musk::io::dump_data_to_file;
use std::path::Path;

/// Converts a fasta file to a database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
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
    let mut database = Database::load(index);
    info!("Database loaded!");

    info!("Beginning test");

    database.update_probabilities_empirically(args.read_length);

    dump_data_to_file(
        serialize(&database).expect("could not serialize database"),
        index,
    )
    .expect(&*format!("could not write database to {:?}", index));

    info!("Done!")
} // end main
