use clap::Parser;
use musk::database::Database;
use musk::io::{create_output_file, dump_data_to_file, load_data_from_file};
use musk::tracing::start_musk_tracing_subscriber;
use std::path::Path;
use tracing::info;

/// Creates a run length encoding database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = std::env::current_dir().unwrap().to_str().unwrap().to_string())]
    /// Where to write the output
    /// If a file, '.musk.db' is added
    /// If a directory, 'musk.db' will be the file name
    /// Name means: musk, (d)ata(b)ase
    output_location: String,

    #[arg()]
    /// Level of compression: one of [1, 2, 3]
    compression_level: usize,

    #[arg()]
    /// The uncompressed database file
    database: String,
}

fn main() {
    // Initialize the tracing subscriber to handle debug, info, warn, and error macro calls
    start_musk_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let database_path = Path::new(&args.database);
    let output_loc_path = Path::new(&args.output_location);

    let compression_level = if args.compression_level >= 1 && args.compression_level <= 3 {
        args.compression_level
    } else {
        panic!("compression level was not 1, 2, or 3!")
    };

    let output_file = create_output_file(output_loc_path, "musk.db");

    info!("loading database at {:?}", database_path);
    let mut database = load_data_from_file::<Database>(database_path);

    info!(
        "compressing database using compression level: {}",
        compression_level
    );
    database.lossy_compression(compression_level);

    info!("dumping to file...");
    dump_data_to_file(&database, output_file).expect("could not output database to file");

    info!("done!");
}
