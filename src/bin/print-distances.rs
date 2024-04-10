use clap::Parser;
use log::{debug, info};
use musk::io::load_data_from_file;
use std::path::Path;

/// Creates an ordering of files based on distances between bitmaps
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg()]
    /// the distances file
    distances: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let distances_file = Path::new(&args.distances);

    info!("loading distances at {}", args.distances);
    let all_distances = load_data_from_file::<Vec<(Vec<u64>, String, u32)>>(distances_file);
    debug!("length of distances: {}", all_distances.len());

    for (index, (distances, id, _taxid)) in all_distances.into_iter().enumerate() {
        let mut print_string = String::from(&*id);
        print_string += "\t";
        let mut distances_iterator = distances.into_iter();
        if let Some(distance) = distances_iterator.next() {
            print_string += &*distance.to_string();
        }
        while let Some(distance) = distances_iterator.next() {
            print_string += " ";
            print_string += &*distance.to_string();
        }
        println!("{}", print_string);
        if index % 1000 == 0 && index != 0 {
            debug!("done with {} sequences", index);
        }
    }
}
