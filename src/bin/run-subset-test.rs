use clap::Parser;
use itertools::Itertools;
use musk::{
    io::load_data_from_file,
    rle::{Run, RunLengthEncoding},
    tracing::start_musk_tracing_subscriber,
};
use std::path::Path;

/// Creates a sample of k-mers from the matrix
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg()]
    /// Location of the subset rle vectors
    subset_rle: String,
}

fn main() {
    start_musk_tracing_subscriber();

    // Parse arguments from the command line
    let args = Args::parse();
    let subset_rles_path = Path::new(&args.subset_rle);

    let subset_rles = load_data_from_file::<Vec<(u32, RunLengthEncoding)>>(subset_rles_path);

    println!("{:?}", subset_rles[0].1.get_vector());
    println!(
        "{:?}",
        subset_rles[0]
            .1
            .get_vector()
            .into_iter()
            .map(|x| Run::from_u16(*x))
            .collect_vec()
    );
}
