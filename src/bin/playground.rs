use std::{collections::HashSet, time::Instant};

use itertools::Itertools;
use musk::{rle::NaiveRunLengthEncoding, tracing::start_musk_tracing_subscriber};
use rand::distributions::{Distribution, Uniform};
use tracing::info;

fn main() {
    start_musk_tracing_subscriber();

    // debug!("This should be captured only by stdout");
    // info!("This should be captured only by stdout");
    // warn!("This should be captured only by stderr");
    // error!("This should be captured only by stderr");

    info!("generating and sorting values");
    let distribution = Uniform::new(0, 1_000_000_000);
    let mut rng = rand::thread_rng();
    let mut value_generator = distribution.sample_iter(&mut rng);
    let mut values = HashSet::new();
    while values.len() < 60_000_000 {
        values.insert(value_generator.next().unwrap());
    }
    let mut values = values.into_iter().collect_vec();
    values.sort();

    info!("inserting into naive RLE");
    let mut naive_rle = NaiveRunLengthEncoding::new();
    for value in values {
        naive_rle.push(value);
    }

    info!("compresseing into standard RLE");
    let rle = naive_rle.to_rle();

    info!("iterating RLE");
    let start = Instant::now();
    for _ in rle.iter() {
        continue;
    }
    let elapsed_time = start.elapsed();

    info!("total time to iterate {:?}", elapsed_time);
}
