use musk::{
    kmer_iter::KmerIter,
    rle::{Run, RunLengthEncodingIter}, tracing::start_musk_tracing_subscriber,
};

use itertools::Itertools;
use musk::rle::BuildRunLengthEncoding;
use tracing::{debug, error, info, warn};

// const CLEAR_BITS: usize = 2_usize.pow((14 * 2) as u32) - 1;

// fn reverse_compliment(kmer: usize) -> usize {
//     let mut buffer = 0;
//     let mut complement_kmer = (!kmer) & CLEAR_BITS;
//     for _ in 0..14 {
//         // Pop the right-most letter
//         let letter = complement_kmer & 3;
//         complement_kmer >>= 2;
//         // Add to the right of the buffer
//         buffer <<= 2;
//         buffer |= letter;
//     }
//     buffer
// }

fn main() {
    start_musk_tracing_subscriber();

    debug!("This should be captured only by stdout");
    info!("This should be captured only by stdout");
    warn!("This should be captured only by stderr");
    error!("This should be captured only by stderr");

    let seq = "ATGCTGA".as_bytes();
    let mut seq_iter = KmerIter::from(seq, 3, false);
    while let Some(kmer) = seq_iter.next() {
        println!("{:06b}", kmer);
        let (forward_kmer, rev_comp_kmer) = seq_iter.get_curr_kmers();
        println!(
            "kmer: {:06b}, rev comp kmer: {:06b}",
            forward_kmer, rev_comp_kmer
        );
    }

    let _maximum = (1_usize << 14) - 1;

    let mut dense_vector = (5_usize..50_usize).into_iter().collect_vec();
    let added_vec: Vec<usize> = vec![55, 58, 60, 64, 65, 66, 80, 120];
    for n in added_vec {
        dense_vector.push(n);
    }

    let mut build_rle_1 = BuildRunLengthEncoding::new();
    for int in &dense_vector {
        build_rle_1.push(*int);
    }
    let rle_1 = build_rle_1.to_rle();

    println!(
        "{:?}",
        rle_1
            .get_vector()
            .into_iter()
            .map(|x| Run::from_u16(*x))
            .collect_vec()
    );

    let ground_truth = dense_vector.into_iter();
    let test_rle = RunLengthEncodingIter::from_runs_vector(rle_1.get_vector());
    for (x, y) in ground_truth.zip(test_rle) {
        assert_eq!(x, y);
    }
}
