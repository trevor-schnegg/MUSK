use musk::{kmer_iter::KmerIter, rle::Run, tracing::start_musk_tracing_subscriber};

use itertools::Itertools;
use musk::rle::NaiveRunLengthEncoding;
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

    let underflow = usize::MAX.overflowing_add(1);
    println!("{:?}", underflow);

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

    let mut bits_set = (5_usize..50_usize).into_iter().collect_vec();
    let added_vec: Vec<usize> = vec![55, 58, 60, 64, 65, 66, 80, 120];
    for n in added_vec {
        bits_set.push(n);
    }

    let mut build_rle_1 = NaiveRunLengthEncoding::new();
    for int in &bits_set {
        build_rle_1.push(*int);
    }
    let rle_1 = build_rle_1.to_rle();

    println!(
        "{:?}",
        rle_1
            .get_raw_runs()
            .into_iter()
            .map(|x| Run::from_u16(*x))
            .collect_vec()
    );

    let bits_set_iter = bits_set.into_iter();
    let rle_iter = rle_1.iter();
    for (x, y) in bits_set_iter.zip(rle_iter) {
        assert_eq!(x, y);
    }
}
