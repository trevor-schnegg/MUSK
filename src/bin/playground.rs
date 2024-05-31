use musk::kmer_iter::KmerIter;

use itertools::Itertools;
use musk::rle::{BuildRunLengthEncoding, Run};

const CLEAR_BITS: usize = 2_usize.pow((14 * 2) as u32) - 1;

fn reverse_compliment(kmer: usize) -> usize {
    let mut buffer = 0;
    let mut complement_kmer = (!kmer) & CLEAR_BITS;
    for _ in 0..14 {
        // Pop the right-most letter
        let letter = complement_kmer & 3;
        complement_kmer >>= 2;
        // Add to the right of the buffer
        buffer <<= 2;
        buffer |= letter;
    }
    buffer
}

fn main() {
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

    let vector = vec![5, 6, 7, 8, 9, 20, 25];

    let mut build_rle_1 = BuildRunLengthEncoding::new();
    for int in vector {
        build_rle_1.push(int);
    }
    let rle_1 = build_rle_1.to_rle();
    println!(
        "{:?}",
        rle_1
            .get_vector()
            .iter()
            .map(|x| Run::from_u16(*x))
            .collect_vec()
    );
    println!(
        "{:?}",
        rle_1
            .get_vector()
            .iter()
            .map(|x| Run::from_u16(*x))
            .collect_vec()
    );

    let mut palendrome_count = 0_usize;
    for i in 0..4_usize.pow(14) {
        if i == reverse_compliment(i) {
            palendrome_count += 1;
        }
    }
    println!("palendromes for 14-mers: {}", palendrome_count);
}
