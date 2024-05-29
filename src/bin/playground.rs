use musk::kmer_iter::KmerIter;

use itertools::Itertools;
use musk::rle::{BuildRunLengthEncoding, Run};

fn main() {
    let seq = "ATGCTGA".as_bytes();
    let mut seq_iter = KmerIter::from(seq, 3);
    while let Some(kmer) = seq_iter.next() {
        println!("{:06b}", kmer);
        let (forward_kmer, rev_comp_kmer) = seq_iter.get_curr_kmers();
        println!("kmer: {:06b}, rev comp kmer: {:06b}", forward_kmer, rev_comp_kmer);
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
}
