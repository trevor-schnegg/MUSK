use musk::kmer_iter::KmerIter;

fn main() {
    let seq = "ATGCTGA".as_bytes();
    let mut seq_iter = KmerIter::from(seq, 3);
    while let Some(kmer) = seq_iter.next() {
        println!("{:08b}", kmer);
    }
}
