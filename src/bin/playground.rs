use musk::kmer_iter::KmerIter;

fn main() {
    let seq = "ATGCTGA".as_bytes();
    let mut seq_iter = KmerIter::from(seq, 4);
    while let Some(kmer) = seq_iter.next() {
        println!("{:016b}", kmer);
    }
    println!("{}", 1_u16);
}
