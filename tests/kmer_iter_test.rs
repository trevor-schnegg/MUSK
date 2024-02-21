use musk::kmer_iter::KmerIter;

#[test]
fn short_iteration() {
    let seq = "ATGCTGA".as_bytes();
    let mut seq_iter = KmerIter::from(seq, 3);
    assert_eq!(seq_iter.next().unwrap(), 0b_00001110);
    assert_eq!(seq_iter.next().unwrap(), 0b_00100100);
    assert_eq!(seq_iter.next().unwrap(), 0b_00001001);
    assert_eq!(seq_iter.next().unwrap(), 0b_00010010);
    assert_eq!(seq_iter.next().unwrap(), 0b_00110100);
}
