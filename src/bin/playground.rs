use musk::sorted_vector_sets::IntersectIterator;
// use musk::kmer_iter::KmerIter;

fn main() {
    // let seq = "ATGCTGA".as_bytes();
    // let mut seq_iter = KmerIter::from(seq, 3);
    // while let Some(kmer) = seq_iter.next() {
    //     println!("{:08b}", kmer);
    // }

    let vector_1 = vec![1, 2, 3, 6, 13, 15, 18, 22, 40];
    let vector_2 = vec![0, 2, 3, 7, 13, 15, 20, 26, 40];
    let mut intersect_iterator = IntersectIterator::from(&vector_1, &vector_2);
    while let Some(x) = intersect_iterator.next() {
        println!("{}", x);
    }
}
