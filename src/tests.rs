use crate::utility::vec_dna_bytes_to_u32;

#[test]
fn dna_to_u32_test() {
    // A = 00, C = 01, G = 10, T = 11
    let kmer = "ATTGGGGGCCTGTCG".as_bytes();
    assert_eq!(
        vec_dna_bytes_to_u32(kmer).unwrap(),
        0b_00_00_11_11_10_10_10_10_10_01_01_11_10_11_01_10
    )
}
