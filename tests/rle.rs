use itertools::Itertools;
use musk::rle::NaiveRunLengthEncoding;

#[test]
fn first_is_set() {
    let test_vec = vec![0, 8, 64, 65];

    // Go through the process to create the RLE
    let mut test_naive_rle = NaiveRunLengthEncoding::new();
    test_vec.iter().for_each(|x| {
        test_naive_rle.push(*x);
    });
    let test_rle = test_naive_rle.to_rle();

    assert_eq!(test_vec, test_rle.iter().collect_vec());
}

#[test]
fn first_isnt_set() {
    let test_vec = vec![1, 36, 65];

    // Go through the process to create the RLE
    let mut test_naive_rle = NaiveRunLengthEncoding::new();
    test_vec.iter().for_each(|x| {
        test_naive_rle.push(*x);
    });
    let test_rle = test_naive_rle.to_rle();

    assert_eq!(test_vec, test_rle.iter().collect_vec());
}

#[test]
fn exactly_15_zeros() {
    let test_vec = vec![15, 16, 17, 18, 19];

    // Go through the process to create the RLE
    let mut test_naive_rle = NaiveRunLengthEncoding::new();
    test_vec.iter().for_each(|x| {
        test_naive_rle.push(*x);
    });
    let test_rle = test_naive_rle.to_rle();

    assert_eq!(test_vec, test_rle.iter().collect_vec());
}
