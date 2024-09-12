use musk::io::{create_output_file, dump_data_to_file, load_data_from_file};
use std::path::Path;

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
    // start_musk_tracing_subscriber();

    // debug!("This should be captured only by stdout");
    // info!("This should be captured only by stdout");
    // warn!("This should be captured only by stderr");
    // error!("This should be captured only by stderr");

    let create = false;
    let file_path = Path::new("./test");

    if create {
        let test: Box<[Box<[u16]>]> =
            vec![vec![0].into_boxed_slice(); 4_usize.pow(14)].into_boxed_slice();
        let file = create_output_file(file_path, "ser");
        dump_data_to_file(&test, file).unwrap();
    } else {
        let file = file_path.with_extension("ser");
        let _vec = load_data_from_file::<Box<[Box<[u16]>]>>(&file);
    }
}
