use clap::Parser;
use itertools::Itertools;
use musk::{
    io::load_data_from_file,
    rle::{Run, RunLengthEncoding},
};
use std::path::Path;

/// Creates a sample of k-mers from the matrix
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg()]
    /// Location of the subset rle vectors
    subset_rle: String,
}

const COMPRESS_MASK: u16 = 0b1000000000000000;
const MAX_ENCODED: u16 = 0b0011111111111111;

const FLIP_TOLERANCE: u16 = 1;

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let subset_rles_path = Path::new(&args.subset_rle);
    let subset_rles = load_data_from_file::<Vec<(u32, RunLengthEncoding)>>(subset_rles_path);
    let subset_vector = subset_rles[0].1.get_vector();

    // surface level statistics - num compressed vs uncompressed blocks compared to total encoded
    let total_elements = subset_vector.len();
    let count_compressed = subset_vector.iter().filter(|&&value| (value >> 15) & 1 == 0).count();
    let count_uncompressed = subset_vector.iter().filter(|&&value| (value >> 15) & 1 == 1).count();
    println!("Total elements: {}", total_elements);
    println!("Count of compressed elements: {}", count_compressed);
    println!("Count of uncompressed elements: {}", count_uncompressed);

    // compute options with cost - focuses on uncompressed blocks composed almost entirely of all-1s or all-0s
    // quick n dirty function to check if a value represents a string of 0s or 1s and the count
    fn decode_neighbor(value: u16) -> Option<(u16, bool)> {
        if (value >> 15) & 1 == 0 {
            let is_one = (value >> 14) & 1 == 1;
            let count = value & 0b0011111111111111;
            return Some((count, is_one));
        }
        None
    }

    let mut output_vector = vec![];
    let mut num_flips = 0;
    let mut uncompressed_bits = 0;
    let mut total_bits = 0;

    let mut i = 0;
    while i < subset_vector.len() {
        // if this is an uncompressed block - attempt to combine it with neighboring blocks
        if subset_vector[i] & COMPRESS_MASK != 0 {
            total_bits += 15;
            uncompressed_bits += 15;
            let current_value = subset_vector[i] & 32767;
            let count_set_bits = current_value.count_ones() as u16;

            let mut merge_prev = false;
            let mut merge_next = false;
            let mut block_count = 15;

            // check left neighbor - as stored in output_vector
            // if compressed, determine if we can combine with it with minimal bit flips
            if i > 0 {
                let prev_value = output_vector[output_vector.len() - 1];
                if let Some((count, is_one)) = decode_neighbor(prev_value) {
                    if count_set_bits <= FLIP_TOLERANCE && !is_one && block_count + count <= MAX_ENCODED { 
                        merge_prev = true;
                        block_count += count;
                    }
                }
            }
            // check right neighbor - still stored in sample_vector
            if i < subset_vector.len() - 1 {
                let next_value = subset_vector[i + 1];
                if let Some((count, is_one)) = decode_neighbor(next_value) {
                    if count_set_bits <= FLIP_TOLERANCE && !is_one && block_count + count <= MAX_ENCODED { 
                        merge_next = true;
                        block_count += count;
                    }
                }
            }

            // combine if possible
            if merge_prev || merge_next {
                if merge_prev {
                    // if previous value is involved, need to update the last value in the output_vector
                    num_flips += count_set_bits;
                    let length = output_vector.len();
                    output_vector[length - 1] = block_count;

                    // if we merged with the next value, skip the next iteration
                    if merge_next {
                        total_bits += subset_vector[i + 1] & 0b0011111111111111;
                        i += 1;
                    }
                } else {
                    total_bits += subset_vector[i + 1] & 0b0011111111111111;
                    num_flips += count_set_bits;
                    output_vector.push(block_count);
                    i += 1;
                }
            } else {
                // we could not combine this block into one of its neighbors - just add it to the output vector
                output_vector.push(subset_vector[i]);
            }
        } else {
            // this block is already compressed - just add it to the output vector
            total_bits += subset_vector[i] & 0b0011111111111111;
            output_vector.push(subset_vector[i]);
        }

        i += 1;
    }

    let count_compressed = output_vector.iter().filter(|&&value| (value >> 15) & 1 == 0).count();
    let count_uncompressed = output_vector.iter().filter(|&&value| (value >> 15) & 1 == 1).count();
    
    println!("\nFinal number of elements: {}", output_vector.len());
    println!("Count of compressed elements: {}", count_compressed);
    println!("Count of uncompressed elements: {}", count_uncompressed);
    println!("Number of bits flipped: {}", num_flips);
   
    println!(
        "Data loss: {} flips from {} uncompressed bits or {} total bits",
        num_flips,
        uncompressed_bits,
        total_bits,
    );

    println!(
        "\nOriginal Compression:\n{:?}",
        subset_rles[0]
            .1
            .get_vector()
            .into_iter()
            .map(|x| Run::from_u16(*x))
            .collect_vec()
    );

    println!(
        "\nNew Compression:\n{:?}",
        output_vector.into_iter()
            .map(|x| Run::from_u16(x))
            .collect_vec()
    );
}
