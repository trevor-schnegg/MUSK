use clap::Parser;
// use itertools::Itertools;
use musk::{
    io::load_data_from_file,
    rle::{RunLengthEncoding},
    // rle::{Run, RunLengthEncoding},
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

    #[clap(long)]
    // invoked with --comp-freq
    comp_freq: bool,
}

const FLIP_TOLERANCE: u16 = 1;
const MAX_ENCODED: u16 = (1 << 14) - 1;

/* 
 * Method computes the number of uncompressed blocks, organized by the number of set bits x the number of 
 * neighboring blocks containing a compressed run of zeros. Outputs the result with the format:
 *
 *                              # of neighboring 0 blocks
 *          set bits    0 neighbors     1 neighbor      2 neighbors
 *              1           ...             ...             ...
 *              2           ...             ...             ...
*/
fn comp_freq_uncompressed(subset_rles_path: &Path) {
    // load data from file and create mutable matrix to store result
    let subset_rles = load_data_from_file::<Vec<(u32, RunLengthEncoding)>>(subset_rles_path);
    let mut freq_matrix: Vec<Vec<u32>> = vec![vec![0; 3]; 14];

    // iterate through each sample in subset_rles
    for subset in &subset_rles {
        let subset_vector = subset.1.get_vector();

        // iterate through each element in the subset_vector
        let mut i = 0;
        while i < subset_vector.len() {
            // determine if subset_vector[i] is an uncompressed element
            let current_element = subset_vector[i];
            let is_uncompressed = current_element & (1 << 15) > 0;
            if is_uncompressed {
                // if we are looking at an uncompressed element, determine its value and # of set bits
                let current_value = current_element & 0b_0111_1111_1111_1111;
                let count_ones = current_value.count_ones() as usize;

                // determine number of neighboring compressed zero blocks and update result
                let mut num_neighbors = 0;
                if i > 0 {
                    let prev_element = subset_vector[i-1];
                    if prev_element & 0b1100_0000_0000_0000 == 0 {
                        num_neighbors += 1;
                    }
                }
                if i < subset_vector.len() - 1 {
                    let next_element = subset_vector[i+1];
                    if next_element & 0b1100_0000_0000_0000 == 0 {
                        num_neighbors += 1;
                    }
                }
                freq_matrix[count_ones - 1][num_neighbors] += 1;
            }

            i += 1;
        }
    }

    // output result
    for row in &freq_matrix {
        println!("{:?}", row);
    }
}

/* 
 * Method performs basic error testing by ensuring that the same number of bits are being represented
 * between the input and compressed vector.
 * 
 * Returns true if the number of bits across both vectors is equal, otherwise false.
 */
fn verify_bits(input: Vec<u16>, compressed: Vec<u16>) -> bool {
    let mut input_bits: u32 = 0;
    let mut result_bits: u32 = 0;

    for block in input {
        if block & (1 << 15) > 0 {
            input_bits += 15;
        } else {
            input_bits += (block & 0b0011_1111_1111_1111) as u32;
        }
    }

    for block in compressed {
        if block & (1 << 15) > 0 {
            result_bits += 15;
        } else {
            result_bits += (block & 0b0011_1111_1111_1111) as u32;
        }
    }

    input_bits == result_bits
}

/*
 * Method identifies uncompressed blocks where the number of set bits is within the range specified by 
 * FLIP_TOLERANCE. If the uncompressed block has neighboring compressed zero block(s), it attempts to flip
 * the bit and combine the uncompressed block into the compressed one(s).
 *
 * In order to combine, the following conditions must be met:
 *      - There can be, at most, FLIP_TOLERANCE many set bits in the uncompressed block
 *      - At least one of the neighboring blocks (at position i - 1 or i + 1) must be a compressed block of zeros.
 *      - The resultant number of zeros across all combined blocks must be <= MAX_ENCODED
 *
 * Method outputs the following statistics:
 *      - Input size (total number of blocks, number of compressed, number of uncompressed)
 *      - Final size (total number of blocks, number of compressed, number of uncompressed)
 *      - Data loss (number of bits flipped vs. number of uncompressed bits vs. total number of bits)
 *      - Summary statistics (percent improved compression, percent data loss)
 */
fn attempt_compression(subset_rles_path: &Path) {
    // load data from file and initialize variables to track statistics
    let subset_rles = load_data_from_file::<Vec<(u32, RunLengthEncoding)>>(subset_rles_path);
    let mut input_num_blocks = 0;
    let mut input_num_compressed = 0;
    let mut input_num_uncompressed = 0;
    
    let mut final_num_blocks = 0;
    let mut final_num_compressed = 0;
    let mut final_num_uncompressed = 0;
    
    let mut req_bit_flips: u32 = 0;
    let mut uncompressed_bits: u32 = 0;
    let mut total_bits: u32 = 0;

    // method checks if a neighboring value is compressed, whether it represents a run of 0s or 1s, and the count
    fn decode_neighbor(value: u16) -> Option<(u16, bool)> {
        if value & (1 << 15) == 0 {
            let is_zero = value & (1 << 14) == 0;
            let count = value & 0b0011_1111_1111_1111;
            return Some((count, is_zero));
        }
        None
    }

    // iterate through each sample in the subset
    for subset in &subset_rles {
        let subset_vector = subset.1.get_vector();

        // update statistics for input file
        input_num_blocks += subset_vector.len();
        input_num_compressed += subset_vector.iter().filter(|&&value| value & (1 << 15) == 0).count();
        input_num_uncompressed += subset_vector.iter().filter(|&&value| value & (1 << 15) > 0).count();

        // initialize vector to hold result for this subset
        let mut result_vector = vec![];

        // iterate through each element in the subset vector
        let mut i = 0;
        while i < subset_vector.len() {
            let current_element = subset_vector[i];

            // determine whether we are looking at an uncompressed block
            let is_uncompressed = current_element & (1 << 15) > 0;
            if is_uncompressed {
                let current_value = current_element & 0b0111_1111_1111_1111;
                let count_ones = current_value.count_ones() as u16;
                uncompressed_bits += 15;
                total_bits += 15;
                
                let mut can_merge_prev = false;
                let mut can_merge_next = false;
                let mut block_count = 15;

                if i > 0 {
                    // retrieve the previous element - stored in the result_vector as the final element
                    let prev_value = result_vector[result_vector.len() - 1];
                    if let Some((count, is_zero)) = decode_neighbor(prev_value) {
                        if count_ones <= FLIP_TOLERANCE && is_zero && count + block_count <= MAX_ENCODED {
                            can_merge_prev = true;
                            block_count += count;
                        }
                    }
                }

                if i < subset_vector.len() - 1 {
                    // retrieve the next element - stored at i + 1 in subset_vector
                    let next_value = subset_vector[i + 1];
                    if let Some((count, is_zero)) = decode_neighbor(next_value) {
                        if count_ones <= FLIP_TOLERANCE && is_zero && count + block_count <= MAX_ENCODED {
                            can_merge_next = true;
                            block_count += count;
                        }
                    }
                }

                // attempt to merge
                if can_merge_prev || can_merge_next {
                    req_bit_flips += count_ones as u32;
                    if can_merge_prev {
                        // need to update the value already saved in result_vector
                        let length = result_vector.len();
                        result_vector[length - 1] = block_count;

                        // if we are merging both previous and next - need to move forward one
                        if can_merge_next {
                            total_bits += (subset_vector[i + 1] & 0b0011_1111_1111_1111) as u32;
                            i += 1;
                        }
                    } else {
                        // only merging with the next block - append to result vector and move forward one
                        total_bits += (subset_vector[i + 1] & 0b0011_1111_1111_1111) as u32;
                        result_vector.push(block_count);
                        i += 1;
                    }
                } else {
                    // could not merge, add it to the result_vector
                    result_vector.push(current_element);
                }
            } else {
                // value is already compressed, add it to the result_vector
                result_vector.push(current_element);
                total_bits += (current_element & 0b0011_1111_1111_1111) as u32;
            }

            i += 1;
        }

        final_num_blocks += result_vector.len();
        final_num_compressed += result_vector.iter().filter(|&&value| value & (1 << 15) == 0).count();
        final_num_uncompressed += result_vector.iter().filter(|&&value| value & (1 << 15) > 0).count();

        // basic error checking
        // subset_vector is a borrowed slice, needs to be converted to an owned vector
        if !verify_bits(subset_vector.to_vec(), result_vector) {
            println!("Error! - number of bits does not match");
        }
    }

    // output result
    println!("\nTotal number of elements: {}", input_num_blocks);
    println!("Count of compressed elements: {}", input_num_compressed);
    println!("Count of uncompressed elements: {}", input_num_uncompressed);

    println!("\nFinal number of elements: {}", final_num_blocks);
    println!("Final count of compressed elements: {}", final_num_compressed);
    println!("Final count of uncompressed elements: {}", final_num_uncompressed);
   
    println!(
        "Data loss: {} flips from {} uncompressed bits or {} total bits",
        req_bit_flips,
        uncompressed_bits,
        total_bits,
    );
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let subset_rles_path = Path::new(&args.subset_rle);

    // if flag set, compute # of uncompressed blocks (# of bits set x # of neighboring 0 blocks)
    if args.comp_freq {
        comp_freq_uncompressed(subset_rles_path);
    }

    // main logic
    attempt_compression(subset_rles_path);
}
