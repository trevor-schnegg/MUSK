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

    #[clap(short, long)]
    // Invoked with --comp-freq or -c
    comp_freq: bool,
}

// Constant defines the number of bits we are willing to flip in an unencoded block
const FLIP_TOLERANCE: u16 = 1;

#[derive(Default)]
struct CompressionStats {
    in_num_elements: u32,
    in_num_compressed: u32,
    in_num_uncompressed: u32,

    out_num_elements: u32,
    out_num_compressed: u32,
    out_num_uncompressed: u32,

    req_bit_flips: u32,
    total_set_bits: u32,
}

impl CompressionStats {
    fn new() -> Self {
        CompressionStats::default()
    }

    fn add_input_stats(&mut self, total: u32, compressed: u32, uncompressed: u32) {
        self.in_num_elements += total;
        self.in_num_compressed += compressed;
        self.in_num_uncompressed += uncompressed;
    }

    fn add_output_stats(&mut self, total: u32, compressed: u32, uncompressed: u32) {
        self.out_num_elements += total;
        self.out_num_compressed += compressed;
        self.out_num_uncompressed += uncompressed;
    }

    fn add_bit_flips(&mut self, flips: u32) {
        self.req_bit_flips += flips;
    }

    fn add_set_bits(&mut self, bits: u32) {
        self.total_set_bits += bits;
    }

    fn print(&self) {
        println!("Input Statistics:");
        println!("Total blocks: {}\nCompressed blocks: {}\nUncompressed blocks: {}",
            self.in_num_elements,
            self.in_num_compressed,
            self.in_num_uncompressed);
        println!("\nOutput Statistics:");
        println!("Total blocks: {}\nCompressed blocks: {}\nUncompressed blocks: {}",
            self.out_num_elements,
            self.out_num_compressed,
            self.out_num_uncompressed);
        println!("\nData Loss:");
        println!("{} bits flipped of {} total set bits", self.req_bit_flips, self.total_set_bits);
        if self.in_num_elements > 0 && self.total_set_bits > 0 {
            let percent_compression = (1.0 - (self.out_num_elements as f64 / self.in_num_elements as f64)) * 100.0;
            let data_lost = (self.req_bit_flips as f64 / self.total_set_bits as f64) * 100.0;
            println!("\nSummary Statistics:");
            println!("{:.2}% reduction in size", percent_compression);
            println!("{:.2}% of data lost", data_lost);
        }
    }
}

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
    for subset in subset_rles.iter() {
        let subset_vector = subset.1.get_vector().into_iter().map(|run| {
            Run::from_u16(*run)
        }).collect_vec();

        // iterate through each element in the subset_vector
        for (i, current_element) in subset_vector.iter().enumerate() {
            let mut num_neighbors = 0;
            // determine if subset_vector[i] is an uncompressed element
            match current_element {
                Run::Uncompressed(bits) => {
                    // determine number of neighboring compressed zero blocks
                    let count_ones = bits.count_ones() as usize;
                    if i > 0 {
                        match subset_vector[i - 1] {
                            Run::Zeros(_) => {num_neighbors += 1},
                            _ => {},
                        }
                    }
                    
                    if i < subset_vector.len() - 1 {
                        match subset_vector[i + 1] {
                            Run::Zeros(_) => {num_neighbors += 1},
                            _ => {},
                        }
                    }
                    // update result
                    freq_matrix[count_ones - 1][num_neighbors] += 1;
                },
                _ => {},
            };
        }
    }

    // output result
    for row in &freq_matrix {
        println!("{:?}", row);
    }

    println!("");
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
 *          - Method does NOT currently handle the case where neighbor blocks can not hold the additional bits 
 *            individually, but could hold the additional bits if they were divided between the two neighbors.
 *
 * Method outputs the following statistics:
 *      - Input size (total number of blocks, number of compressed, number of uncompressed)
 *      - Final size (total number of blocks, number of compressed, number of uncompressed)
 *      - Data loss (number of bits flipped vs. total number of set bits)
 *      - Summary statistics (percent improved compression, percent data loss)
 */
 fn attempt_compression(subset_rles_path: &Path) {
    // load data from file and initialize variables to track statistics
    let subset_rles = load_data_from_file::<Vec<(u32, RunLengthEncoding)>>(subset_rles_path);

    let mut stats = CompressionStats::new();

    // interate through each sample in the subset
    for subset in subset_rles.iter() {
        let subset_vector = subset.1.get_vector().into_iter().map(|run| {
            Run::from_u16(*run)
        }).collect_vec();

        // update compression statistics for input vector
        let num_elements = subset_vector.len() as u32;
        let num_compressed = subset_vector.iter().filter(|&value| Run::to_u16(value) & (1 << 15) == 0).count() as u32;
        let num_uncompressed = subset_vector.iter().filter(|&value| Run::to_u16(value) & (1 << 15) > 0).count() as u32;
        stats.add_input_stats(num_elements, num_compressed, num_uncompressed);

        let mut result_vector = vec![];
        let mut skip_next = false;

        // iterate through each element in the subset_vector, attempting to compress it
        for (i, current_element) in subset_vector.iter().enumerate() {
            if skip_next {
                skip_next = false;
                continue;
            }
        
            match current_element {
                Run::Uncompressed(bits) => {
                    // dealing with uncompressed block, determine number of set bits 
                    let count_ones = bits.count_ones() as u16;
                    stats.add_set_bits(count_ones as u32);
                    if count_ones <= FLIP_TOLERANCE {
                        // within the flip tolerance, determine if a flip would allow for merge(s)
                        let max_encoded_value = (1 << 14) - 1;
                        let mut can_merge_prev = false;
                        let mut can_merge_next = false;
                        let mut merged_count = 15;

                        // check previous neighbor for potential merge
                        if i > 0 {
                            can_merge_prev = match result_vector[result_vector.len() - 1] {
                                Run::Zeros(count) => {
                                    if merged_count + count <= max_encoded_value {
                                        merged_count += count;
                                        true
                                    } else {
                                        false
                                    }
                                },
                                _ => {false},
                            };
                        }

                        // check following neighbor for potential merge
                        if i < subset_vector.len() - 1 {
                            can_merge_next = match subset_vector[i + 1] {
                                Run::Zeros(count) => {
                                    if merged_count + count <= max_encoded_value {
                                        merged_count += count;
                                        true
                                    } else {
                                        false
                                    }
                                },
                                _ => {false},
                            };
                        }

                        skip_next = can_merge_next;
                        if can_merge_prev || can_merge_next {
                            stats.add_bit_flips(count_ones as u32);
                        }

                        // perform merge
                        if can_merge_prev {
                            let result_length = result_vector.len();
                            result_vector[result_length - 1] = Run::from_u16(merged_count);
                        } else if can_merge_next {
                            result_vector.push(Run::from_u16(merged_count));
                        } else {
                            let current_value = Run::to_u16(current_element);
                            result_vector.push(Run::from_u16(current_value));
                        }
                    } else {
                        // block has more than FLIP_TOLERANCE set bits, add it to result_vector
                        let current_value = Run::to_u16(current_element);
                        result_vector.push(Run::from_u16(current_value));
                    }
                },
                Run::Ones(count) => {
                    stats.add_set_bits(*count as u32);
                    let current_value = Run::to_u16(current_element);
                    result_vector.push(Run::from_u16(current_value));
                }
                Run::Zeros(_) => {
                    let current_value = Run::to_u16(current_element);
                    result_vector.push(Run::from_u16(current_value));
                },
            };
        }

        // update compression statistics for final vector
        let num_elements = result_vector.len() as u32;
        let num_compressed = result_vector.iter().filter(|&value| Run::to_u16(value) & (1 << 15) == 0).count() as u32;
        let num_uncompressed = result_vector.iter().filter(|&value| Run::to_u16(value) & (1 << 15) > 0).count() as u32;
        stats.add_output_stats(num_elements, num_compressed, num_uncompressed);

        if !verify_bits(subset_vector, result_vector) {
            println!("Error! number of bits does not match");
        }
    }

    // output result
    stats.print();
}

/*
 * Method performs basic error testing by ensuring that the same number of bits are being represented
 * between the input and compressed vector.
 *
 * Returns true if the number of bits across both vectors is equal, otherwise false.
 */
 fn verify_bits(input: Vec<Run>, compressed: Vec<Run>) -> bool {
    let mut input_bits: u32 = 0;
    let mut result_bits: u32 = 0;

    for input_value in input {
        match input_value {
            Run::Zeros(count) => {
                input_bits += count as u32;
            },
            Run::Ones(count) => {
                input_bits += count as u32;
            },
            Run::Uncompressed(_) => {
                input_bits += 15;
            },
        };
    }

    for compressed_value in compressed {
        match compressed_value {
            Run::Zeros(count) => {
                result_bits += count as u32;
            },
            Run::Ones(count) => {
                result_bits += count as u32;
            },
            Run::Uncompressed(_) => {
                result_bits += 15;
            },
        };
    }

    input_bits == result_bits
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
