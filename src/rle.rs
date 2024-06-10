use bit_iter::BitIter;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::{slice::Iter, vec::IntoIter};
use tracing::warn;

const MAX_UNCOMPRESSED_BITS: usize = 15;

#[derive(Debug)]
pub enum Run {
    Zeros(u16),
    Ones(u16),
    Uncompressed(u16),
}

/// Zeros = 0b_00_14-bit-count
/// Ones = 0b_01_14-bit-count
/// Uncompressed = 0b_1_bits

const MAX_RUN: u16 = (1 << 14) - 1;

impl Run {
    pub fn to_u16(&self) -> u16 {
        match *self {
            Run::Zeros(count) => count,
            Run::Ones(count) => count | (1 << 14),
            Run::Uncompressed(bits) => bits | (1 << 15),
        }
    }

    pub fn from_u16(value: u16) -> Run {
        if value & 0b_1000_0000_0000_0000 > 0 {
            Run::Uncompressed(value & 0b_0111_1111_1111_1111)
        } else {
            // the first bit is 0
            if value & 0b_0100_0000_0000_0000 > 0 {
                Run::Ones(value & 0b_0011_1111_1111_1111)
            } else {
                Run::Zeros(value & 0b_0011_1111_1111_1111)
            }
        }
    }
}

#[derive(Clone)]
pub struct NaiveRunLengthEncoding {
    highest_index: usize,
    runs: Vec<u16>,
}

#[derive(Serialize, Deserialize)]
pub struct RunLengthEncoding {
    runs: Vec<u16>,
}

impl NaiveRunLengthEncoding {
    pub fn new() -> Self {
        NaiveRunLengthEncoding {
            highest_index: usize::MAX,
            runs: vec![],
        }
    }

    pub fn get_raw_runs(&self) -> &Vec<u16> {
        &self.runs
    }

    pub fn push(&mut self, index: usize) -> () {
        // This is the smallest value I am allowed to insert
        let next_sequential_index = self.highest_index.overflowing_add(1).0;

        if index < next_sequential_index {
            warn!(
                "Tried to insert {} into an RLE with highest value {}, skipping...",
                index, self.highest_index
            );
            return;
        } else if index == next_sequential_index {
            // In this case, I should try to extend the last run of ones (if possible)
            if self.runs.is_empty() {
                // Handle the special case where this is index 0
                self.runs.push(Run::Ones(1).to_u16());
            } else {
                // Otherwise, the last run should be a run of ones
                let last_run_raw = self.runs.last_mut().unwrap();

                let current_num_ones = match Run::from_u16(*last_run_raw) {
                    Run::Ones(num) => num,
                    _ => {
                        panic!("The last value was not a run of ones but it should have been");
                    }
                };

                if current_num_ones == MAX_RUN {
                    // Cannot increment the run because it is at max capacity
                    // Push a new block with a run length of 1
                    self.runs.push(Run::Ones(1).to_u16());
                } else {
                    *last_run_raw += 1;
                }
            }
        } else {
            // Need to pad with zeros until we reach the desired value
            let mut zeros_needed = if self.runs.is_empty() {
                index
            } else {
                index - self.highest_index - 1
            };

            while zeros_needed > 0 {
                if zeros_needed <= MAX_RUN as usize {
                    self.runs.push(Run::Zeros(zeros_needed as u16).to_u16());
                    zeros_needed = 0;
                } else {
                    self.runs.push(Run::Zeros(MAX_RUN).to_u16());
                    zeros_needed -= MAX_RUN as usize;
                }
            }

            // After padding with zeros, insert the set bit
            self.runs.push(Run::Ones(1).to_u16());
        }

        // Finally, set the highest index as the index we just inserted
        self.highest_index = index;
    }

    pub fn to_rle(self) -> RunLengthEncoding {
        RunLengthEncoding::compress_from(self.runs)
    }
}

impl RunLengthEncoding {
    pub fn get_raw_runs(&self) -> &Vec<u16> {
        &self.runs
    }

    fn compress_from(runs: Vec<u16>) -> Self {
        // The compressed vector that composes the new run length encoding
        let mut compressed_runs = vec![];
        // A buffer that represents a bit set of at most MAX_UNCOMPRESSED_BITS
        let mut bits_buffer = vec![];
        // A variable to store the current number of bits in the buffer
        // This is NOT equal to the length of the buffer
        let mut num_bits_in_buffer = 0_usize;

        for run in runs.into_iter().map(|x| Run::from_u16(x)) {
            // Get the number of bits in the current run
            let run_bits = match run {
                Run::Ones(count) => count,
                Run::Zeros(count) => count,
                Run::Uncompressed(_) => panic!("tried to call compress on an rle twice"),
            } as usize;

            if num_bits_in_buffer + run_bits < MAX_UNCOMPRESSED_BITS {
                // If the run can be stored in the buffer, simply push it to the buffer
                bits_buffer.push(run);
                num_bits_in_buffer += run_bits;
            } else if num_bits_in_buffer + run_bits == MAX_UNCOMPRESSED_BITS {
                // If this run is added the buffer will be exactly at capacity

                // If the buffer started as empty, there is no advantage to storing this run as uncompressed
                // Therefore, this is done for code simplicity and nothing else
                bits_buffer.push(run);
                let decompressed_run = decompress_buffer(&mut bits_buffer);

                assert!(bits_buffer.is_empty());
                num_bits_in_buffer = 0;

                compressed_runs.push(decompressed_run);
            } else {
                // Adding the next run would be strictly greater than MAX_UNCOMPRESSED_BITS
                if bits_buffer.is_empty() {
                    // If the buffer is empty, this run overhangs MAX_UNCOMPRESSED_BITS
                    // Simply push it to the compressed vector
                    compressed_runs.push(run);
                } else {
                    // The buffer is not empty AND adding the next run overhangs MAX_UNCOMPRESSED_BITS
                    // Therefore, borrow some bits from the run we are trying to add to the buffer to fill it
                    let num_bits_to_fill_buffer =
                        (MAX_UNCOMPRESSED_BITS - num_bits_in_buffer) as u16;
                    let leftover_num_bits = run_bits as u16 - num_bits_to_fill_buffer;

                    // Split the run into one to fill the buffer and one leftover
                    let (buffer_run_to_add, leftover_run) = if let Run::Ones(_) = run {
                        (
                            Run::Ones(num_bits_to_fill_buffer),
                            Run::Ones(leftover_num_bits),
                        )
                    } else {
                        (
                            Run::Zeros(num_bits_to_fill_buffer),
                            Run::Zeros(leftover_num_bits),
                        )
                    };

                    bits_buffer.push(buffer_run_to_add);
                    let decompressed_run = decompress_buffer(&mut bits_buffer);

                    assert!(bits_buffer.is_empty());
                    num_bits_in_buffer = 0;

                    // Push the decompressed bits to the compressed vector
                    compressed_runs.push(decompressed_run);

                    // Handle the leftover bits
                    if leftover_num_bits < MAX_UNCOMPRESSED_BITS as u16 {
                        bits_buffer.push(leftover_run);
                        num_bits_in_buffer += leftover_num_bits as usize;
                    } else {
                        compressed_runs.push(leftover_run);
                    }
                }
            }
        }

        // After the last run is processed, handle the leftover bits in the buffer (if any)
        if bits_buffer.len() <= 1 {
            compressed_runs.append(&mut bits_buffer);
        } else {
            compressed_runs.push(decompress_buffer(&mut bits_buffer))
        }

        RunLengthEncoding {
            runs: compressed_runs
                .into_iter()
                .map(|run| run.to_u16())
                .collect_vec(),
        }
    }
}

// Takes a buffer of exactly MAX_UNCOMPRESSED_BITS and converts it to a bit set
fn decompress_buffer(buffer: &mut Vec<Run>) -> Run {
    let mut decompressed = 0;
    let mut current_index = 0;
    for run in buffer.iter() {
        match *run {
            Run::Uncompressed(_) => panic!("impossible case reached"),
            Run::Ones(count) => {
                for i in current_index..current_index + count {
                    decompressed |= 1 << i;
                }
                current_index += count;
            }
            Run::Zeros(count) => current_index += count,
        }
    }
    buffer.clear();
    Run::Uncompressed(decompressed)
}

pub struct RunLengthEncodingIter<'a> {
    curr_i: usize,
    runs_iter: Iter<'a, u16>,
    curr_run_iter: IntoIter<usize>,
}

impl<'a> RunLengthEncodingIter<'a> {
    pub fn from_runs_vector(vec: &'a Vec<u16>) -> Self {
        let mut curr_i = 0_usize;
        let mut runs_iter = vec.into_iter();
        while let Some(run) = runs_iter.next() {
            match Run::from_u16(*run) {
                Run::Zeros(count) => {
                    curr_i += count as usize;
                }
                Run::Ones(count) => {
                    return RunLengthEncodingIter {
                        curr_i: curr_i + count as usize,
                        runs_iter,
                        curr_run_iter: (curr_i..curr_i + count as usize).collect_vec().into_iter(),
                    }
                }
                Run::Uncompressed(bits) => {
                    return RunLengthEncodingIter {
                        curr_i: curr_i + 15,
                        runs_iter,
                        curr_run_iter: BitIter::from(bits)
                            .map(|x| x + curr_i)
                            .collect_vec()
                            .into_iter(),
                    }
                }
            }
        }
        RunLengthEncodingIter {
            curr_i,
            runs_iter,
            curr_run_iter: (0_usize..0_usize).collect_vec().into_iter(),
        }
    }
}

impl<'a> Iterator for RunLengthEncodingIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        match self.curr_run_iter.next() {
            Some(n) => Some(n),
            None => {
                while let Some(run) = self.runs_iter.next() {
                    match Run::from_u16(*run) {
                        Run::Zeros(count) => {
                            self.curr_i += count as usize;
                        }
                        Run::Ones(count) => {
                            self.curr_run_iter = (self.curr_i..self.curr_i + count as usize)
                                .collect_vec()
                                .into_iter();
                            self.curr_i += count as usize;
                            return self.curr_run_iter.next();
                        }
                        Run::Uncompressed(bits) => {
                            self.curr_run_iter = BitIter::from(bits)
                                .map(|x| x + self.curr_i)
                                .collect_vec()
                                .into_iter();
                            self.curr_i += 15;
                            return self.curr_run_iter.next();
                        }
                    }
                }
                None
            }
        }
    }
}
