use bit_iter::BitIter;
use serde::{Deserialize, Serialize};
use std::{ops::Range, slice::Iter};
use tracing::warn;

pub const MAX_RUN: u16 = (1 << 14) - 1;
pub const MAX_UNCOMPRESSED_BITS: usize = 15;

#[derive(Debug, Clone, Copy)]
pub enum Block {
    Zeros(u16),
    Ones(u16),
    Uncompressed(u16),
}

/// Zeros = 0b_00_14-bit-count
/// Ones = 0b_01_14-bit-count
/// Uncompressed = 0b_1_bits

impl Block {
    pub fn to_u16(&self) -> u16 {
        match *self {
            Block::Zeros(count) => count,
            Block::Ones(count) => count | (1 << 14),
            Block::Uncompressed(bits) => bits | (1 << 15),
        }
    }

    pub fn from_u16(value: u16) -> Block {
        if value & 0b_1000_0000_0000_0000 > 0 {
            Block::Uncompressed(value & 0b_0111_1111_1111_1111)
        } else {
            // the first bit is 0
            if value & 0b_0100_0000_0000_0000 > 0 {
                Block::Ones(value & 0b_0011_1111_1111_1111)
            } else {
                Block::Zeros(value & 0b_0011_1111_1111_1111)
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
    blocks: Box<[u16]>,
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

    pub fn num_of_blocks(&self) -> usize {
        self.runs.len()
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
                self.runs.push(Block::Ones(1).to_u16());
            } else {
                // Otherwise, the last run should be a run of ones
                let last_run_raw = self.runs.last_mut().unwrap();

                let current_num_ones = match Block::from_u16(*last_run_raw) {
                    Block::Ones(num) => num,
                    _ => {
                        panic!("The last value was not a run of ones but it should have been");
                    }
                };

                if current_num_ones == MAX_RUN {
                    // Cannot increment the run because it is at max capacity
                    // Push a new block with a run length of 1
                    self.runs.push(Block::Ones(1).to_u16());
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
                    self.runs.push(Block::Zeros(zeros_needed as u16).to_u16());
                    zeros_needed = 0;
                } else {
                    self.runs.push(Block::Zeros(MAX_RUN).to_u16());
                    zeros_needed -= MAX_RUN as usize;
                }
            }

            // After padding with zeros, insert the set bit
            self.runs.push(Block::Ones(1).to_u16());
        }

        // Before ending the insert, set the highest index
        self.highest_index = index;
    }

    pub fn to_rle(self) -> RunLengthEncoding {
        RunLengthEncoding::compress_from(self.runs)
    }
}

impl RunLengthEncoding {
    pub fn num_of_blocks(&self) -> usize {
        self.blocks.len()
    }

    pub fn get_raw_blocks(&self) -> &[u16] {
        &self.blocks
    }

    pub fn into_raw_blocks(self) -> Box<[u16]> {
        self.blocks
    }

    pub fn iter(&self) -> RunLengthEncodingIter {
        RunLengthEncodingIter::from_blocks(&self.blocks)
    }

    pub fn from(blocks: Box<[u16]>) -> RunLengthEncoding {
        RunLengthEncoding { blocks }
    }

    fn compress_from(runs: Vec<u16>) -> Self {
        // The compressed vector that composes the new run length encoding
        let mut compressed_blocks = vec![];

        // A buffer that represents a bit set of at most MAX_UNCOMPRESSED_BITS
        let mut bits_buffer = vec![];

        // A variable to store the current number of bits in the buffer
        // This is NOT equal to the length of the buffer
        let mut num_bits_in_buffer = 0_usize;

        for run in runs.into_iter().map(|x| Block::from_u16(x)) {
            // Get the number of bits in the current run
            let run_bits = match run {
                Block::Ones(count) => count,
                Block::Zeros(count) => count,
                Block::Uncompressed(_) => panic!("tried to call compress on an rle twice"),
            } as usize;

            if num_bits_in_buffer + run_bits < MAX_UNCOMPRESSED_BITS {
                // If the run can be stored in the buffer, simply push it to the buffer
                bits_buffer.push(run);
                num_bits_in_buffer += run_bits;
            } else if num_bits_in_buffer + run_bits == MAX_UNCOMPRESSED_BITS {
                // If this run is added the buffer will be exactly at capacity

                // If the buffer started as empty, there is no advantage (or disadvantage) to storing this run as uncompressed
                // Therefore, this is done for code simplicity and nothing else
                bits_buffer.push(run);
                let decompressed_run = decompress_buffer(&mut bits_buffer);

                assert!(bits_buffer.is_empty());
                num_bits_in_buffer = 0;

                compressed_blocks.push(decompressed_run);
            } else {
                // Adding the next run would be strictly greater than MAX_UNCOMPRESSED_BITS
                if bits_buffer.is_empty() {
                    // If the buffer is empty, this run overhangs MAX_UNCOMPRESSED_BITS
                    // Simply push it to the compressed vector
                    compressed_blocks.push(run);
                } else {
                    // The buffer is not empty AND adding the next run overhangs MAX_UNCOMPRESSED_BITS
                    // Therefore, borrow some bits from the run we are trying to add to the buffer to fill it
                    let num_bits_to_fill_buffer =
                        (MAX_UNCOMPRESSED_BITS - num_bits_in_buffer) as u16;
                    let leftover_num_bits = run_bits as u16 - num_bits_to_fill_buffer;

                    // Split the run into one to fill the buffer and one leftover
                    let (buffer_run_to_add, leftover_run) = if let Block::Ones(_) = run {
                        (
                            Block::Ones(num_bits_to_fill_buffer),
                            Block::Ones(leftover_num_bits),
                        )
                    } else {
                        (
                            Block::Zeros(num_bits_to_fill_buffer),
                            Block::Zeros(leftover_num_bits),
                        )
                    };

                    bits_buffer.push(buffer_run_to_add);
                    let decompressed_run = decompress_buffer(&mut bits_buffer);

                    assert!(bits_buffer.is_empty());
                    num_bits_in_buffer = 0;

                    // Push the decompressed bits to the compressed vector
                    compressed_blocks.push(decompressed_run);

                    // Handle the leftover bits
                    if leftover_num_bits < MAX_UNCOMPRESSED_BITS as u16 {
                        bits_buffer.push(leftover_run);
                        num_bits_in_buffer += leftover_num_bits as usize;
                    } else {
                        compressed_blocks.push(leftover_run);
                    }
                }
            }
        }

        // After the last run is processed, handle the leftover bits in the buffer (if any)
        if bits_buffer.len() <= 1 {
            compressed_blocks.append(&mut bits_buffer);
        } else {
            compressed_blocks.push(decompress_buffer(&mut bits_buffer))
        }

        let blocks = compressed_blocks
            .into_iter()
            .map(|run| run.to_u16())
            .collect::<Box<[u16]>>();

        RunLengthEncoding { blocks }
    }
}

// Takes a buffer of exactly MAX_UNCOMPRESSED_BITS and converts it to a bit set
fn decompress_buffer(buffer: &mut Vec<Block>) -> Block {
    let mut decompressed = 0;
    let mut current_index = 0;
    for run in buffer.iter() {
        match *run {
            Block::Uncompressed(_) => panic!("impossible case reached"),
            Block::Ones(count) => {
                for i in current_index..current_index + count {
                    decompressed |= 1 << i;
                }
                current_index += count;
            }
            Block::Zeros(count) => current_index += count,
        }
    }
    buffer.clear();
    Block::Uncompressed(decompressed)
}

enum CurrBlockIter {
    None,
    Range(Range<usize>),
    BitIter(BitIter<u16>),
}

pub struct RunLengthEncodingIter<'a> {
    curr_i: usize,
    blocks_iter: Iter<'a, u16>,
    curr_block_iter: CurrBlockIter,
}

impl<'a> RunLengthEncodingIter<'a> {
    pub fn from_blocks(blocks: &'a [u16]) -> Self {
        RunLengthEncodingIter {
            curr_i: 0,
            blocks_iter: blocks.iter(),
            curr_block_iter: CurrBlockIter::None,
        }
    }

    fn find_next(&mut self) -> Option<usize> {
        while let Some(run) = self.blocks_iter.next() {
            match Block::from_u16(*run) {
                Block::Zeros(zeroes_count) => {
                    // If the run is zeros, add it to curr_i and look for the next run
                    self.curr_i += zeroes_count as usize;
                }
                Block::Ones(ones_count) => {
                    // If the run is ones, create a new curr_run_iter
                    self.curr_block_iter =
                        CurrBlockIter::Range(self.curr_i..self.curr_i + ones_count as usize);

                    // Increment curr_i
                    self.curr_i += ones_count as usize;

                    // Return the next value
                    match &mut self.curr_block_iter {
                        CurrBlockIter::Range(range) => return Some(range.next().unwrap()),
                        _ => panic!("impossible case"),
                    }
                }
                Block::Uncompressed(bits) => {
                    // Because of design choices, it is possible to have an empty
                    // uncompressed run -- have to handle this case
                    if bits == 0 {
                        self.curr_i += MAX_UNCOMPRESSED_BITS;
                        continue;
                    }

                    // Otherwise, create new curr_run_iter
                    self.curr_block_iter = CurrBlockIter::BitIter(BitIter::from(bits));

                    // Do NOT increment curr_i here. BitIter is 0 indexed and needs the curr_i to
                    // return the correct value

                    // Return the next value
                    match &mut self.curr_block_iter {
                        CurrBlockIter::BitIter(bit_iter) => {
                            return Some(bit_iter.next().unwrap() + self.curr_i)
                        }
                        _ => panic!("impossible case"),
                    }
                }
            }
        }
        // If there is no next block, reutrn None
        None
    }
}

impl<'a> Iterator for RunLengthEncodingIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        match &mut self.curr_block_iter {
            CurrBlockIter::None => self.find_next(),
            CurrBlockIter::Range(range) => match range.next() {
                None => self.find_next(),
                x => x,
            },
            CurrBlockIter::BitIter(bit_iter) => match bit_iter.next() {
                Some(value) => Some(value + self.curr_i),
                None => {
                    // Increment curr_i after a BitIter is exhausted
                    self.curr_i += MAX_UNCOMPRESSED_BITS;
                    self.find_next()
                }
            },
        }
    }
}
