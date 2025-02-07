use bit_iter::BitIter;
use serde::{Deserialize, Serialize};
use std::slice::Iter;
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
        RunLengthEncoding::allow_uncompressed_from(self.runs)
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

    pub fn block_iters(&self) -> RunLengthEncodingBlockIter {
        RunLengthEncodingBlockIter::from_blocks(&self.blocks)
    }

    pub fn from(blocks: Box<[u16]>) -> RunLengthEncoding {
        RunLengthEncoding { blocks }
    }

    fn allow_uncompressed_from(runs: Vec<u16>) -> Self {
        // The compressed vector that composes the new run length encoding
        let mut blocks_w_uncompressed_allowed = Vec::with_capacity(runs.len());

        // A buffer that represents a bit map of at most MAX_UNCOMPRESSED_BITS
        let mut bits_buffer = vec![];

        // A variable to store the current number of bits in the buffer
        // This is NOT equal to the length of the buffer
        let mut num_bits_in_buffer = 0_usize;

        runs.into_iter().for_each(|naive_u16| {
            let naive_block = Block::from_u16(naive_u16);

            // Get the number of bits in the current run
            let naive_block_bits = match naive_block {
                Block::Ones(count) => count,
                Block::Zeros(count) => count,
                Block::Uncompressed(_) => panic!("tried to allow uncompressed blocks in an RLE that already has uncompressed blocks"),
            } as usize;

            if num_bits_in_buffer + naive_block_bits < MAX_UNCOMPRESSED_BITS {
                // If the run can be stored in the buffer, simply push it to the buffer
                bits_buffer.push(naive_block);
                num_bits_in_buffer += naive_block_bits;

            } else if num_bits_in_buffer + naive_block_bits == MAX_UNCOMPRESSED_BITS {
                // If this run is added the buffer will be exactly at capacity

                if bits_buffer.is_empty() {
                    // If the bits buffer is empty, then a single run exactly filled the buffer
                    // Therefore, it does not matter if we store as a run or as uncompressed
                    // Simply push the original u16 to the new vector
                    blocks_w_uncompressed_allowed.push(naive_u16);

                    // Assert that we've maintained our invariant
                    assert_eq!(num_bits_in_buffer, 0);
                } else {
                    // Otherwise, there was something already in the buffer
                    bits_buffer.push(naive_block);

                    // Now at least 2 items are in the buffer, we can save by using uncompressed
                    blocks_w_uncompressed_allowed.push(create_uncompressed_from(&bits_buffer));

                    // Clear the bits buffer and the number of bits in it
                    bits_buffer.clear(); num_bits_in_buffer = 0;
                }
            } else {
                // Adding the next run would be strictly greater than MAX_UNCOMPRESSED_BITS
                if bits_buffer.is_empty() {
                    // If the buffer is empty, this run overhangs MAX_UNCOMPRESSED_BITS
                    // Simply push the original u16 to the new vector
                    blocks_w_uncompressed_allowed.push(naive_u16);

                    // Assert that we've maintained our invariant
                    assert_eq!(num_bits_in_buffer, 0);
                } else {
                    // The buffer is not empty AND adding the next run overhangs MAX_UNCOMPRESSED_BITS
                    // Therefore, borrow some bits from the run we are trying to add to the buffer to fill it
                    let num_bits_to_fill_buffer =
                        (MAX_UNCOMPRESSED_BITS - num_bits_in_buffer) as u16;
                    let leftover_num_bits = naive_block_bits as u16 - num_bits_to_fill_buffer;

                    // Split into 2 blocks
                    let (block_to_add_to_buffer, leftover_block) = if let Block::Ones(_) = naive_block {
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

                    // Push the new block to fill the buffer
                    bits_buffer.push(block_to_add_to_buffer);

                    // Create an uncompressed block from the buffer and push to the new vector
                    blocks_w_uncompressed_allowed.push(create_uncompressed_from(&bits_buffer));

                    // Clear the bits buffer and the number of bits in it
                    bits_buffer.clear(); num_bits_in_buffer = 0;

                    // Handle the leftover bits
                    if leftover_num_bits < MAX_UNCOMPRESSED_BITS as u16 {
                        // If the run can be stored in the buffer, simply push it to the buffer
                        bits_buffer.push(leftover_block);
                        num_bits_in_buffer += leftover_num_bits as usize;
                    } else {
                        // Otherwise, push the run as a u16 to the new vector
                        blocks_w_uncompressed_allowed.push(leftover_block.to_u16());
                    }
                }
            }
        });

        // After the last run is processed, handle the leftover bits in the buffer (if any)
        if bits_buffer.len() <= 1 {
            bits_buffer
                .into_iter()
                .for_each(|run| blocks_w_uncompressed_allowed.push(run.to_u16()));
        } else {
            blocks_w_uncompressed_allowed.push(create_uncompressed_from(&bits_buffer))
        }

        RunLengthEncoding {
            blocks: blocks_w_uncompressed_allowed.into_boxed_slice(),
        }
    }

    pub fn collect_indices(&self) -> Vec<u32> {
        // Create the blocks iterator
        let mut blocks_iter = self
            .blocks
            .iter()
            .map(|block_u16| Block::from_u16(*block_u16));

        // Initialize curr_i and the return value
        let mut curr_i = 0_u32;
        let mut indices = vec![];

        while let Some(block) = blocks_iter.next() {
            match block {
                Block::Zeros(zeroes_count) => curr_i += zeroes_count as u32,
                Block::Ones(ones_count) => {
                    let ones_count = ones_count as u32;
                    indices.extend(curr_i..curr_i + ones_count);
                    curr_i += ones_count;
                }
                Block::Uncompressed(bits) => {
                    indices.extend(BitIter::from(bits).map(|i| i as u32 + curr_i));
                    curr_i += MAX_UNCOMPRESSED_BITS as u32;
                }
            }
        }
        indices
    }
}

// Takes a buffer of exactly MAX_UNCOMPRESSED_BITS and converts it to a bit set
fn create_uncompressed_from(buffer: &Vec<Block>) -> u16 {
    let mut uncompressed = 0;
    let mut current_index = 0;
    buffer.iter().for_each(|run| match *run {
        Block::Uncompressed(_) => panic!("impossible case"),
        Block::Ones(count) => {
            assert!((current_index + count - 1) as usize <= MAX_UNCOMPRESSED_BITS);
            for i in current_index..current_index + count {
                uncompressed |= 1 << i;
            }
            current_index += count;
        }
        Block::Zeros(count) => current_index += count,
    });

    // Return the uncompressed block
    Block::Uncompressed(uncompressed).to_u16()
}

pub enum BlockIter {
    Range((usize, usize)),
    BitIter((BitIter<u16>, usize)),
}

pub struct RunLengthEncodingBlockIter<'a> {
    curr_i: usize,
    blocks_iter: Iter<'a, u16>,
}

impl<'a> RunLengthEncodingBlockIter<'a> {
    pub fn from_blocks(blocks: &'a [u16]) -> Self {
        RunLengthEncodingBlockIter {
            curr_i: 0,
            blocks_iter: blocks.iter(),
        }
    }

    fn find_next(&mut self) -> Option<BlockIter> {
        while let Some(run) = self.blocks_iter.next() {
            match Block::from_u16(*run) {
                Block::Zeros(zeroes_count) => {
                    // If the run is zeros, add it to curr_i and look for the next run
                    self.curr_i += zeroes_count as usize;
                }
                Block::Ones(ones_count) => {
                    // If the run is ones, create a new curr_run_iter
                    let return_block =
                        BlockIter::Range((self.curr_i, self.curr_i + ones_count as usize));

                    // Increment curr_i
                    self.curr_i += ones_count as usize;

                    // Return the next value
                    return Some(return_block);
                }
                Block::Uncompressed(bits) => {
                    // Because of design choices, it is possible to have an empty
                    // uncompressed run -- have to handle this case
                    if bits == 0 {
                        self.curr_i += MAX_UNCOMPRESSED_BITS;
                        continue;
                    }

                    // Otherwise, create new curr_run_iter
                    let return_block = BlockIter::BitIter((BitIter::from(bits), self.curr_i));

                    // Do NOT increment curr_i here. BitIter is 0 indexed and needs the curr_i to
                    // return the correct value
                    self.curr_i += MAX_UNCOMPRESSED_BITS;

                    // Return the next value
                    return Some(return_block);
                }
            }
        }
        // If there is no next block, reutrn None
        None
    }
}

impl<'a> Iterator for RunLengthEncodingBlockIter<'a> {
    type Item = BlockIter;

    fn next(&mut self) -> Option<Self::Item> {
        self.find_next()
    }
}
