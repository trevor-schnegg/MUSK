use bit_iter::BitIter;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::{slice::Iter, vec::IntoIter};
use tracing::warn;

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
pub struct BuildRunLengthEncoding {
    highest: usize,
    vector: Vec<u16>,
}

#[derive(Serialize, Deserialize)]
pub struct RunLengthEncoding {
    vector: Vec<u16>,
}

impl BuildRunLengthEncoding {
    pub fn new() -> Self {
        BuildRunLengthEncoding {
            highest: 0,
            vector: vec![],
        }
    }

    pub fn get_vector(&self) -> &Vec<u16> {
        &self.vector
    }

    pub fn push(&mut self, value: usize) -> () {
        if self.vector.is_empty() {
            if value == 0 {
                self.vector.push(Run::Ones(1).to_u16());
            } else {
                // need to add zeros until we reach the correct value
                let mut zeros_needed = value;
                while zeros_needed > 0 {
                    if zeros_needed <= MAX_RUN as usize {
                        self.vector.push(Run::Zeros(zeros_needed as u16).to_u16());
                        zeros_needed = 0;
                    } else {
                        self.vector.push(Run::Zeros(MAX_RUN).to_u16());
                        zeros_needed -= MAX_RUN as usize;
                    }
                }
                self.vector.push(Run::Ones(1).to_u16());
                self.highest = value;
            }
        } else if value <= self.highest {
            warn!(
                "Tried to insert {} into an RLE with highest value {}, skipping...",
                value, self.highest
            );
        } else if value == self.highest + 1 {
            let current_run = self.vector.last_mut().unwrap();
            if let Run::Ones(count) = Run::from_u16(*current_run) {
                if count == MAX_RUN {
                    // Push block with a run length of 1
                    self.vector.push(Run::Ones(1).to_u16());
                } else {
                    *current_run += 1;
                }
            } else {
                panic!("The last value was not a run of ones but it should have been");
            }
            self.highest += 1;
        } else {
            // need to add zeros until we get reach the correct value
            let mut zeros_needed = value - self.highest - 1;
            while zeros_needed > 0 {
                if zeros_needed <= MAX_RUN as usize {
                    self.vector.push(Run::Zeros(zeros_needed as u16).to_u16());
                    zeros_needed = 0;
                } else {
                    self.vector.push(Run::Zeros(MAX_RUN).to_u16());
                    zeros_needed -= MAX_RUN as usize;
                }
            }
            self.vector.push(Run::Ones(1).to_u16());
            self.highest = value;
        }
    }

    pub fn to_rle(self) -> RunLengthEncoding {
        let mut rle = RunLengthEncoding {
            vector: self.vector,
        };
        rle.compress();
        rle
    }
}

impl RunLengthEncoding {
    pub fn get_vector(&self) -> &Vec<u16> {
        &self.vector
    }

    fn compress(&mut self) -> () {
        let mut runs_iterator = self.vector.iter().map(|x| Run::from_u16(*x));
        let mut compressed_vector = vec![];
        let mut buffer = vec![];
        let mut num_bits_in_buffer = 0_usize;

        while let Some(run) = runs_iterator.next() {
            let run_bits = match run {
                Run::Ones(count) => count,
                Run::Zeros(count) => count,
                Run::Uncompressed(_) => panic!("tried to call compress on a vector twice"),
            } as usize;

            if num_bits_in_buffer + run_bits < 15 {
                buffer.push(run);
                num_bits_in_buffer += run_bits;
            } else if num_bits_in_buffer + run_bits == 15 {
                // check if it should be uncompressed
                if num_bits_in_buffer == 0 {
                    compressed_vector.push(run);
                } else {
                    buffer.push(run);
                    let decompressed = decompress_buffer(&mut buffer);
                    buffer.clear();
                    num_bits_in_buffer = 0;
                    compressed_vector.push(Run::Uncompressed(decompressed));
                }
            } else {
                // adding the next run would be strictly greater than 15
                if buffer.is_empty() {
                    compressed_vector.push(run);
                } else {
                    // buffer.len() >= 1
                    let fill_buffer_size = 15 - num_bits_in_buffer as u16;
                    let leftover_size = run_bits as u16 - fill_buffer_size;

                    let (run_to_add, leftover_run) = if let Run::Ones(_) = run {
                        (Run::Ones(fill_buffer_size), Run::Ones(leftover_size))
                    } else {
                        (Run::Zeros(fill_buffer_size), Run::Zeros(leftover_size))
                    };

                    buffer.push(run_to_add);
                    let decompressed = decompress_buffer(&mut buffer);

                    compressed_vector.push(Run::Uncompressed(decompressed));

                    buffer.clear();
                    num_bits_in_buffer = 0;

                    if leftover_size < 15 {
                        buffer.push(leftover_run);
                        num_bits_in_buffer += leftover_size as usize;
                    } else {
                        compressed_vector.push(leftover_run);
                    }
                }
            }
        }

        // Handle the leftover bits in the buffer
        if buffer.len() <= 1 {
            compressed_vector.append(&mut buffer);
        } else {
            compressed_vector.push(Run::Uncompressed(decompress_buffer(&buffer)))
        }

        self.vector = compressed_vector
            .into_iter()
            .map(|x| x.to_u16())
            .collect_vec();
    }
}

fn decompress_buffer(buffer: &Vec<Run>) -> u16 {
    let mut decompressed = 0;
    let mut current_index = 0;
    for run in buffer {
        match *run {
            Run::Uncompressed(_) => panic!("impossible case"),
            Run::Ones(count) => {
                for i in current_index..current_index + count {
                    decompressed |= 1 << i;
                }
                current_index += count;
            }
            Run::Zeros(count) => current_index += count,
        }
    }
    decompressed
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
