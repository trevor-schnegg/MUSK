use itertools::Itertools;
use log::warn;


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
            Run::Zeros(count) => {count},
            Run::Ones(count) => {count | (1 << 14)},
            Run::Uncompressed(bits) => {bits | (1 << 15)},
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
            warn!("Tried to insert {} into an RLE with highest value {}, skipping...", value, self.highest);
        } else if value == self.highest + 1 {
            let current_run = self.vector.last_mut().unwrap();
            if let Run::Ones(count) = Run::from_u16(*current_run) {
                if count == MAX_RUN {
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
        let mut current_buffer_size = 0_usize;
        while let Some(run) = runs_iterator.next() {
            let run_size = match run {
                Run::Ones(count) => count,
                Run::Zeros(count) => count,
                Run::Uncompressed(_) => panic!("tried to compress an already compressed vector"),
            } as usize;
            if current_buffer_size + run_size < 15 {
                buffer.push(run);
                current_buffer_size += run_size;
            } else if current_buffer_size + run_size == 15 {
                // check if it should be uncompressed
                if buffer.len() == 0 {
                    compressed_vector.push(run);
                } else {
                    buffer.push(run);
                    let decompressed = decompress(&mut buffer);
                    buffer.clear();
                    current_buffer_size = 0;
                    compressed_vector.push(Run::Uncompressed(decompressed));
                }
            } else {
                // adding the next run would be strictly greater than 15
                if buffer.len() == 0 {
                    compressed_vector.push(run);
                } else if buffer.len() == 1 {
                    compressed_vector.append(&mut buffer);
                    current_buffer_size = 0;
                    if run_size < 15 {
                        buffer.push(run);
                        current_buffer_size = run_size;
                    } else {
                        compressed_vector.push(run);
                    }
                } else {
                    // buffer.len() >= 2
                    let fill_buffer_size = 15 - current_buffer_size as u16;
                    let leftover_size = run_size as u16 - fill_buffer_size;
                    let (run_to_push, leftover_run) = if let Run::Ones(_) = run {
                        (Run::Ones(fill_buffer_size), Run::Ones(leftover_size))
                    } else {
                        (Run::Zeros(fill_buffer_size), Run::Zeros(leftover_size))
                    };
                    buffer.push(run_to_push);
                    let decompressed = decompress(&mut buffer);
                    buffer.clear();
                    current_buffer_size = 0;
                    compressed_vector.push(Run::Uncompressed(decompressed));
                    if leftover_size < 15 {
                        buffer.push(leftover_run);
                        current_buffer_size += leftover_size as usize;
                    } else {
                        compressed_vector.push(leftover_run);
                    }
                }
            }
        }
        if buffer.len() <= 1 {
            compressed_vector.append(&mut buffer);
        } else {
            compressed_vector.push(Run::Uncompressed(decompress(&buffer)))
        }
        self.vector = compressed_vector.into_iter().map(|x| x.to_u16()).collect_vec();
    }
}

fn decompress(buffer: &Vec<Run>) -> u16 {
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
            },
            Run::Zeros(count) => current_index += count,
        }
    }
    decompressed
}
