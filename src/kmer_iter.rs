use crate::utility::two_bit_dna_representation;
use std::slice::Iter;

pub struct KmerIter<'a> {
    char_iter: Iter<'a, u8>,
    kmer_len: usize,
    clear_bits: usize,
    curr_kmer: usize,
    initialized: bool,
}

impl<'a> KmerIter<'a> {
    pub fn from(seq: &'a [u8], kmer_len: usize) -> Self {
        KmerIter {
            char_iter: seq.iter(),
            kmer_len,
            clear_bits: 2_usize.pow((kmer_len * 2) as u32) - 1,
            curr_kmer: 0,
            initialized: false,
        }
    }

    fn initialize(&mut self) -> Option<usize> {
        let mut buf = 0;
        let mut kmers_included = 0_usize;
        while kmers_included < self.kmer_len {
            match self.char_iter.next() {
                None => {
                    return None;
                }
                Some(c) => {
                    match two_bit_dna_representation(c) {
                        None => {
                            // Encountered a character that isn't A (a), C (c), G (g), or T (t)
                            buf = 0;
                            kmers_included = 0;
                        }
                        Some(n) => {
                            buf <<= 2;
                            buf |= n;
                            kmers_included += 1;
                        }
                    }
                }
            }
        }
        self.curr_kmer = buf;
        Some(self.curr_kmer)
    }
}

impl<'a> Iterator for KmerIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.initialized {
            self.initialized = true;
            self.initialize()
        } else {
            match self.char_iter.next() {
                None => {
                    return None;
                }
                Some(c) => match two_bit_dna_representation(c) {
                    None => self.initialize(),
                    Some(n) => {
                        self.curr_kmer <<= 2;
                        self.curr_kmer |= n;
                        self.curr_kmer &= self.clear_bits;
                        Some(self.curr_kmer)
                    }
                },
            }
        }
    }
}
