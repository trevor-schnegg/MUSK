use std::cmp::min;
use std::collections::HashMap;
use std::slice::Iter;

const COMPLEMENT: [usize; 4] = [3, 2, 1, 0];

pub struct KmerIter<'a> {
    base2int: HashMap<u8, usize>,
    char_iter: Iter<'a, u8>,
    clear_bits: usize,
    curr_kmer: usize,
    curr_rev_comp_kmer: usize,
    first_letter_shift: usize,
    initialized: bool,
    kmer_len: usize,
}

impl<'a> KmerIter<'a> {
    pub fn from(seq: &'a [u8], kmer_len: usize) -> Self {
        let base2int = HashMap::from([
            (b'A', 0_usize),
            (b'a', 0_usize),
            (b'C', 1_usize),
            (b'c', 1_usize),
            (b'G', 2_usize),
            (b'g', 2_usize),
            (b'T', 3_usize),
            (b't', 3_usize),
        ]);
        KmerIter {
            base2int,
            char_iter: seq.iter(),
            clear_bits: 2_usize.pow((kmer_len * 2) as u32) - 1,
            curr_kmer: 0,
            curr_rev_comp_kmer: 0,
            first_letter_shift: (kmer_len - 1) * 2,
            initialized: false,
            kmer_len,
        }
    }

    fn find_next_kmer(&mut self) -> Option<usize> {
        let mut buf = 0;
        let mut pos = 0_usize;
        while pos < self.kmer_len {
            match self.char_iter.next() {
                None => {
                    return None;
                }
                Some(c) => {
                    match self.base2int.get(c) {
                        None => {
                            // Encountered a character that isn't A (a), C (c), G (g), or T (t)
                            buf = 0;
                            pos = 0;
                        }
                        Some(i) => {
                            buf <<= 2;
                            buf |= *i;
                            pos += 1;
                        }
                    }
                }
            }
        }
        self.curr_kmer = buf;
        self.curr_rev_comp_kmer = self.reverse_compliment(buf);
        Some(min(self.curr_kmer, self.curr_rev_comp_kmer))
    }

    /// Only call this if I already have an actual k-mer
    fn reverse_compliment(&self, kmer: usize) -> usize {
        let mut buf = 0;
        let mut comp_kmer = (!kmer) & self.clear_bits;
        for _ in 0..self.kmer_len {
            // Pop the right-most letter
            let letter = comp_kmer & 3;
            comp_kmer >>= 2;
            // Add to the right of the buffer
            buf <<= 2;
            buf |= letter;
        }
        buf
    }
}

impl<'a> Iterator for KmerIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.initialized {
            self.initialized = true;
            self.find_next_kmer()
        } else {
            match self.char_iter.next() {
                None => {
                    return None;
                }
                Some(c) => {
                    match self.base2int.get(c) {
                        None => {
                            // Encountered a character that isn't A (a), C (c), G (g), or T (t)
                            self.find_next_kmer()
                        }
                        Some(i) => {
                            self.curr_kmer <<= 2;
                            self.curr_kmer |= *i;
                            self.curr_kmer &= self.clear_bits;

                            self.curr_rev_comp_kmer >>= 2;
                            self.curr_rev_comp_kmer |= COMPLEMENT[*i] << self.first_letter_shift;

                            Some(min(self.curr_kmer, self.curr_rev_comp_kmer))
                        }
                    }
                }
            }
        }
    }
}
