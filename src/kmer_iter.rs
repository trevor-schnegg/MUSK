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
    kmer_length: usize,
}

impl<'a> KmerIter<'a> {
    pub fn from(sequence: &'a [u8], kmer_length: usize) -> Self {
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
            char_iter: sequence.iter(),
            clear_bits: 2_usize.pow((kmer_length * 2) as u32) - 1,
            curr_kmer: 0,
            curr_rev_comp_kmer: 0,
            first_letter_shift: (kmer_length - 1) * 2,
            initialized: false,
            kmer_length,
        }
    }

    fn find_next_kmer(&mut self) -> Option<usize> {
        let mut buffer = 0;
        let mut position = 0_usize;
        while position < self.kmer_length {
            match self.char_iter.next() {
                None => {
                    return None;
                }
                Some(c) => {
                    match self.base2int.get(c) {
                        None => {
                            // Encountered a character that isn't A (a), C (c), G (g), or T (t)
                            buffer = 0;
                            position = 0;
                        }
                        Some(i) => {
                            buffer <<= 2;
                            buffer |= *i;
                            position += 1;
                        }
                    }
                }
            }
        }
        self.curr_kmer = buffer;
        self.curr_rev_comp_kmer = self.reverse_compliment(buffer);
        Some(min(self.curr_kmer, self.curr_rev_comp_kmer))
    }

    /// Only call this if I already have an actual k-mer
    fn reverse_compliment(&self, kmer: usize) -> usize {
        let mut buffer = 0;
        let mut complement_kmer = (!kmer) & self.clear_bits;
        for _ in 0..self.kmer_length {
            // Pop the right-most letter
            let letter = complement_kmer & 3;
            complement_kmer >>= 2;
            // Add to the right of the buffer
            buffer <<= 2;
            buffer |= letter;
        }
        buffer
    }

    pub fn get_curr_kmers(&self) -> (usize, usize) {
        (self.curr_kmer, self.curr_rev_comp_kmer)
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
                Some(character) => {
                    match self.base2int.get(character) {
                        None => {
                            // Encountered a character that isn't A (a), C (c), G (g), or T (t)
                            self.find_next_kmer()
                        }
                        Some(integer) => {
                            self.curr_kmer <<= 2;
                            self.curr_kmer |= *integer;
                            self.curr_kmer &= self.clear_bits;

                            self.curr_rev_comp_kmer >>= 2;
                            self.curr_rev_comp_kmer |=
                                COMPLEMENT[*integer] << self.first_letter_shift;

                            Some(min(self.curr_kmer, self.curr_rev_comp_kmer))
                        }
                    }
                }
            }
        }
    }
}
