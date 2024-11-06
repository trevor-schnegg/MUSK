use std::cmp::min;
use std::slice::Iter;

const COMPLEMENT: [usize; 4] = [3, 2, 1, 0];

fn base2int(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'a' => Some(0),
        b'C' => Some(1),
        b'c' => Some(1),
        b'G' => Some(2),
        b'g' => Some(2),
        b'T' => Some(3),
        b't' => Some(3),
        _ => None,
    }
}

pub struct KmerIter<'a> {
    canonical: bool,
    char_iter: Iter<'a, u8>,
    clear_bits: usize,
    curr_kmer: usize,
    curr_rev_comp_kmer: usize,
    first_letter_shift: usize,
    initialized: bool,
    kmer_length: usize,
}

impl<'a> KmerIter<'a> {
    pub fn from(sequence: &'a [u8], kmer_length: usize, canonical: bool) -> Self {
        KmerIter {
            canonical,
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
                Some(char) => {
                    match base2int(*char) {
                        Some(bit_representation) => {
                            buffer <<= 2;
                            buffer |= bit_representation;
                            position += 1;
                        }
                        None => {
                            // Encountered a character that isn't A (a), C (c), G (g), or T (t)
                            buffer = 0;
                            position = 0;
                        }
                    }
                }
            }
        }
        self.curr_kmer = buffer;
        self.curr_rev_comp_kmer = self.reverse_compliment(buffer);
        if self.canonical {
            Some(min(self.curr_kmer, self.curr_rev_comp_kmer))
        } else {
            Some(self.curr_kmer)
        }
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
                Some(char) => {
                    match base2int(*char) {
                        Some(bit_representation) => {
                            self.curr_kmer <<= 2;
                            self.curr_kmer |= bit_representation;
                            self.curr_kmer &= self.clear_bits;

                            self.curr_rev_comp_kmer >>= 2;
                            self.curr_rev_comp_kmer |=
                                COMPLEMENT[bit_representation] << self.first_letter_shift;

                            if self.canonical {
                                Some(min(self.curr_kmer, self.curr_rev_comp_kmer))
                            } else {
                                Some(self.curr_kmer)
                            }
                        }
                        None => {
                            // Encountered a character that isn't A (a), C (c), G (g), or T (t)
                            self.find_next_kmer()
                        }
                    }
                }
            }
        }
    }
}
