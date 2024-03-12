use std::cmp::min;
use std::collections::HashMap;
use std::slice::Iter;

const COMPLEMENT: [usize; 4] = [3, 2, 1, 0];

pub struct KmerIter<'a> {
    base2integer: HashMap<u8, usize>,
    character_iterator: Iter<'a, u8>,
    clear_bits: usize,
    current_kmer: usize,
    current_reverse_complement_kmer: usize,
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
            base2integer: base2int,
            character_iterator: sequence.iter(),
            clear_bits: 2_usize.pow((kmer_length * 2) as u32) - 1,
            current_kmer: 0,
            current_reverse_complement_kmer: 0,
            first_letter_shift: (kmer_length - 1) * 2,
            initialized: false,
            kmer_length,
        }
    }

    fn find_next_kmer(&mut self) -> Option<usize> {
        let mut buffer = 0;
        let mut position = 0_usize;
        while position < self.kmer_length {
            match self.character_iterator.next() {
                None => {
                    return None;
                }
                Some(c) => {
                    match self.base2integer.get(c) {
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
        self.current_kmer = buffer;
        self.current_reverse_complement_kmer = self.reverse_compliment(buffer);
        Some(min(self.current_kmer, self.current_reverse_complement_kmer))
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
}

impl<'a> Iterator for KmerIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.initialized {
            self.initialized = true;
            self.find_next_kmer()
        } else {
            match self.character_iterator.next() {
                None => {
                    return None;
                }
                Some(character) => {
                    match self.base2integer.get(character) {
                        None => {
                            // Encountered a character that isn't A (a), C (c), G (g), or T (t)
                            self.find_next_kmer()
                        }
                        Some(integer) => {
                            self.current_kmer <<= 2;
                            self.current_kmer |= *integer;
                            self.current_kmer &= self.clear_bits;

                            self.current_reverse_complement_kmer >>= 2;
                            self.current_reverse_complement_kmer |=
                                COMPLEMENT[*integer] << self.first_letter_shift;

                            Some(min(self.current_kmer, self.current_reverse_complement_kmer))
                        }
                    }
                }
            }
        }
    }
}
