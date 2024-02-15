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
                    if *c == b'A' || *c == b'a' {
                        // binary is 00
                        buf <<= 2;
                        kmers_included += 1;
                    } else if *c == b'C' || *c == b'c' {
                        // binary is 01
                        buf <<= 2;
                        buf |= 1;
                        kmers_included += 1;
                    } else if *c == b'G' || *c == b'g' {
                        // binary is 10
                        buf <<= 2;
                        buf |= 2;
                        kmers_included += 1;
                    } else if *c == b'T' || *c == b't' {
                        // binary is 11
                        buf <<= 2;
                        buf |= 3;
                        kmers_included += 1;
                    } else {
                        // Encountered a character that isn't A (a), C (c), G (g), or T (t)
                        buf = 0;
                        kmers_included = 0;
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
                Some(c) => {
                    if *c == b'A' || *c == b'a' {
                        // binary is 00
                        self.curr_kmer <<= 2;
                        self.curr_kmer &= self.clear_bits;
                        Some(self.curr_kmer)
                    } else if *c == b'C' || *c == b'c' {
                        // binary is 01
                        self.curr_kmer <<= 2;
                        self.curr_kmer |= 1;
                        self.curr_kmer &= self.clear_bits;
                        Some(self.curr_kmer)
                    } else if *c == b'G' || *c == b'g' {
                        // binary is 10
                        self.curr_kmer <<= 2;
                        self.curr_kmer |= 2;
                        self.curr_kmer &= self.clear_bits;
                        Some(self.curr_kmer)
                    } else if *c == b'T' || *c == b't' {
                        // binary is 11
                        self.curr_kmer <<= 2;
                        self.curr_kmer |= 3;
                        self.curr_kmer &= self.clear_bits;
                        Some(self.curr_kmer)
                    } else {
                        // Encountered a character that isn't A (a), C (c), G (g), or T (t)
                        self.initialize()
                    }
                }
            }
        }
    }
}
