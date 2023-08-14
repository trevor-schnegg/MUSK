use std::collections::{HashMap, HashSet};
use rug::Float;
use crate::intervals::Intervals;

pub struct Database {
    accession2probability: HashMap<String, Float>,
    kmer_len: usize,
    intervals: Intervals
}

impl Database {

    pub fn new(kmer_len: usize) -> Self {
        Database {
            accession2probability: HashMap::new(),
            kmer_len,
            intervals: Intervals::new(),
        }
    }

    pub fn push_interval(&mut self, interval: (usize, usize, String)) -> () {
        self.intervals.push(interval)
    }

    pub fn get_accession_of_index(&self, index: usize) -> String {
        self.intervals.get_accession_of_index(index)
    }

    pub fn calculate_probability(&self, sequences: (&str, &str)) -> Float {
        let mut set = HashSet::new();
        for index in 0..(sequences.0.len() - self.kmer_len) {
            let kmer_1 = &sequences.0[index..(index + self.kmer_len)];
            if !kmer_1.contains("N") {
                set.insert(kmer_1.as_bytes());
            }
            let kmer_2 = &sequences.1[index..(index + self.kmer_len)];
            if !kmer_2.contains("N") {
                set.insert(kmer_2.as_bytes());
            }
        }
        Float::with_val(256, set.len()) / Float::with_val(256, 4_usize.pow(self.kmer_len as u32))
    }

    pub fn set_probabilities(&mut self, accession2probabilities: HashMap<String, Float>) -> () {
        self.accession2probability = accession2probabilities
    }

    pub fn get_probability(&self, accession: &str) -> Float {
        self.accession2probability.get(accession).unwrap().clone()
    }

}