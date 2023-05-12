use std::collections::HashMap;
use std::slice::Iter;
use f128::f128;
use log::debug;

pub struct Database {
    kmer2index: HashMap<Vec<u8>, Vec<u32>>,
    index2count: Vec<usize>,
    index2accession: Vec<String>,
    index2probability: Vec<f128>,
    kmer_len: usize,
}

impl Database {
    pub fn new(kmer_len: usize) -> Self {
        Database {
            kmer2index: HashMap::new(),
            index2count: Vec::new(),
            index2accession: Vec::new(),
            index2probability: Vec::new(),
            kmer_len,
        }
    }

    pub fn insert_records(&mut self, accession: String, records: Vec<String>) -> () {
        self.index2accession.push(accession);
        self.index2count.push(0);
        let index_of_accession = self.index2accession.len() - 1;
        for record in records {
            self.insert_kmers(record, index_of_accession);
        }
    }

    fn insert_kmers(&mut self, seq: String, index_of_accession: usize) -> () {
        let count = self.index2count.get_mut(index_of_accession).unwrap();
        let index_of_accession = index_of_accession as u32;
        for index in 0..(seq.len() - self.kmer_len) {
            let kmer = &seq[index..(index + self.kmer_len)];
            if kmer.contains("N") {
                continue;
            }
            let kmer_bytes = kmer.as_bytes();
            match self.kmer2index.get_mut(kmer_bytes) {
                None => {
                    self.kmer2index.insert(Vec::from(kmer_bytes), vec![index_of_accession]);
                    *count += 1;
                }
                Some(vec) => {
                    if !vec.contains(&index_of_accession) {
                        vec.push(index_of_accession);
                        *count += 1;
                    }
                }
            }
        }
    }

    pub fn query(&self, kmer: &[u8]) -> Option<Iter<u32>> {
        match self.kmer2index.get(kmer) {
            None => {None}
            Some(vec) => {Some(vec.iter())}
        }
    }

    pub fn init_probabilities(&mut self) -> () {
        for count in &self.index2count {
            let probability = f128::from(*count) / f128::from(4_usize.pow(self.kmer_len as u32));
            self.index2probability.push(probability);
        }
        debug!("Probabilities: {:?}", self.index2probability);
    }

    pub fn get_probability_of_index(&self, index: usize) -> f128 {
        *self.index2probability.get(index).unwrap()
    }

    pub fn get_accession_of_index(&self, index: usize) -> &str {
        self.index2accession.get(index).unwrap()
    }
}
