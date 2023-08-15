use crate::binomial::Binomial;
use rug::Float;
use std::collections::{HashMap, HashSet};
use log::debug;

fn base_to_binary(character: char) -> i8 {
    if character == 'A' {
        0
    } else if character == 'C' {
        1
    } else if character == 'G' {
        2
    } else if character == 'T' {
        3
    } else {
        -1
    }
}

fn convert_vec_i8_to_u32(kmer: &[i8]) -> u32 {
    let mut acc = String::new();
    for  n in kmer {
        if *n == 0_i8 {
            acc += "00"
        } else if *n == 1_i8 {
            acc += "01"
        } else if *n == 2_i8 {
            acc += "10"
        } else if *n == 3_i8 {
            acc += "11"
        } else {
            panic!("impossible case")
        }
    }
    u32::from_str_radix(&*acc, 2).expect("could not convert string to u32")
}

pub struct Database {
    index2probability: Vec<Float>,
    kmer_len: usize,
    index2accession: Vec<String>,
    point2occ: HashMap<u32, Vec<usize>>,
}

impl Database {
    pub fn new(kmer_len: usize) -> Self {
        Database {
            index2probability: Vec::new(),
            kmer_len,
            index2accession: Vec::new(),
            point2occ: HashMap::new(),
        }
    }

    pub fn insert_record(
        &mut self,
        forward_seq: String,
        reverse_seq: String,
        accession: String,
    ) -> () {
        let accession_index = self.index2accession.len();
        self.index2accession.push(accession);

        let mut kmer_set = HashSet::new();
        let forward_chars = forward_seq
            .chars()
            .map(|x| base_to_binary(x))
            .collect::<Vec<i8>>();
        let reverse_chars = reverse_seq
            .chars()
            .map(|x| base_to_binary(x))
            .collect::<Vec<i8>>();

        for (kmer_1, kmer_2) in forward_chars
            .windows(self.kmer_len)
            .zip(reverse_chars.windows(self.kmer_len))
        {
            if kmer_1.contains(&-1) || kmer_2.contains(&-1) {
                continue;
            }
            let (kmer_1, kmer_2) = (convert_vec_i8_to_u32(kmer_1), convert_vec_i8_to_u32(kmer_2));
            kmer_set.insert(kmer_1);
            kmer_set.insert(kmer_2);
        }

        let insert_count = kmer_set.len();
        for kmer in kmer_set {
            self.insert_kmer(kmer, accession_index);
        }

        self.index2probability
            .push(self.calculate_probability(insert_count));
    }

    pub fn query_read(&self, read: String) -> Option<String> {
        let mut kmer_set = HashSet::new();
        let read = read.chars().map(|x| base_to_binary(x)).collect::<Vec<i8>>();
        for kmer in read.windows(self.kmer_len) {
            kmer_set.insert(convert_vec_i8_to_u32(kmer));
        }
        let num_queries = kmer_set.len() as u64;

        let mut index_to_hit_counts = HashMap::new();
        for kmer in kmer_set {
            match self.point2occ.get(&kmer) {
                None => continue,
                Some(vec) => {
                    for index in vec {
                        match index_to_hit_counts.get_mut(index) {
                            None => {
                                index_to_hit_counts.insert(*index, 1_u64);
                            },
                            Some(count) => *count += 1,
                        }
                    }
                }
            }
        }
        self.calculate_result(index_to_hit_counts, num_queries)
    }

    fn calculate_result(
        &self,
        index_to_hit_counts: HashMap<usize, u64>,
        num_queries: u64,
    ) -> Option<String> {
        let needed_probability = Float::with_val(256, 1.0e-100);
        let mut best_prob = Float::with_val(256, 1.0);
        let mut best_prob_index = None;
        for (accession_index, num_hits) in index_to_hit_counts {
            let binomial = Binomial::new(
                Float::with_val(256, self.get_probability_of_index(accession_index)),
                num_queries,
            )
            .unwrap();
            let prob = binomial.sf(num_hits).abs();
            if prob < best_prob {
                best_prob = prob;
                best_prob_index = Some(accession_index);
            }
        }
        match best_prob_index {
            None => None,
            Some(index) => if best_prob < needed_probability { Some(self.get_accession_of_index(index)) } else { None },
        }
    }

    fn insert_kmer(&mut self, kmer: u32, accession_index: usize) -> () {
        match self.point2occ.get_mut(&kmer) {
            None => {self.point2occ.insert(kmer, Vec::from(vec![accession_index]));}
            Some(vec) => { vec.push(accession_index) }
        }
    }

    fn calculate_probability(&self, count: usize) -> Float {
        Float::with_val(256, count) / Float::with_val(256, 4_usize.pow(self.kmer_len as u32))
    }

    fn get_probability_of_index(&self, accession_index: usize) -> Float {
        self.index2probability.get(accession_index).unwrap().clone()
    }

    fn get_accession_of_index(&self, accession_index: usize) -> String {
        self.index2accession.get(accession_index).unwrap().clone()
    }
}
