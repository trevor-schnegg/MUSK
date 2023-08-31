use crate::binomial::Binomial;
use rug::Float;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

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
    let mut acc = 0;
    for (kmer_position, n) in kmer.iter().rev().enumerate() {
        if *n == 0_i8 {
            // binary is 00, don't need to add anything
            continue
        } else if *n == 1_i8 {
            // binary is 01
            acc += 1_u32 * 2_u32.pow(2_u32 * kmer_position as u32)
        } else if *n == 2_i8 {
            // binary is 10
            acc += 2_u32 * 2_u32.pow(2_u32 * kmer_position as u32)
        } else if *n == 3_i8 {
            // binary is 11
            acc += 3_u32 * 2_u32.pow(2_u32 * kmer_position as u32)
        } else {
            panic!("impossible case")
        }
    }
    acc
}

#[derive(Serialize, Deserialize)]
pub struct Database {
    index2probability: HashMap<usize, f64>,
    multi_accession2index: HashMap<String, usize>,
    kmer_len: usize,
    index2accession: Vec<String>,
    point2occ: HashMap<u32, usize>,
}

impl Database {
    pub fn new(kmer_len: usize) -> Self {
        Database {
            index2probability: HashMap::new(),
            multi_accession2index: HashMap::new(),
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
            .insert(accession_index, self.calculate_probability(insert_count));
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
                Some(index) => {
                    let accession = self.index2accession.get(*index).unwrap();
                    if self.multi_accession2index.contains_key(accession) {
                        for index in accession.split("&").map(|x| x.parse::<usize>().unwrap()) {
                            match index_to_hit_counts.get_mut(&index) {
                                None => {index_to_hit_counts.insert(index, 1);}
                                Some(count) => {*count += 1}
                            }
                        }
                    } else {
                        match index_to_hit_counts.get_mut(index) {
                            None => {index_to_hit_counts.insert(*index, 1);}
                            Some(count) => {*count += 1}
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
            Some(index) => {
                if best_prob < needed_probability {
                    Some(self.get_accession_of_index(index))
                } else {
                    None
                }
            }
        }
    }

    fn insert_kmer(&mut self, kmer: u32, accession_index: usize) -> () {
        match self.point2occ.get_mut(&kmer) {
            None => {
                self.point2occ
                    .insert(kmer, accession_index);
            }
            Some(index) => {
                // since we only search a kmer once per sequence, a collision immediately means we
                // might need a new accession (if this collision hasn't happened before)
                let current_accession = self.index2accession.get(index.clone()).unwrap().clone();
                let new_accession = if current_accession.contains("&") {
                    current_accession + "&" + &*accession_index.to_string()
                } else {
                    index.to_string() + "&" + &*accession_index.to_string()
                };
                match self.multi_accession2index.get(&*new_accession) {
                    None => {
                        let new_accession_index = self.index2accession.len();
                        self.index2accession.push(new_accession.clone());
                        self.multi_accession2index.insert(new_accession, new_accession_index);
                        *index = new_accession_index
                    }
                    Some(new_accession_index) => {
                        *index = *new_accession_index;
                    }
                }
            },
        }
    }

    fn calculate_probability(&self, count: usize) -> f64 {
        count as f64 / 4_usize.pow(self.kmer_len as u32) as f64
    }

    fn get_probability_of_index(&self, accession_index: usize) -> f64 {
        *self.index2probability.get(&accession_index).expect(&*format!("accession index {} does not have a probability", accession_index))
    }

    fn get_accession_of_index(&self, accession_index: usize) -> String {
        self.index2accession.get(accession_index).unwrap().clone()
    }
}
