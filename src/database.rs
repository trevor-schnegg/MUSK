use crate::accession_tree::AccessionTree;
use crate::accession_tree::AccessionTreeNode::{Accession, Branch};
use crate::binomial::Binomial;
use rug::Float;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

fn convert_vec_i8_to_u32(kmer: &[u8]) -> Option<u32> {
    let mut acc = 0;
    for (kmer_position, n) in kmer.iter().rev().enumerate() {
        if *n == b'A' {
            // binary is 00, don't need to add anything
            continue;
        } else if *n == b'C' {
            // binary is 01
            acc += 1_u32 * 2_u32.pow(2_u32 * kmer_position as u32)
        } else if *n == b'G' {
            // binary is 10
            acc += 2_u32 * 2_u32.pow(2_u32 * kmer_position as u32)
        } else if *n == b'T' {
            // binary is 11
            acc += 3_u32 * 2_u32.pow(2_u32 * kmer_position as u32)
        } else {
            return None;
        }
    }
    Some(acc)
}

#[derive(Serialize, Deserialize)]
pub struct Database {
    index2probability: HashMap<i32, f64>,
    kmer_len: usize,
    accessions: AccessionTree,
    point2occ: Vec<i32>,
}

impl Database {
    pub fn new(kmer_len: usize) -> Self {
        Database {
            index2probability: HashMap::new(),
            accessions: AccessionTree::new(),
            kmer_len,
            point2occ: vec![-1; 4_usize.pow(kmer_len as u32)],
        }
    }

    pub fn insert_record(
        &mut self,
        forward_seq: String,
        reverse_seq: String,
        accession: String,
    ) -> () {
        let accession_index = self.accessions.push_new_node(Accession(accession));

        let mut kmer_set = HashSet::new();
        let forward_bytes = forward_seq.as_bytes();
        let reverse_bytes = reverse_seq.as_bytes();

        for (kmer_1, kmer_2) in forward_bytes
            .windows(self.kmer_len)
            .zip(reverse_bytes.windows(self.kmer_len))
        {
            let (kmer_1, kmer_2) = (convert_vec_i8_to_u32(kmer_1), convert_vec_i8_to_u32(kmer_2));
            match (kmer_1, kmer_2) {
                (Some(int_1), Some(int_2)) => {
                    kmer_set.insert(int_1);
                    kmer_set.insert(int_2);
                }
                (Some(int), None) => {
                    kmer_set.insert(int);
                }
                (None, Some(int)) => {
                    kmer_set.insert(int);
                }
                (None, None) => continue,
            }
        }

        let insert_count = kmer_set.len();
        self.insert_kmers(kmer_set, accession_index);

        self.index2probability
            .insert(accession_index, self.calculate_probability(insert_count));
    }

    pub fn query_read(&self, read: String) -> Option<&str> {
        let mut kmer_set = HashSet::new();
        let read = read.as_bytes();
        for kmer in read.windows(self.kmer_len) {
            match convert_vec_i8_to_u32(kmer) {
                None => continue,
                Some(int) => {
                    kmer_set.insert(int);
                }
            }
        }
        let num_queries = kmer_set.len() as u64;

        let mut index_to_hit_counts = HashMap::new();
        for kmer in kmer_set {
            match self.point2occ.get(kmer as usize).unwrap() {
                -1 => continue,
                index => match index_to_hit_counts.get_mut(index) {
                    None => {
                        index_to_hit_counts.insert(*index, 1_usize);
                    }
                    Some(count) => *count += 1,
                },
            }
        }
        let mut accession2hit_counts = HashMap::new();
        for (index, index_count) in index_to_hit_counts {
            let mut accessions = self.accessions.get_all_accession_indices(index).into_iter();
            while let Some(accession_index) = accessions.next() {
                match accession2hit_counts.get_mut(&accession_index) {
                    None => {
                        accession2hit_counts.insert(accession_index, index_count);
                    }
                    Some(count) => *count += index_count,
                }
            }
        }
        self.calculate_result(accession2hit_counts, num_queries)
    }

    fn calculate_result(
        &self,
        index_to_hit_counts: HashMap<i32, usize>,
        num_queries: u64,
    ) -> Option<&str> {
        let needed_probability = Float::with_val(256, 1.0e-100);
        let mut best_prob = Float::with_val(256, 1.0);
        let mut best_prob_index = None;
        for (accession_index, num_hits) in index_to_hit_counts {
            let accession_probability = self.get_probability_of_index(accession_index);
            if (num_hits as f64) < accession_probability * num_queries as f64 {
                continue;
            }
            let binomial = Binomial::new(
                Float::with_val(256, self.get_probability_of_index(accession_index)),
                num_queries,
            )
            .unwrap();
            let prob = binomial.sf(num_hits as u64).abs();
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

    fn insert_kmers(&mut self, kmers: HashSet<u32>, accession_index: i32) -> () {
        let mut previous_mappings = HashMap::new();
        for kmer in kmers {
            let current_index = self.point2occ.get_mut(kmer as usize).unwrap();
            match current_index {
                -1 => {
                    *current_index = accession_index;
                }
                index => {
                    if previous_mappings.contains_key(index) {
                        *index = *previous_mappings.get(index).unwrap()
                    } else {
                        let new_index = self
                            .accessions
                            .push_new_node(Branch(index.clone(), accession_index));
                        previous_mappings.insert(index.clone(), new_index);
                        *index = new_index
                    }
                }
            }
        }
    }

    fn calculate_probability(&self, count: usize) -> f64 {
        let prob = count as f64 / 4_usize.pow(self.kmer_len as u32) as f64;
        prob
    }

    fn get_probability_of_index(&self, accession_index: i32) -> f64 {
        *self
            .index2probability
            .get(&accession_index)
            .expect(&*format!(
                "accession index {} does not have a probability",
                accession_index
            ))
    }

    fn get_accession_of_index(&self, accession_index: i32) -> &str {
        self.accessions.get_accession_of_index(accession_index)
    }
}
