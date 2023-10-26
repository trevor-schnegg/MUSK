use crate::utility::Sequence::Double;
use crate::utility::{get_kmers_as_u32, vec_dna_bytes_to_u32};
use bit_iter::BitIter;
use serde::{Deserialize, Serialize};
use statrs::distribution::{DiscreteCDF, Hypergeometric};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;
use std::path::Path;

#[derive(Serialize, Deserialize)]
pub struct Database<T> {
    index2kmer_count: Vec<u64>,
    kmer_len: usize,
    accessions: Vec<String>,
    point2occ: Vec<T>,
    num_kmers: u64,
}

impl Database<u32> {
    pub fn new(kmer_len: usize) -> Self {
        Database {
            index2kmer_count: Vec::new(),
            accessions: Vec::new(),
            kmer_len,
            point2occ: vec![0_u32; 4_usize.pow(kmer_len as u32)],
            num_kmers: 4_usize.pow(kmer_len as u32) as u64,
        }
    }

    pub fn load(file: &Path) -> Self {
        let mut f = File::open(file).expect(&*format!("could not open database file at {:?}", file));
        let mut buf: Vec<u8> = vec![];
        f.read_to_end(&mut buf).unwrap();
        bincode::deserialize(&*buf).expect(&*format!(
            "could not deserialize database file at {:?}",
            file
        ))
    }

    pub fn insert_record(
        &mut self,
        forward_seq: String,
        reverse_seq: String,
        accession: String,
    ) -> () {
        let accession_index = {
            let len = self.accessions.len();
            self.accessions.push(accession);
            len
        };
        let record_kmer_set = get_kmers_as_u32(Double(forward_seq, reverse_seq), self.kmer_len);
        let insert_count = record_kmer_set.len();

        self.insert_kmers(record_kmer_set, accession_index);

        self.index2kmer_count.push(insert_count as u64);
    }

    pub fn query_read(
        &self,
        read: String,
        num_queries: usize,
        required_probability_exponent: i32,
    ) -> Option<&str> {
        let mut already_queried = HashSet::new();
        let mut index_to_hit_counts = HashMap::new();
        let mut kmer_iter = read
            .as_bytes()
            .windows(self.kmer_len)
            .filter_map(|kmer_bytes| vec_dna_bytes_to_u32(kmer_bytes));
        while let Some(kmer) = kmer_iter.next() {
            if already_queried.len() >= num_queries {
                break;
            } else if already_queried.contains(&kmer) {
                continue;
            }
            already_queried.insert(kmer);
            for index_of_one in BitIter::from(*self.point2occ.get(kmer as usize).unwrap()) {
                match index_to_hit_counts.get_mut(&index_of_one) {
                    None => {
                        index_to_hit_counts.insert(index_of_one, 1_u64);
                    }
                    Some(count) => {
                        *count += 1;
                    }
                }
            }
        }
        self.calculate_result(
            index_to_hit_counts,
            already_queried.len() as u64,
            required_probability_exponent,
        )
    }

    fn calculate_result(
        &self,
        index_to_hit_counts: HashMap<usize, u64>,
        num_queries: u64,
        required_probability_exponent: i32,
    ) -> Option<&str> {
        let needed_probability = 10.0_f64.powi(required_probability_exponent);
        let mut best_prob = 1.0;
        let mut best_prob_index = None;
        for (accession_index, num_hits) in index_to_hit_counts {
            let num_accession_kmers = self.get_num_kmers_of_index(accession_index);
            let prob = Hypergeometric::new(self.num_kmers, num_accession_kmers, num_queries)
                .unwrap()
                .sf(num_hits);
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

    fn insert_kmers(&mut self, kmers: HashSet<u32>, accession_index: usize) -> () {
        for kmer in kmers {
            *self.point2occ.get_mut(kmer as usize).unwrap() |= 1_u32 << accession_index
        }
    }

    fn get_num_kmers_of_index(&self, accession_index: usize) -> u64 {
        *self.index2kmer_count.get(accession_index).expect(&*format!(
            "accession index {} does not have a probability",
            accession_index
        ))
    }

    fn get_accession_of_index(&self, accession_index: usize) -> &str {
        self.accessions.get(accession_index).unwrap()
    }
}
