use crate::big_exp_float::BigExpFloat;
use crate::binomial_sf::sf;
use crate::consts::Consts;
use crate::io::load_data_from_file;
use crate::kmer_iter::KmerIter;
use bio::io::fasta::Record;
use bit_iter::BitIter;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use statrs::distribution::{Binomial, DiscreteCDF};
use std::collections::HashSet;
use std::path::Path;

#[derive(Serialize, Deserialize)]
pub struct Database<T> {
    probabilities: Vec<f64>,
    kmer_len: usize,
    accessions: Vec<String>,
    kmer2occ: Vec<T>,
    num_kmers: usize,
    consts: Consts,
}

impl Database<u16> {
    pub fn new(kmer_len: usize) -> Self {
        Database {
            probabilities: Vec::new(),
            accessions: Vec::new(),
            kmer_len,
            kmer2occ: vec![0_u16; 4_usize.pow(kmer_len as u32)],
            num_kmers: 4_usize.pow(kmer_len as u32),
            consts: Consts::new(),
        }
    }

    pub fn load(file: &Path) -> Self {
        load_data_from_file::<Self>(file)
    }

    pub fn insert_record(&mut self, record: Record) -> () {
        let accession_index = self.accessions.len();
        self.accessions.push(record.id().to_string());

        let all_kmers = KmerIter::from(record.seq(), self.kmer_len, false);
        let kmer_set = HashSet::from_iter(all_kmers);
        let insert_count = kmer_set.len();

        self.insert_kmers(kmer_set, accession_index);

        self.probabilities
            .push(insert_count as f64 / self.num_kmers as f64);
    }

    pub fn classify_read(&self, read: Record, num_queries: usize) -> Option<(BigExpFloat, &str)> {
        let kmers = KmerIter::from(read.seq(), self.kmer_len, false)
            .take(num_queries)
            .collect::<Vec<usize>>();
        let num_queries = kmers.len();

        let mut hits = vec![0_u64; self.accessions.len()];
        for kmer in kmers {
            for bit_index in BitIter::from(*self.kmer2occ.get(kmer).unwrap()) {
                *hits.get_mut(bit_index).unwrap() += 1;
            }
        }
        self.get_lowest_probability(hits, num_queries as u64)
    }

    fn get_lowest_probability(
        &self,
        index_to_hit_counts: Vec<u64>,
        num_queries: u64,
    ) -> Option<(BigExpFloat, &str)> {
        let (mut lowest_prob, mut best_prob_index) = (BigExpFloat::one(), None);
        for (accession_index, num_hits) in index_to_hit_counts.into_iter().enumerate() {
            let accession_probability = *self.probabilities.get(accession_index).unwrap();
            if num_hits.is_zero()
                || (num_hits as f64) < (accession_probability * num_queries as f64)
            {
                continue;
            }
            let prob = {
                let prob_f64 = Binomial::new(accession_probability, num_queries)
                    .unwrap()
                    .sf(num_hits);
                if prob_f64.is_zero() {
                    sf(accession_probability, num_queries, num_hits, &self.consts)
                } else {
                    BigExpFloat::from_f64(prob_f64)
                }
            };
            if prob < lowest_prob {
                lowest_prob = prob;
                best_prob_index = Some(accession_index);
            }
        }
        match best_prob_index {
            None => None,
            Some(index) => Some((lowest_prob, self.get_accession_of_index(index))),
        }
    }

    fn insert_kmers(&mut self, kmers: HashSet<usize>, accession_index: usize) -> () {
        let bit_to_set = 1_u16 << accession_index;
        for kmer in kmers {
            *self.kmer2occ.get_mut(kmer).unwrap() |= bit_to_set
        }
    }

    fn get_accession_of_index(&self, accession_index: usize) -> &str {
        self.accessions.get(accession_index).unwrap()
    }
}
