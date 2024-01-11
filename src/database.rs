use crate::big_exp_float::BigExpFloat;
use crate::binomial_sf::sf;
use crate::consts::Consts;
use crate::kmer_iter::KmerIter;
use crate::utility::reverse_complement;
use bio::io::fasta::Record;
use bit_iter::BitIter;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use statrs::distribution::{Binomial, DiscreteCDF};
use std::collections::HashSet;
use std::fs::File;
use std::io::Read;
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
        let mut f =
            File::open(file).expect(&*format!("could not open database file at {:?}", file));
        let mut buf: Vec<u8> = vec![];
        f.read_to_end(&mut buf).unwrap();
        bincode::deserialize(&*buf).expect(&*format!(
            "could not deserialize database file at {:?}",
            file
        ))
    }

    pub fn insert_record(&mut self, record: Record) -> () {
        let accession_index = self.accessions.len();
        self.accessions.push(record.id().to_string());

        let all_kmers = KmerIter::from(record.seq(), self.kmer_len);
        let kmer_set = HashSet::from_iter(all_kmers);
        let insert_count = kmer_set.len();

        self.insert_kmers(kmer_set, accession_index);

        self.probabilities
            .push(insert_count as f64 / self.num_kmers as f64);
    }

    pub fn classify_read(
        &self,
        read: Record,
        num_queries: usize,
    ) -> Option<(BigExpFloat, &str)> {
        let f_kmers = KmerIter::from(read.seq(), self.kmer_len)
            .take(num_queries)
            .collect::<Vec<usize>>();
        let num_queries = f_kmers.len();
        let rc = reverse_complement(read.seq());
        let rc_kmers = KmerIter::from(&rc, self.kmer_len)
            .take(num_queries)
            .collect::<Vec<usize>>();
        assert_eq!(f_kmers.len(), rc_kmers.len());

        let mut f_hits = vec![0_u64; self.accessions.len()];
        let mut rc_hits = vec![0_u64; self.accessions.len()];
        for (f_kmer, rc_kmer) in f_kmers.into_iter().zip(rc_kmers) {
            for bit_index in BitIter::from(*self.kmer2occ.get(f_kmer).unwrap()) {
                *f_hits.get_mut(bit_index).unwrap() += 1;
            }
            for bit_index in BitIter::from(*self.kmer2occ.get(rc_kmer).unwrap()) {
                *rc_hits.get_mut(bit_index).unwrap() += 1;
            }
        }
        match (self.get_lowest_probability(f_hits, num_queries as u64), self.get_lowest_probability(rc_hits, num_queries as u64)) {
            (Some(t1), Some(t2)) => {
                if t1.0 < t2.0 {
                    Some(t1)
                } else {
                    Some(t2)
                }
            },
            _ => None
        }
    }

    fn get_lowest_probability(
        &self,
        index_to_hit_counts: Vec<u64>,
        num_queries: u64,
    ) -> Option<(BigExpFloat, &str)> {
        let (mut lowest_prob, mut best_prob_index) = (BigExpFloat::one(), None);
        for (accession_index, num_hits) in index_to_hit_counts.into_iter().enumerate() {
            let accession_probability = *self.probabilities.get(accession_index).unwrap();
            if num_hits.is_zero() || (num_hits as f64) < (accession_probability * num_queries as f64) {
                continue;
            }
            let prob_f64 = Binomial::new(accession_probability, num_queries)
                .unwrap()
                .sf(num_hits);
            let prob = if prob_f64.is_zero() {
                sf(accession_probability, num_queries, num_hits, &self.consts)
            } else {
                BigExpFloat::from_f64(prob_f64)
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
