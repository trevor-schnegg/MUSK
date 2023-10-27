use crate::utility::{reverse_complement, vec_dna_bytes_to_u32};
use bio::io::fasta::Record;
use bit_iter::BitIter;
use serde::{Deserialize, Serialize};
use statrs::distribution::{DiscreteCDF, Hypergeometric};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;
use std::path::Path;

#[derive(Serialize, Deserialize)]
pub struct Database<T> {
    index2num_kmers: Vec<u64>,
    index2required_hits: Vec<u64>,
    kmer_len: usize,
    accessions: Vec<String>,
    point2occ: Vec<T>,
    num_kmers: u64,
    max_num_queries: usize,
}

impl Database<u32> {
    pub fn new(kmer_len: usize, num_queries: usize) -> Self {
        Database {
            index2num_kmers: Vec::new(),
            index2required_hits: Vec::new(),
            accessions: Vec::new(),
            kmer_len,
            point2occ: vec![0_u32; 4_usize.pow(kmer_len as u32)],
            num_kmers: 4_usize.pow(kmer_len as u32) as u64,
            max_num_queries: num_queries,
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

        let record_kmer_set = self.get_kmers_as_u32(record.seq(), false);
        let insert_count = record_kmer_set.len();

        self.insert_kmers(record_kmer_set, accession_index);

        self.index2num_kmers.push(insert_count as u64);
        self.index2required_hits.push(self.get_needed_num_queries(insert_count as u64));
    }

    pub fn query_read(
        &self,
        read: Record,
        num_queries: Option<usize>,
        required_probability_exponent: Option<i32>,
    ) -> Option<&str> {
        let max_num_queries = match num_queries {
            None => {self.max_num_queries}
            Some(n) => {n}
        };
        let mut already_queried = HashSet::new();
        let mut index_to_hit_counts = HashMap::new();
        let mut kmer_iter = read
            .seq()
            .windows(self.kmer_len)
            .filter_map(|kmer_bytes| vec_dna_bytes_to_u32(kmer_bytes));
        while let Some(kmer) = kmer_iter.next() {
            if already_queried.len() >= max_num_queries {
                break;
            } else if already_queried.contains(&kmer) {
                continue;
            }
            already_queried.insert(kmer);
            for index in BitIter::from(*self.point2occ.get(kmer as usize).unwrap()) {
                match index_to_hit_counts.get_mut(&index) {
                    None => {
                        index_to_hit_counts.insert(index, 1_u64);
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
        required_probability_exponent: Option<i32>,
    ) -> Option<&str> {
        let needed_probability = { match required_probability_exponent {
            None => {1e-12}
            Some(exp) => {10.0_f64.powi(exp)}
        }};
        let mut best_prob = 1.0;
        let mut best_prob_index = None;
        for (accession_index, num_hits) in index_to_hit_counts {
            if num_queries == self.max_num_queries as u64 && num_hits < *self.index2required_hits.get(accession_index).unwrap() {
                continue;
            }
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

    fn get_needed_num_queries(&self, successes: u64) -> u64 {
        let population = self.num_kmers;
        let trials = self.max_num_queries;
        let hypergeometric = Hypergeometric::new(population, successes, trials as u64).unwrap();
        hypergeometric.inverse_cdf(1.0 - 1e-6)
    }

    fn get_kmers_as_u32(&self, sequence: &[u8], do_reverse_compliment: bool) -> HashSet<u32> {
        let mut kmer_set = HashSet::new();
        for kmer in sequence
            .windows(self.kmer_len)
            .filter_map(|kmer_bytes| vec_dna_bytes_to_u32(kmer_bytes))
        {
            kmer_set.insert(kmer);
        }
        if do_reverse_compliment {
            let reverse_compliment = reverse_complement(sequence);
            for kmer in reverse_compliment
                .windows(self.kmer_len)
                .filter_map(|kmer_bytes| vec_dna_bytes_to_u32(kmer_bytes))
            {
                kmer_set.insert(kmer);
            }
        }
        kmer_set
    }

    fn insert_kmers(&mut self, kmers: HashSet<u32>, accession_index: usize) -> () {
        let bit_to_set = 1_u32 << accession_index;
        for kmer in kmers {
            *self.point2occ.get_mut(kmer as usize).unwrap() |= bit_to_set
        }
    }

    fn get_num_kmers_of_index(&self, accession_index: usize) -> u64 {
        *self.index2num_kmers.get(accession_index).expect(&*format!(
            "accession index {} does not have a probability",
            accession_index
        ))
    }

    fn get_accession_of_index(&self, accession_index: usize) -> &str {
        self.accessions.get(accession_index).unwrap()
    }
}
