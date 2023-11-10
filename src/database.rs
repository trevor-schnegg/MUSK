use crate::utility::{reverse_complement, vec_dna_bytes_to_u32};
use bio::io::fasta::Record;
use bit_iter::BitIter;
use serde::{Deserialize, Serialize};
use statrs::distribution::{Binomial, DiscreteCDF};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;
use std::path::Path;

#[derive(Serialize, Deserialize)]
pub struct Database<T> {
    open_probabilities: Vec<f64>,
    extend_probabilities: Option<Vec<f64>>,
    kmer_len: usize,
    accessions: Vec<String>,
    kmer2occ: Vec<T>,
    num_kmers: u64,
    max_num_queries: usize,
}

impl Database<u16> {
    pub fn new(kmer_len: usize, num_queries: usize) -> Self {
        Database {
            open_probabilities: Vec::new(),
            extend_probabilities: None,
            accessions: Vec::new(),
            kmer_len,
            kmer2occ: vec![0_u16; 4_usize.pow(kmer_len as u32)],
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

        let record_kmer_set = self.get_kmers_as_u32(record.seq(), true);
        let insert_count = record_kmer_set.len();

        self.insert_kmers(record_kmer_set, accession_index);

        self.open_probabilities
            .push(insert_count as f64 / self.num_kmers as f64);
    }

    pub fn classify_read(
        &self,
        read: Record,
        num_queries: Option<usize>,
        required_probability_exponent: Option<i32>,
    ) -> Option<&str> {
        let max_num_queries = match num_queries {
            None => self.max_num_queries,
            Some(n) => n,
        };
        let kmers = read
            .seq()
            .windows(self.kmer_len)
            .filter_map(|kmer_bytes| vec_dna_bytes_to_u32(kmer_bytes))
            .take(max_num_queries)
            .collect::<Vec<u32>>();
        let num_queries = kmers.len();
        self.calculate_result(
            self.collect_hits(kmers),
            num_queries as u64,
            required_probability_exponent,
        )
    }

    fn collect_hits(&self, kmers: Vec<u32>) -> HashMap<usize, (u64, u64, usize)> {
        let mut index_to_hit_counts = HashMap::new();
        let mut kmer_iter = kmers.into_iter().enumerate();
        while let Some((idx, kmer)) = kmer_iter.next() {
            for index in BitIter::from(*self.kmer2occ.get(kmer as usize).unwrap()) {
                match index_to_hit_counts.get_mut(&index) {
                    None => {
                        index_to_hit_counts.insert(index, (1_u64, 0_u64, idx));
                    }
                    Some(counts) => {
                        if counts.2 == idx - 1 {
                            counts.1 += 1;
                        } else {
                            counts.0 += 1;
                        }
                        counts.2 = idx;
                    }
                }
            }
        }
        index_to_hit_counts
    }

    fn calculate_result(
        &self,
        index_to_hit_counts: HashMap<usize, (u64, u64, usize)>,
        num_queries: u64,
        required_probability_exponent: Option<i32>,
    ) -> Option<&str> {
        let needed_probability = {
            match required_probability_exponent {
                None => 1e-3,
                Some(exp) => 10.0_f64.powi(exp),
            }
        };
        let (mut best_prob, mut best_prob_index) = (1.0, None);
        for (accession_index, (num_starts, num_extends, _)) in index_to_hit_counts {
            let (open_prob, extend_prob) = self.get_probs_of_index(accession_index);
            let extend_prob = match extend_prob {
                None => 0.25,
                Some(prob) => prob,
            };
            let prob = {
                let prob_starts =
                    Binomial::new(open_prob, num_queries - (num_extends + num_starts))
                        .unwrap()
                        .sf(num_starts);
                let prob_extends = Binomial::new(extend_prob, num_extends + num_starts)
                    .unwrap()
                    .sf(num_extends);
                prob_starts * prob_extends
            };
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
        let bit_to_set = 1_u16 << accession_index;
        for kmer in kmers {
            *self.kmer2occ.get_mut(kmer as usize).unwrap() |= bit_to_set
        }
    }

    fn get_probs_of_index(&self, accession_index: usize) -> (f64, Option<f64>) {
        let open_prob = *self
            .open_probabilities
            .get(accession_index)
            .expect(&*format!(
                "accession index {} does not have a probability",
                accession_index
            ));
        let extend_prob = match &self.extend_probabilities {
            None => None,
            Some(map) => Some(*map.get(accession_index).unwrap()),
        };
        (open_prob, extend_prob)
    }

    fn get_accession_of_index(&self, accession_index: usize) -> &str {
        self.accessions.get(accession_index).unwrap()
    }
}
