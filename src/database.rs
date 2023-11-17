use crate::utility::{reverse_complement, vec_dna_bytes_to_u32};
use bio::io::fasta::Record;
use bit_iter::BitIter;
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
    num_kmers: u64,
    max_num_queries: usize,
}

impl Database<u16> {
    pub fn new(kmer_len: usize, num_queries: usize) -> Self {
        Database {
            probabilities: Vec::new(),
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

        self.probabilities
            .push(insert_count as f64 / self.num_kmers as f64);
    }

    pub fn classify_read(
        &self,
        read: Record,
        num_queries: Option<usize>,
        required_probability_exponent: Option<i32>,
    ) -> (Option<&str>, f64) {
        let max_num_queries = match num_queries {
            None => self.max_num_queries,
            Some(n) => n,
        };
        let mut collected_hits = vec![0_u64;self.accessions.len()];
        let kmers = read
            .seq()
            .windows(self.kmer_len)
            .filter_map(|kmer_bytes| vec_dna_bytes_to_u32(kmer_bytes))
            .take(max_num_queries)
            .collect::<Vec<u32>>();
        let num_queries = kmers.len();
        for kmer in kmers {
            for bit_index in BitIter::from(*self.kmer2occ.get(kmer as usize).unwrap()) {
                *collected_hits.get_mut(bit_index).unwrap() += 1;
            }
        }
        self.classify_hits(
            collected_hits,
            num_queries as u64,
            required_probability_exponent,
        )
    }

    fn classify_hits(
        &self,
        index_to_hit_counts: Vec<u64>,
        num_queries: u64,
        required_probability_exponent: Option<i32>,
    ) -> (Option<&str>, f64) {
        let needed_probability = {
            match required_probability_exponent {
                None => 1e-3,
                Some(exp) => 10.0_f64.powi(exp),
            }
        };
        let (mut best_prob, mut best_prob_index) = (1.0, None);
        for (accession_index, num_hits) in index_to_hit_counts.into_iter().enumerate() {
            let accession_probability = *self.probabilities.get(accession_index).unwrap();
            let prob = Binomial::new(accession_probability, num_queries).unwrap().sf(num_hits);
            if prob < best_prob {
                best_prob = prob;
                best_prob_index = Some(accession_index);
            }
        }
        match best_prob_index {
            None => (None, best_prob),
            Some(index) => {
                if best_prob < needed_probability {
                    (Some(self.get_accession_of_index(index)), best_prob)
                } else {
                    (None, best_prob)
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

    fn get_accession_of_index(&self, accession_index: usize) -> &str {
        self.accessions.get(accession_index).unwrap()
    }
}
