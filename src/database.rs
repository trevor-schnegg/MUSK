use crate::utility::{reverse_complement, vec_dna_bytes_to_u32};
use bio::io::fasta::Record;
use bit_iter::BitIter;
use serde::{Deserialize, Serialize};
use statrs::distribution::{Binomial, DiscreteCDF};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;
use std::path::Path;
use crate::hit_counter::HitCounter;
use crate::generator::create_random_read;

#[derive(Serialize, Deserialize)]
pub struct Database<T> {
    start_probabilities: Vec<f64>,
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
            start_probabilities: Vec::new(),
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

        self.start_probabilities
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
        self.classify_hits(
            self.collect_hits(kmers),
            num_queries as u64,
            required_probability_exponent,
        )
    }

    pub fn update_probabilities_empirically(&mut self, read_len: usize) -> () {
        let random_read = create_random_read(read_len);
        let kmers = random_read
            .as_bytes()
            .windows(self.kmer_len)
            .map(|kmer_bytes| vec_dna_bytes_to_u32(kmer_bytes).unwrap())
            .collect::<Vec<u32>>();
        let mut index_to_hit_counts = vec![HitCounter::new();9];
        for kmer in kmers.into_iter() {
            let set_bits = {
                let mut set = HashSet::new();
                for bit in BitIter::from(*self.kmer2occ.get(kmer as usize).unwrap()) {
                    set.insert(bit);
                }
                set
            };
            for i in 0_usize..9_usize {
                if set_bits.contains(&i) {
                    index_to_hit_counts.get_mut(i).unwrap().hit();
                } else {
                    index_to_hit_counts.get_mut(i).unwrap().miss();
                }
            }
        }
        for (index, (m2m, m2nm, nm2m, nm2nm)) in index_to_hit_counts.into_iter().map(|x| x.get_counts()).enumerate() {
            println!("match2match: {}, match2nonmatch: {}, nonmatch2match: {}, nonmatch2nonmatch: {}", m2m, m2nm, nm2m, nm2nm);
            let (orig_start_prob, orig_extend_prob) = self.get_probs_of_index(index);
            let observed_start_prob = nm2m as f64 / (nm2m + nm2nm) as f64;
            let observed_extend_prob = m2m as f64 / (m2m + m2nm) as f64;
            println!("start probability was originally: {}, observed start probability was: {}", orig_start_prob, observed_start_prob);
            println!("extend probability was originally: {:?}, observed extend probability was: {}\n", orig_extend_prob, observed_extend_prob);
            *self.start_probabilities.get_mut(index).unwrap() = observed_start_prob;
            match &mut self.extend_probabilities {
                Some(vec) => {
                    *vec.get_mut(index).unwrap() = observed_extend_prob;
                },
                none=> {
                    *none = Some(vec![0.25; self.start_probabilities.len()]);
                    *none.as_mut().unwrap().get_mut(index).unwrap() = observed_extend_prob;
                }
            }
        }
    }

    pub fn quick_test(&self, read_len: usize) -> () {
        let random_read = create_random_read(read_len);
        let kmers = random_read
            .as_bytes()
            .windows(self.kmer_len)
            .map(|kmer_bytes| vec_dna_bytes_to_u32(kmer_bytes).unwrap())
            .collect::<Vec<u32>>();
        let total_kmers = kmers.len() as f64;
        let mut hit_counter = vec![0_usize;9];
        for kmer in kmers {
            for bit in BitIter::from(*self.kmer2occ.get(kmer as usize).unwrap()) {
                *hit_counter.get_mut(bit).unwrap() += 1;
            }
        }
        println!("{:?}", hit_counter.into_iter().map(|x| x as f64 / total_kmers).collect::<Vec<f64>>());
    }

    pub fn completely_random_kmer_test(&self, num_kmers: usize) -> () {
        let mut hit_counts = vec![0_usize;9];
        let mut curr_num = 0_usize;
        while curr_num < num_kmers {
            let kmer = vec_dna_bytes_to_u32(create_random_read(self.kmer_len).as_bytes()).unwrap();
            for bit_index in BitIter::from(*self.kmer2occ.get(kmer as usize).unwrap()) {
                *hit_counts.get_mut(bit_index).unwrap() += 1;
            }
            curr_num += 1;
        }
        println!("{:?}", hit_counts.into_iter().map(|x| x as f64 / curr_num as f64).collect::<Vec<f64>>())
    }

    pub fn expected_hit_percentages(&self) -> Vec<f64> {
        let mut collection = Vec::new();
        for i in 0_usize..self.accessions.len() {
            collection.push(*self.start_probabilities.get(i).unwrap())
        }
        collection
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

    fn classify_hits(
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
            .start_probabilities
            .get(accession_index)
            .unwrap();
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
