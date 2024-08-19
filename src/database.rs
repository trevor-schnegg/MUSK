use indicatif::ProgressIterator;
use itertools::Itertools;
use num_traits::One;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use serde::{Deserialize, Serialize};
use statrs::distribution::{Binomial, DiscreteCDF};
use tracing::{debug, info};

use crate::{
    big_exp_float::BigExpFloat,
    binomial_sf::sf,
    consts::Consts,
    kmer_iter::KmerIter,
    rle::{NaiveRunLengthEncoding, Run, RunLengthEncoding},
};

#[derive(Serialize, Deserialize)]
pub struct Database {
    canonical: bool,
    consts: Consts,
    cutoff_threshold: BigExpFloat,
    file2taxid: Vec<(String, usize)>,
    kmer_len: usize,
    kmer_runs: Vec<RunLengthEncoding>,
    n_queries: u64,
    p_values: Vec<f64>,
    significant_hits: Vec<u64>,
}

impl Database {
    pub fn from(
        bitmaps: Vec<RoaringBitmap>,
        canonical: bool,
        cutoff_threshold: f64,
        file2taxid: Vec<(String, usize)>,
        kmer_len: usize,
        num_queries: u64,
    ) -> Self {
        let total_num_kmers = 4_usize.pow(kmer_len as u32);

        let p_values = bitmaps
            .par_iter()
            .map(|bitmap| bitmap.len() as f64 / total_num_kmers as f64)
            .collect::<Vec<f64>>();

        let significant_hits = p_values
            .par_iter()
            .map(|p| {
                Binomial::new(*p, num_queries)
                    .unwrap()
                    .inverse_cdf(cutoff_threshold)
            })
            .collect::<Vec<u64>>();

        let mut naive_kmer_runs = vec![NaiveRunLengthEncoding::new(); total_num_kmers];
        info!("constructing naive runs...");
        // Insert all indicies into the kmer runs
        for (index, bitmap) in bitmaps.into_iter().progress().enumerate() {
            for kmer in bitmap {
                naive_kmer_runs[kmer as usize].push(index);
            }
        }

        // Log information about the number of naive runs
        let naive_run_num = naive_kmer_runs
            .iter()
            .map(|build_rle| build_rle.get_raw_runs().len())
            .sum::<usize>();
        debug!("number of naive rle runs: {}", naive_run_num);

        info!("naive runs constructed! compressing...");
        // Compress the database using uncompressed bit sets
        let kmer_runs = naive_kmer_runs
            .into_par_iter()
            .map(|build_rle| build_rle.to_rle())
            .collect::<Vec<RunLengthEncoding>>();

        // Log information about the number of compressed runs
        let compressed_run_num = kmer_runs
            .iter()
            .map(|rle| rle.get_raw_runs().len())
            .sum::<usize>();
        debug!("number of compressed rle runs: {}", compressed_run_num);

        Database {
            canonical,
            consts: Consts::new(),
            cutoff_threshold: BigExpFloat::from_f64(cutoff_threshold),
            file2taxid,
            kmer_len,
            kmer_runs,
            n_queries: num_queries,
            p_values,
            significant_hits,
        }
    }

    pub fn lossy_compression(&mut self, compression_level: usize) -> () {
        let total_set_bits_before = self.kmer_runs.par_iter().map(|runs| runs.iter().count()).sum::<usize>();
        info!("total set bits before compression {}", total_set_bits_before);

        let mut compressed_encoding: Vec<RunLengthEncoding> = vec![];

        fn can_compress(comp_level: usize, set_bits: usize, comp_gain: usize) -> bool {
            if comp_gain < 1 {
                false
            } else {
                match comp_level {
                    1 => set_bits == 1 && comp_gain == 2,
                    2 => set_bits == 1 || set_bits == 2 && comp_gain == 2,
                    3 => set_bits == 1 || set_bits == 2 || set_bits <= 4 && comp_gain == 2,
                    _ => false,
                }
            }
        }

        for encoding in self.kmer_runs.iter() {
            let runs = encoding
                .get_raw_runs()
                .into_iter()
                .map(|run| Run::from_u16(*run))
                .collect_vec();

            let mut compressed_runs: Vec<u16> = vec![];
            let mut skip_iter = false;

            let max_encoded_value = (1 << 14) - 1;

            for (i, curr) in runs.iter().enumerate() {
                if skip_iter {
                    skip_iter = false;
                    continue;
                }

                match curr {
                    Run::Uncompressed(bits) => {
                        let set_bits = bits.count_ones() as usize;
                        let mut tot_merged_bits: u16 = 15;

                        // determine whether neighbors can incorporate this unencoded block
                        let merge_prev = {
                            if i > 0 {
                                match Run::from_u16(compressed_runs[compressed_runs.len() - 1]) {
                                    Run::Zeros(count) => {
                                        if tot_merged_bits + count <= max_encoded_value {
                                            tot_merged_bits += count;
                                            true
                                        } else {
                                            false
                                        }
                                    }
                                    _ => false,
                                }
                            } else {
                                false
                            }
                        };

                        let merge_next = {
                            if i < runs.len() - 1 {
                                match runs[i + 1] {
                                    Run::Zeros(count) => {
                                        if tot_merged_bits + count <= max_encoded_value {
                                            tot_merged_bits += count;
                                            true
                                        } else {
                                            false
                                        }
                                    }
                                    _ => false,
                                }
                            } else {
                                false
                            }
                        };

                        let comp_gain: usize = match (merge_prev, merge_next) {
                            (false, false) => 0,
                            (true, false) | (false, true) => 1,
                            (true, true) => 2,
                        };

                        // test whether merge should be completed
                        if can_compress(compression_level, set_bits, comp_gain) {
                            skip_iter = merge_next;
                            if merge_prev {
                                let target_index = compressed_runs.len() - 1;
                                compressed_runs[target_index] = tot_merged_bits;
                            } else if merge_next {
                                compressed_runs.push(tot_merged_bits);
                            } else {
                                panic!("compression attempted when not possible");
                            }
                        } else {
                            compressed_runs.push(Run::to_u16(curr));
                        }
                    }
                    _ => {
                        compressed_runs.push(Run::to_u16(curr));
                    }
                }
            }

            compressed_encoding.push(RunLengthEncoding::new(compressed_runs));
        }

        // update kmer_runs
        self.kmer_runs = compressed_encoding;

        let total_set_bits_after = self.kmer_runs.par_iter().map(|runs| runs.iter().count()).sum::<usize>();
        info!("total set bits after compression {}", total_set_bits_after);

        // Recompute the p_values and significant hits after
        self.recompute_statistics();
    }

    fn recompute_statistics(&mut self) -> () {
        let total_num_kmers = 4_usize.pow(self.kmer_len as u32) as f64;
        let cutoff_threshold = self.cutoff_threshold.as_f64();

        let mut file2kmer_num = vec![0_usize; self.file2taxid.len()];

        for kmer_runs in &self.kmer_runs {
            for kmer in kmer_runs.iter() {
                file2kmer_num[kmer] += 1;
            }
        }

        let p_values = file2kmer_num
            .into_par_iter()
            .map(|kmer_num| kmer_num as f64 / total_num_kmers)
            .collect::<Vec<f64>>();

        let significant_hits = p_values
            .par_iter()
            .map(|p| {
                Binomial::new(*p, self.n_queries)
                    .unwrap()
                    .inverse_cdf(cutoff_threshold)
            })
            .collect::<Vec<u64>>();

        self.p_values = p_values;
        self.significant_hits = significant_hits;
    }

    pub fn classify(&self, read: &[u8], cutoff_threshold: BigExpFloat) -> usize {
        let mut collected_hits = vec![0_u64; self.file2taxid.len()];

        // Find the hits for all kmers
        let mut max_kmer_index = 0;
        for (index, kmer) in KmerIter::from(read, self.kmer_len, self.canonical).enumerate() {
            for sequence in self.kmer_runs[kmer].iter() {
                collected_hits[sequence] += 1;
                max_kmer_index = index;
            }
        }

        let n_total = (max_kmer_index + 1) as f64;

        // Classify the hits
        // Would do this using min_by_key but the Ord trait is difficult to implement for float types
        let (mut lowest_prob_index, mut lowest_prob) = (0, BigExpFloat::one());
        for (index, probability) in collected_hits
            .into_iter()
            .zip(self.p_values.iter())
            .enumerate()
            .filter_map(|(index, (n_hits, p))| {
                let x = (n_hits as f64 / n_total) * self.n_queries as f64;
                // Only compute if the number of hits is more than significant
                if x as f64 > (self.n_queries as f64 * p) {
                    // Perform the computation using f64
                    let prob_f64 = Binomial::new(*p, self.n_queries)
                        .unwrap()
                        .sf(x.round() as u64);
                    // If the probability is greater than 0.0, use it
                    let prob_big_exp = if prob_f64 > 0.0 {
                        BigExpFloat::from_f64(prob_f64)
                    } else {
                        // Otherwise, compute the probability using big exp
                        sf(*p, self.n_queries, x.round() as u64, &self.consts)
                    };
                    Some((index, prob_big_exp))
                } else {
                    // If there were less than a significant number of hits, don't compute
                    None
                }
            })
        {
            // For each index that we computed, compare to find the lowest probability
            // If, for whatever reason, two probabilities are the same, this will use the first one
            if probability < lowest_prob {
                (lowest_prob_index, lowest_prob) = (index, probability);
            }
        }

        if lowest_prob < cutoff_threshold {
            self.file2taxid[lowest_prob_index].1
        } else {
            0
        }
    }
}
