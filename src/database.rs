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

const MAX_RUN: u16 = (1 << 14) - 1;

#[derive(Serialize, Deserialize)]
pub struct Database {
    canonical: bool,
    consts: Consts,
    cutoff_threshold: BigExpFloat,
    file2taxid: Vec<(String, usize)>,
    kmer_len: usize,
    kmer_rles: Vec<RunLengthEncoding>,
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
            kmer_rles: kmer_runs,
            n_queries: num_queries,
            p_values,
            significant_hits,
        }
    }

    pub fn lossy_compression(&mut self, compression_level: usize) -> () {
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

        info!(
            "total set bits before compression {}",
            self.kmer_rles
                .par_iter()
                .map(|runs| runs.iter().count())
                .sum::<usize>()
        );

        self.kmer_rles = self
            .kmer_rles
            .par_iter()
            .map(|runs_vec| {
                let mut compressed_runs: Vec<Run> = vec![];

                let collected_runs = runs_vec
                    .get_raw_runs()
                    .into_iter()
                    .map(|run| Run::from_u16(*run))
                    .collect_vec();

                let final_window = collected_runs.len() - 1;

                let mut runs_iter = collected_runs.windows(2).enumerate();

                while let Some((index, window)) = runs_iter.next() {
                    let curr_run = window[0];
                    let next_run = window[1];

                    match curr_run {
                        Run::Uncompressed(bits) => {
                            let set_bits = bits.count_ones() as usize;
                            let mut total_merged_bits: u16 = 15;

                            let merge_prev = {
                                if let Some(&prev_run) = compressed_runs.last() {
                                    match prev_run {
                                        Run::Zeros(count) => {
                                            if total_merged_bits + count <= MAX_RUN {
                                                total_merged_bits += count;
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
                                match next_run {
                                    Run::Zeros(count) => {
                                        if total_merged_bits + count <= MAX_RUN {
                                            total_merged_bits += count;
                                            true
                                        } else {
                                            false
                                        }
                                    }
                                    _ => false,
                                }
                            };

                            let comp_gain: usize = match (merge_prev, merge_next) {
                                (false, false) => 0,
                                (true, false) | (false, true) => 1,
                                (true, true) => 2,
                            };

                            if can_compress(compression_level, set_bits, comp_gain) {
                                if merge_prev {
                                    if let Some(last) = compressed_runs.last_mut() {
                                        *last = Run::from_u16(total_merged_bits);
                                    } else {
                                        panic!(
                                            "merge with previous attempted when no previous exists"
                                        );
                                    }
                                } else {
                                    compressed_runs.push(Run::from_u16(total_merged_bits));
                                }

                                if merge_next {
                                    runs_iter.next();
                                }
                            } else {
                                compressed_runs.push(curr_run);
                            }

                            if index >= final_window && !merge_next {
                                compressed_runs.push(next_run);
                            }
                        }
                        _ => {
                            compressed_runs.push(curr_run);

                            if index >= final_window {
                                compressed_runs.push(next_run);
                            }
                        }
                    }
                }

                RunLengthEncoding::new(
                    compressed_runs
                        .iter()
                        .map(|run| Run::to_u16(run))
                        .collect_vec(),
                )
            })
            .collect::<Vec<RunLengthEncoding>>();

        info!(
            "total set bits after compression {}",
            self.kmer_rles
                .par_iter()
                .map(|runs| runs.iter().count())
                .sum::<usize>()
        );

        // Recompute the p_values and significant hits after
        self.recompute_statistics();
    }

    fn recompute_statistics(&mut self) -> () {
        let total_num_kmers = 4_usize.pow(self.kmer_len as u32) as f64;
        let cutoff_threshold = self.cutoff_threshold.as_f64();

        let mut file2kmer_num = vec![0_usize; self.file2taxid.len()];

        for kmer_runs in &self.kmer_rles {
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

    pub fn classify(&self, read: &[u8], cutoff_threshold: BigExpFloat, full_read: bool) -> usize {
        let mut collected_hits = vec![0_u64; self.file2taxid.len()];

        // Find the hits for all kmers
        let mut max_kmer_index = 0;
        for (index, kmer) in KmerIter::from(read, self.kmer_len, self.canonical).enumerate() {
            for sequence in self.kmer_rles[kmer].iter() {
                collected_hits[sequence] += 1;
                max_kmer_index = index;
            }
        }

        let n_total = max_kmer_index + 1;

        // Classify the hits
        // Would do this using min_by_key but the Ord trait is difficult to implement for float types
        let (mut lowest_prob_index, mut lowest_prob) = (0, BigExpFloat::one());
        for (index, probability) in collected_hits
            .into_iter()
            .zip(self.p_values.iter())
            .enumerate()
            .filter_map(|(index, (n_hits, p))| {
                let x = if full_read {
                    n_hits as f64
                } else {
                    (n_hits as f64 / n_total as f64) * self.n_queries as f64
                };
                let n = if full_read {
                    n_total as u64
                } else {
                    self.n_queries
                };

                // Only compute if the number of hits is more than significant
                if x > (n as f64 * p) {
                    // Perform the computation using f64
                    let prob_f64 = Binomial::new(*p, n).unwrap().sf(x.round() as u64);
                    // If the probability is greater than 0.0, use it
                    let prob_big_exp = if prob_f64 > 0.0 {
                        BigExpFloat::from_f64(prob_f64)
                    } else {
                        // Otherwise, compute the probability using big exp
                        sf(*p, n, x.round() as u64, &self.consts)
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
