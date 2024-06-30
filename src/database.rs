use indicatif::ProgressIterator;
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
    rle::{NaiveRunLengthEncoding, RunLengthEncoding},
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

    pub fn lossy_compression(&mut self, _compression_level: usize) -> () {
        // Modify the kmer_runs variable in any desired way

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

        // Collect the hits from the read
        let mut query_count = 0;
        for kmer in KmerIter::from(read, self.kmer_len, self.canonical) {
            query_count += 1;
            for sequence in self.kmer_runs[kmer].iter() {
                collected_hits[sequence] += 1;
            }
        }

        // Classify the hits
        // Would do this using min_by_key but the Ord trait is difficult to implement for float types
        let (mut lowest_prob_index, mut lowest_prob) = (0, BigExpFloat::one());
        for (index, probability) in
            collected_hits
                .into_iter()
                .enumerate()
                .filter_map(|(index, hit_count)| {
                    let hit_count = ((hit_count as f64 / query_count as f64)
                        * self.n_queries as f64)
                        .round() as u64;

                    // Only compute if the number of hits is more than significant
                    if hit_count > self.significant_hits[index] {
                        let p = self.p_values[index];
                        let n = query_count;
                        let x = hit_count;
                        // Perform the computation using f64
                        let prob = Binomial::new(p, n).unwrap().sf(x);
                        // If the probability is greater than 0.0, use it
                        let big_exp_float_prob = if prob > 0.0 {
                            BigExpFloat::from_f64(prob)
                        } else {
                            // Otherwise, compute the probability using higher precision
                            sf(p, n, x, &self.consts)
                        };
                        Some((index, big_exp_float_prob))
                    } else {
                        // If there were 0 hits, don't compute
                        None
                    }
                })
        {
            // For each index that we computed, compare to find the lowest probability
            // If, for whatever reason, the probabilities are the same, this will use the first one
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
