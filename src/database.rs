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
    rle::{NaiveRunLengthEncoding, Run, RunLengthEncoding, MAX_RUN, MAX_UNCOMPRESSED_BITS},
};

#[derive(Serialize, Deserialize)]
pub struct Database {
    canonical: bool,
    consts: Consts,
    file2taxid: Vec<(String, usize)>,
    kmer_len: usize,
    kmer_rles: Vec<RunLengthEncoding>,
    p_values: Vec<f64>,
}

impl Database {
    pub fn from(
        bitmaps: Vec<RoaringBitmap>,
        canonical: bool,
        file2taxid: Vec<(String, usize)>,
        kmer_len: usize,
    ) -> Self {
        let total_num_kmers = 4_usize.pow(kmer_len as u32);

        let p_values = bitmaps
            .par_iter()
            .map(|bitmap| bitmap.len() as f64 / total_num_kmers as f64)
            .collect::<Vec<f64>>();

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
            file2taxid,
            kmer_len,
            kmer_rles: kmer_runs,
            p_values,
        }
    }

    pub fn lossy_compression(&mut self, compression_level: usize) -> () {
        fn should_compress(compression_level: usize, set_bits: u32, run_reduction: usize) -> bool {
            if run_reduction < 1 {
                false
            } else {
                match compression_level {
                    1 => set_bits == 1 && run_reduction == 2,
                    2 => set_bits == 1 || set_bits == 2 && run_reduction == 2,
                    3 => set_bits == 1 || set_bits == 2 || set_bits <= 4 && run_reduction == 2,
                    _ => false,
                }
            }
        }

        fn get_zeros(run: Option<&Run>) -> Option<usize> {
            match run {
                Some(Run::Zeros(x)) => Some(*x as usize),
                _ => None,
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
            .map(|current_runs| {
                // variable to hold the new lossy compressed runs as u16s
                let mut compressed_runs = vec![];

                // peekable iterator over the current runs
                let mut runs_iter = current_runs
                    .get_raw_runs()
                    .into_iter()
                    .map(|run| Run::from_u16(*run))
                    .peekable();

                while let Some(curr_run) = runs_iter.next() {
                    match curr_run {
                        // if the current run is an uncompressed run, lossy compression may apply
                        Run::Uncompressed(bits) => {
                            // the following variables are `Option` types
                            // if `None`, the previous or next run either did not exist or was not a run of zeros
                            // otherwise, it will be `Some(<length of zero run>)`
                            let prev_zeros = get_zeros(compressed_runs.last());
                            let next_zeros = get_zeros(runs_iter.peek());

                            let set_bits = bits.count_ones();

                            // decide what to do based on the previous and next run
                            match (prev_zeros, next_zeros) {
                                (None, None) => {
                                    // no lossy compression possible, push and continue
                                    compressed_runs.push(curr_run);
                                }
                                (Some(length), None) => {
                                    // try to merge with the previous
                                    let run_reduction =
                                        if length + MAX_UNCOMPRESSED_BITS <= MAX_RUN as usize {
                                            1
                                        } else {
                                            0
                                        };
                                    if should_compress(compression_level, set_bits, run_reduction) {
                                        match compressed_runs.last_mut().unwrap() {
                                            // intentionally overloading the variable `length` because they are the same thing
                                            Run::Zeros(length) => {
                                                *length += MAX_UNCOMPRESSED_BITS as u16
                                            }
                                            _ => panic!("impossible case"),
                                        }
                                    } else {
                                        // if we shouldn't compress, push and continue
                                        compressed_runs.push(curr_run);
                                    }
                                }
                                (None, Some(length)) => {
                                    // try to merge with the next
                                    let run_reduction =
                                        if length + MAX_UNCOMPRESSED_BITS <= MAX_RUN as usize {
                                            1
                                        } else {
                                            0
                                        };
                                    if should_compress(compression_level, set_bits, run_reduction) {
                                        // get/burn the next run in the iterator which should be zeros
                                        match runs_iter.next().unwrap() {
                                            // intentionally overloading the variable `length` because they are the same thing
                                            Run::Zeros(length) => compressed_runs.push(Run::Zeros(
                                                MAX_UNCOMPRESSED_BITS as u16 + length,
                                            )),
                                            _ => panic!("impossible case"),
                                        }
                                    } else {
                                        // if we shouldn't compress, push and continue
                                        compressed_runs.push(curr_run);
                                    }
                                }
                                (Some(prev_length), Some(next_length)) => {
                                    // try to merge with both
                                    let total_length =
                                        prev_length + MAX_UNCOMPRESSED_BITS + next_length;
                                    let run_reduction = if total_length <= MAX_RUN as usize {
                                        2
                                    } else if total_length <= (MAX_RUN as usize) * 2 {
                                        1
                                    } else {
                                        0
                                    };
                                    if should_compress(compression_level, set_bits, run_reduction) {
                                        match (
                                            compressed_runs.last_mut().unwrap(),
                                            runs_iter.next().unwrap(),
                                        ) {
                                            // once again, intentionally overloading these variables because they are the same
                                            (Run::Zeros(prev_length), Run::Zeros(next_length)) => {
                                                if total_length <= MAX_RUN as usize {
                                                    // if the total length can fit as one run, just add it to the existing run
                                                    *prev_length +=
                                                        MAX_UNCOMPRESSED_BITS as u16 + next_length;
                                                } else {
                                                    // otherwise, max out the previous run and insert the leftover as a new run of zeros
                                                    *prev_length = MAX_RUN;
                                                    compressed_runs.push(Run::Zeros(
                                                        (total_length - MAX_RUN as usize) as u16,
                                                    ))
                                                }
                                            }
                                            _ => panic!("impossible case"),
                                        }
                                    } else {
                                        // if we shouldn't compress, push and continue
                                        compressed_runs.push(curr_run);
                                    }
                                }
                            };
                        }
                        // if the current run isn't an uncompressed run, push it and continue
                        _ => {
                            compressed_runs.push(curr_run);
                        }
                    } // end match
                } // end while

                RunLengthEncoding::new(
                    compressed_runs
                        .into_iter()
                        .map(|run| run.to_u16())
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

        self.p_values = p_values;
    }

    pub fn classify(&self, read: &[u8], cutoff_threshold: BigExpFloat, n_max: u64) -> usize {
        let mut collected_hits = vec![0_u64; self.file2taxid.len()];

        // Find the hits for all kmers
        let mut n_total = 0;
        for kmer in KmerIter::from(read, self.kmer_len, self.canonical) {
            for sequence in self.kmer_rles[kmer].iter() {
                collected_hits[sequence] += 1;
            }
            n_total += 1;
        }

        // Classify the hits
        // Would do this using min_by_key but the Ord trait is difficult to implement for float types
        let (mut lowest_prob_index, mut lowest_prob) = (0, BigExpFloat::one());
        for (index, probability) in collected_hits
            .into_iter()
            .zip(self.p_values.iter())
            .enumerate()
            .filter_map(|(index, (n_hits, p))| {
                let x = if n_total <= n_max {
                    n_hits
                } else {
                    ((n_hits as f64 / n_total as f64) * n_max as f64).round() as u64
                };
                let n = if n_total <= n_max { n_total } else { n_max };

                // This check saves runtime in practice
                // Only do probability computation if the classification probability is going to be < 0.5
                if x as f64 > (n as f64 * p) {
                    // Perform the computation using f64
                    let prob_f64 = Binomial::new(*p, n).unwrap().sf(x);
                    // If the probability is greater than 0.0, use it
                    let prob_big_exp = if prob_f64 > 0.0 {
                        BigExpFloat::from_f64(prob_f64)
                    } else {
                        // Otherwise, compute the probability using big exp
                        sf(*p, n, x, &self.consts)
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
