use indicatif::ProgressIterator;
use itertools::Itertools;
use num_traits::One;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use serde::{Deserialize, Serialize};
use statrs::distribution::{Binomial, DiscreteCDF};
use std::{collections::HashMap, fs::File, path::Path};
use tracing::{debug, info};

use crate::{
    big_exp_float::BigExpFloat,
    binomial_sf::sf,
    consts::Consts,
    io::{dump_data_to_file, load_data_from_file},
    kmer_iter::KmerIter,
    rle::{Block, NaiveRunLengthEncoding, RunLengthEncoding, MAX_RUN, MAX_UNCOMPRESSED_BITS},
};

#[derive(Serialize, Deserialize)]
pub struct Database {
    canonical: bool,
    consts: Consts,
    file2taxid: Vec<(String, usize)>,
    kmer_len: usize,
    kmer_rles: HashMap<u32, RunLengthEncoding>,
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

        let mut naive_kmer_rles: HashMap<u32, NaiveRunLengthEncoding> = HashMap::new();
        info!("constructing naive runs...");
        // Insert all indicies into the kmer runs
        for (index, bitmap) in bitmaps.into_iter().progress().enumerate() {
            for kmer in bitmap {
                match naive_kmer_rles.get_mut(&kmer) {
                    Some(naive_rle) => naive_rle.push(index),
                    None => {
                        let mut new_naive_rle = NaiveRunLengthEncoding::new();
                        new_naive_rle.push(index);
                        naive_kmer_rles.insert(kmer, new_naive_rle);
                    }
                }
            }
        }

        // Log information about the number of naive runs
        let naive_run_num = naive_kmer_rles
            .par_iter()
            .map(|(_kmer, naive_rle)| naive_rle.get_raw_runs().len())
            .sum::<usize>();
        debug!("number of naive rle runs: {}", naive_run_num);

        info!("naive runs constructed! compressing...");
        // Compress the database using uncompressed bit sets
        // Process in parallel and then insert into a new hashmap
        let kmer_rles_vec = naive_kmer_rles
            .into_par_iter()
            .map(|(kmer, naive_rle)| (kmer, naive_rle.to_rle()))
            .collect::<Vec<(u32, RunLengthEncoding)>>();
        let kmer_rles = HashMap::from_iter(kmer_rles_vec.into_iter());

        // Log information about the number of compressed runs
        let compressed_block_num = kmer_rles
            .par_iter()
            .map(|(_kmer, rle)| rle.get_raw_blocks().len())
            .sum::<usize>();
        debug!("number of compressed rle runs: {}", compressed_block_num);

        Database {
            canonical,
            consts: Consts::new(),
            file2taxid,
            kmer_len,
            kmer_rles,
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

        fn get_zeros(run: Option<&Block>) -> Option<usize> {
            match run {
                Some(Block::Zeros(x)) => Some(*x as usize),
                _ => None,
            }
        }

        info!(
            "total set bits before compression {}",
            self.kmer_rles
                .par_iter()
                .map(|(_kmer, rle)| rle.iter().count())
                .sum::<usize>()
        );

        self.kmer_rles
            .par_iter_mut()
            .for_each(|(_kmer, current_rle)| {
                // variable to hold the new lossy compressed blocks as u16s
                let mut compressed_blocks = vec![];

                // peekable iterator over the current runs
                let mut block_iter = current_rle
                    .get_raw_blocks()
                    .into_iter()
                    .map(|block| Block::from_u16(*block))
                    .peekable();

                while let Some(curr_block) = block_iter.next() {
                    match curr_block {
                        // if the current run is an uncompressed run, lossy compression may apply
                        Block::Uncompressed(bits) => {
                            // the following variables are `Option` types
                            // if `None`, the previous or next run either did not exist or was not a run of zeros
                            // otherwise, it will be `Some(<length of zero run>)`
                            let prev_zeros = get_zeros(compressed_blocks.last());
                            let next_zeros = get_zeros(block_iter.peek());

                            let set_bits = bits.count_ones();

                            // decide what to do based on the previous and next run
                            match (prev_zeros, next_zeros) {
                                (None, None) => {
                                    // no lossy compression possible, push and continue
                                    compressed_blocks.push(curr_block);
                                }
                                (Some(length), None) => {
                                    // try to merge with the previous
                                    let blocks_saved =
                                        if length + MAX_UNCOMPRESSED_BITS <= MAX_RUN as usize {
                                            1
                                        } else {
                                            0
                                        };
                                    if should_compress(compression_level, set_bits, blocks_saved) {
                                        match compressed_blocks.last_mut().unwrap() {
                                            // intentionally overloading the variable `length` because they are the same thing
                                            Block::Zeros(length) => {
                                                *length += MAX_UNCOMPRESSED_BITS as u16
                                            }
                                            _ => panic!("impossible case"),
                                        }
                                    } else {
                                        // if we shouldn't compress, push and continue
                                        compressed_blocks.push(curr_block);
                                    }
                                }
                                (None, Some(length)) => {
                                    // try to merge with the next
                                    let blocks_saved =
                                        if length + MAX_UNCOMPRESSED_BITS <= MAX_RUN as usize {
                                            1
                                        } else {
                                            0
                                        };
                                    if should_compress(compression_level, set_bits, blocks_saved) {
                                        // get/burn the next run in the iterator which should be zeros
                                        match block_iter.next().unwrap() {
                                            // intentionally overloading the variable `length` because they are the same thing
                                            Block::Zeros(length) => compressed_blocks.push(
                                                Block::Zeros(MAX_UNCOMPRESSED_BITS as u16 + length),
                                            ),
                                            _ => panic!("impossible case"),
                                        }
                                    } else {
                                        // if we shouldn't compress, push and continue
                                        compressed_blocks.push(curr_block);
                                    }
                                }
                                (Some(prev_length), Some(next_length)) => {
                                    // try to merge with both
                                    let total_length =
                                        prev_length + MAX_UNCOMPRESSED_BITS + next_length;
                                    let blocks_saved = if total_length <= MAX_RUN as usize {
                                        2
                                    } else if total_length <= (MAX_RUN as usize) * 2 {
                                        1
                                    } else {
                                        0
                                    };
                                    if should_compress(compression_level, set_bits, blocks_saved) {
                                        match (
                                            compressed_blocks.last_mut().unwrap(),
                                            block_iter.next().unwrap(),
                                        ) {
                                            // once again, intentionally overloading these variables because they are the same
                                            (
                                                Block::Zeros(prev_length),
                                                Block::Zeros(next_length),
                                            ) => {
                                                if total_length <= MAX_RUN as usize {
                                                    // if the total length can fit as one run, just add it to the existing run
                                                    *prev_length +=
                                                        MAX_UNCOMPRESSED_BITS as u16 + next_length;
                                                } else {
                                                    // otherwise, max out the previous run and insert the leftover as a new run of zeros
                                                    *prev_length = MAX_RUN;
                                                    compressed_blocks.push(Block::Zeros(
                                                        (total_length - MAX_RUN as usize) as u16,
                                                    ))
                                                }
                                            }
                                            _ => panic!("impossible case"),
                                        }
                                    } else {
                                        // if we shouldn't compress, push and continue
                                        compressed_blocks.push(curr_block);
                                    }
                                }
                            };
                        }
                        // if the current run isn't an uncompressed run, push it and continue
                        _ => {
                            compressed_blocks.push(curr_block);
                        }
                    } // end match
                } // end while

                // At this point, we have created the variable `compressed_blocks` with the lossy compressed blocks in it
                // We can then just update the original rle to the new lossy compressed blocks

                *current_rle = RunLengthEncoding::from(
                    compressed_blocks
                        .into_iter()
                        .map(|run| run.to_u16())
                        .collect_vec(),
                );
            }); // end mut for each

        info!(
            "total set bits after compression {}",
            self.kmer_rles
                .par_iter()
                .map(|(_kmer, rle)| rle.iter().count())
                .sum::<usize>()
        );

        // Recompute the p_values after
        self.recompute_p_values();
    }

    fn recompute_p_values(&mut self) -> () {
        let total_num_kmers = 4_usize.pow(self.kmer_len as u32) as f64;

        let mut file2kmer_num = vec![0_usize; self.file2taxid.len()];

        for (_kmer, rle) in self.kmer_rles.iter() {
            for kmer in rle.iter() {
                file2kmer_num[kmer] += 1;
            }
        }

        let p_values = file2kmer_num
            .into_par_iter()
            .map(|kmer_num| kmer_num as f64 / total_num_kmers)
            .collect::<Vec<f64>>();

        self.p_values = p_values;
    }

    pub fn classify(
        &self,
        read: &[u8],
        cutoff_threshold: BigExpFloat,
        n_max: u64,
    ) -> Option<&(String, usize)> {
        let mut collected_hits = vec![0_u64; self.file2taxid.len()];

        // Find the hits for all kmers
        let mut n_total = 0_u64;
        for kmer in KmerIter::from(read, self.kmer_len, self.canonical).map(|kmer| kmer as u32) {
            if let Some(rle) = self.kmer_rles.get(&kmer) {
                for sequence in rle.iter() {
                    collected_hits[sequence] += 1;
                }
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
            Some(&self.file2taxid[lowest_prob_index])
        } else {
            None
        }
    }

    pub fn serialize_to(mut self, file: File, metadata_file: File) -> () {
        dump_data_to_file(&self.kmer_rles, file).unwrap();
        self.kmer_rles = HashMap::new();
        dump_data_to_file(&self, metadata_file).unwrap();
    }

    pub fn deserialize_from(file: &Path, metadata_file: &Path) -> Self {
        debug!("loading metadata for database");
        let mut database = load_data_from_file::<Self>(metadata_file);
        debug!("loading kmer rles for database");
        let kmer_rles = load_data_from_file::<HashMap<u32, RunLengthEncoding>>(file);
        debug!("done loading database");
        database.kmer_rles = kmer_rles;
        database
    }
}
