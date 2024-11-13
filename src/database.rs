use itertools::Itertools;
use num_traits::One;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use serde::{Deserialize, Serialize};
use statrs::distribution::{Binomial, DiscreteCDF};
use std::{collections::HashMap, u16, u32};
use tracing::{debug, error, info};

use crate::{
    big_exp_float::BigExpFloat,
    binomial_sf::sf,
    consts::Consts,
    kmer_iter::KmerIter,
    rle::{Block, NaiveRunLengthEncoding, RunLengthEncodingIter, MAX_RUN, MAX_UNCOMPRESSED_BITS},
};

#[derive(Serialize, Deserialize)]
pub struct Database {
    canonical: bool,
    consts: Consts,
    files: Vec<String>,
    flat_rles: Box<[u16]>,
    tax_ids: Vec<usize>,
    kmer_len: usize,
    kmer_to_rle_range: HashMap<u32, (u32, u16)>,
    p_values: Vec<f64>,
}

fn create_flat_index(
    kmers_and_rles: Vec<(u32, Vec<u16>)>,
) -> (HashMap<u32, (u32, u16)>, Box<[u16]>) {
    // Create the true underlying data structures from previous variables
    let mut kmer_to_rle_range = HashMap::with_capacity(kmers_and_rles.len());
    let mut current_index = 0_u32;
    let flat_rles = kmers_and_rles
        .into_iter()
        .map(|(kmer, rle)| {
            // Check to make sure the block can be indexed correctly
            if rle.len() > u16::MAX as usize {
                error!("the size a single run-length encoding exceeded what the current implementation allows");
                panic!();
            }

            // Check to make sure the u32 index doesn't overflow
            if current_index as usize + rle.len() > u32::MAX as usize {
                error!("the total size of the run-length encodings exceeded what the current implementation allows");
                panic!();
            }

            // Insert into the range and return the raw u16 blocks
            kmer_to_rle_range.insert(kmer, (current_index, rle.len() as u16));
            current_index += rle.len() as u32;
            rle
        })
        .flatten()
        .collect::<Box<[u16]>>();
    (kmer_to_rle_range, flat_rles)
}

impl Database {
    pub fn from(
        bitmaps: Vec<RoaringBitmap>,
        canonical: bool,
        files: Vec<String>,
        tax_ids: Vec<usize>,
        kmer_len: usize,
    ) -> Self {
        let total_canonical_kmers =
            (4_usize.pow(kmer_len as u32) - 4_usize.pow(kmer_len.div_ceil(2) as u32)) / 2;

        let p_values = bitmaps
            .par_iter()
            .map(|bitmap| bitmap.len() as f64 / total_canonical_kmers as f64)
            .collect::<Vec<f64>>();

        // Initialize the naive RLEs to be the maximum size
        let mut kmer_to_naive_rle: HashMap<u32, NaiveRunLengthEncoding> =
            HashMap::with_capacity(total_canonical_kmers);

        // Construct all naive kmer RLEs from the bitmaps
        info!("constructing naive runs...");
        for (index, bitmap) in bitmaps.into_iter().enumerate() {
            for kmer in bitmap {
                match kmer_to_naive_rle.get_mut(&kmer) {
                    Some(naive_rle) => naive_rle.push(index),
                    None => {
                        let mut new_naive_rle = NaiveRunLengthEncoding::new();
                        new_naive_rle.push(index);
                        kmer_to_naive_rle.insert(kmer, new_naive_rle);
                    }
                }
            }
        }

        // Log information about the number of naive runs
        let naive_run_count = kmer_to_naive_rle
            .par_iter()
            .map(|(_kmer, naive_rle)| naive_rle.get_raw_runs().len())
            .sum::<usize>();
        debug!("number of naive rle runs: {}", naive_run_count);

        // Compress the database by allowing the use of uncompressed bit sets
        info!("naive runs constructed! allowing uncompressed bit sets...");
        let kmers_and_rles = kmer_to_naive_rle
            .into_par_iter()
            .map(|(kmer, naive_rle)| (kmer, naive_rle.to_rle().into_raw_blocks()))
            .collect::<Vec<(u32, Vec<u16>)>>();

        // Log information about the number of compressed runs
        let compressed_block_num = kmers_and_rles
            .par_iter()
            .map(|(_kmer, rle)| rle.len())
            .sum::<usize>();
        debug!(
            "number of rle runs after allowing uncompressed bit sets: {}",
            compressed_block_num
        );

        info!("uncompressed bit sets added! flattening database...");
        let (kmer_to_rle_range, flat_rles) = create_flat_index(kmers_and_rles);

        Database {
            canonical,
            consts: Consts::new(),
            files,
            flat_rles,
            tax_ids,
            kmer_len,
            kmer_to_rle_range,
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

        info!("unflattening database...");
        let kmers_and_rles = self
            .kmer_to_rle_range
            .par_iter()
            .map(|(kmer, (range_start, offset))| {
                let start = *range_start as usize;
                let end = start + *offset as usize;
                (*kmer, self.flat_rles[start..end].to_vec())
            })
            .collect::<Vec<(u32, Vec<u16>)>>();

        info!("database unflattened! performing lossy compression...");

        let total_set_bits = RunLengthEncodingIter::from_blocks(&self.flat_rles).count();
        debug!("total set bits before compression {}", total_set_bits);

        let kmers_and_compressed_rles = kmers_and_rles
            .into_par_iter()
            .map(|(kmer, current_rle)| {
                // variable to hold the new lossy compressed blocks as u16s
                let mut compressed_blocks = vec![];

                // peekable iterator over the current runs
                let mut block_iter = current_rle
                    .into_iter()
                    .map(|block| Block::from_u16(block))
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
                // We can "return" it as a tuple to collect after the map statement
                (
                    kmer,
                    compressed_blocks
                        .into_iter()
                        .map(|block| block.to_u16())
                        .collect_vec(),
                )
            })
            .collect::<Vec<(u32, Vec<u16>)>>(); // end map

        info!("re-flattening database...");
        let (kmer_to_rle_range, flat_rles) = create_flat_index(kmers_and_compressed_rles);
        info!("database re-flattened!");
        self.flat_rles = flat_rles;
        self.kmer_to_rle_range = kmer_to_rle_range;

        let total_set_bits = RunLengthEncodingIter::from_blocks(&self.flat_rles).count();
        info!("total set bits after compression {}", total_set_bits);

        // Recompute the p_values after
        self.recompute_p_values();
    }

    fn recompute_p_values(&mut self) -> () {
        let total_canonical_kmers =
            (4_usize.pow(self.kmer_len as u32) - 4_usize.pow(self.kmer_len.div_ceil(2) as u32)) / 2;

        let mut file2kmer_num = vec![0_usize; self.files.len()];

        for (_kmer, (range_start, offset)) in self.kmer_to_rle_range.iter() {
            let start = *range_start as usize;
            let end = start + *offset as usize;
            for kmer in RunLengthEncodingIter::from_blocks(&self.flat_rles[start..end]) {
                file2kmer_num[kmer] += 1;
            }
        }

        let p_values = file2kmer_num
            .into_par_iter()
            .map(|kmer_num| kmer_num as f64 / total_canonical_kmers as f64)
            .collect::<Vec<f64>>();

        self.p_values = p_values;
    }

    pub fn classify(
        &self,
        read: &[u8],
        cutoff_threshold: BigExpFloat,
        n_max: u64,
    ) -> Option<(&str, usize)> {
        let mut collected_hits = vec![0_u64; self.files.len()];

        // Find the hits for all kmers
        let mut n_total = 0_u64;
        for kmer in KmerIter::from(read, self.kmer_len, self.canonical).map(|kmer| kmer as u32) {
            if let Some((range_start, offset)) = self.kmer_to_rle_range.get(&kmer) {
                let start = *range_start as usize;
                let end = start + *offset as usize;
                for sequence in RunLengthEncodingIter::from_blocks(&self.flat_rles[start..end]) {
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
                // This check tries to save runtime in practice
                // Only do probability computation if the p-value is going to be < 0.5
                if n_hits as f64 > (n_total as f64 * p) {
                    let x = if n_total <= n_max {
                        n_hits
                    } else {
                        ((n_hits as f64 / n_total as f64) * n_max as f64).round() as u64
                    };

                    let n = if n_total <= n_max { n_total } else { n_max };

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
            // If (for whatever reason) two probabilities are the same, this will use the first one
            if probability < lowest_prob {
                (lowest_prob_index, lowest_prob) = (index, probability);
            }
        }

        if lowest_prob < cutoff_threshold {
            Some((
                &*self.files[lowest_prob_index],
                self.tax_ids[lowest_prob_index],
            ))
        } else {
            None
        }
    }
}
