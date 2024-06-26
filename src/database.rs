use indicatif::ProgressIterator;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use serde::{Deserialize, Serialize};
use statrs::distribution::{Binomial, DiscreteCDF};
use tracing::{debug, info};

use crate::rle::{NaiveRunLengthEncoding, RunLengthEncoding};

#[derive(Serialize, Deserialize)]
pub struct Database {
    canonical: bool,
    cutoff_threshold: f64,
    file2taxid: Vec<(String, usize)>,
    kmer_len: usize,
    kmer_runs: Vec<RunLengthEncoding>,
    num_queries: u64,
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
            cutoff_threshold,
            file2taxid,
            kmer_len,
            kmer_runs,
            num_queries,
            p_values,
            significant_hits,
        }
    }

    pub fn lossy_compression(&mut self, _num_ones: usize) -> () {
        // Modify the kmer_runs variable in any desired way

        // Recompute the p_values and significant hits after
        self.recompute_statistics();
    }

    fn recompute_statistics(&mut self) -> () {
        let total_num_kmers = 4_usize.pow(self.kmer_len as u32) as f64;

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
                Binomial::new(*p, self.num_queries)
                    .unwrap()
                    .inverse_cdf(self.cutoff_threshold)
            })
            .collect::<Vec<u64>>();

        self.p_values = p_values;
        self.significant_hits = significant_hits;
    }
}
