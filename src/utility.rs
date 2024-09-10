use bio::io::{fasta, fastq};
use bio::utils::TextSlice;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::fs::File;
use std::fs::{self, DirEntry};
use std::io::BufReader;
use std::path::Path;
use std::{collections::HashSet, path::PathBuf};
use tracing::{error, info, warn};

use crate::kmer_iter::KmerIter;

pub const XOR_NUMBER: usize = 188_888_881;

fn is_fasta_file(entry: &DirEntry) -> bool {
    let entry_file_name = entry.file_name().to_str().unwrap().to_string();
    entry_file_name.ends_with(".fna")
        || entry_file_name.ends_with(".fasta")
        || entry_file_name.ends_with(".fa")
}

pub fn get_fasta_files(reference_loc: &Path) -> Vec<PathBuf> {
    let dir_content = fs::read_dir(reference_loc).expect("could not read reference directory");
    dir_content
        .par_bridge()
        .into_par_iter()
        .filter_map(|dir_entry| match dir_entry {
            Ok(entry) => {
                if entry.path().is_file() && is_fasta_file(&entry) {
                    Some(entry.path())
                } else {
                    warn!(
                        "directory entry {:?} did not end with '.fna', '.fasta', or '.fa', skipping...",
                        entry
                    );
                    None
                }
            }
            Err(e) => {
                error!(
                    "error encountered while reading reference directory {:?}",
                    reference_loc
                );
                error!("{}", e);
                warn!("attempting to continue execution...");
                None
            }
        })
        .collect::<Vec<PathBuf>>()
}

pub fn convert_to_uppercase(sequence: TextSlice) -> String {
    match std::str::from_utf8(sequence) {
        Err(error) => panic!("{}", error),
        Ok(some_str) => some_str.to_uppercase(),
    }
}

pub fn get_fasta_iter_of_file(file_path: &Path) -> fasta::Records<BufReader<File>> {
    match fasta::Reader::from_file(file_path) {
        Ok(reader) => reader.records(),
        Err(error) => panic!("{}", error),
    }
}

pub fn get_fastq_iter_of_file(file_path: &Path) -> fastq::Records<BufReader<File>> {
    match fastq::Reader::from_file(file_path) {
        Ok(reader) => reader.records(),
        Err(error) => panic!("{}", error),
    }
}

pub fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|x| {
            if *x == b'A' || *x == b'a' {
                b'T'
            } else if *x == b'C' || *x == b'c' {
                b'G'
            } else if *x == b'G' || *x == b'g' {
                b'C'
            } else if *x == b'T' || *x == b't' {
                b'A'
            } else {
                *x
            }
        })
        .collect()
}

pub fn two_bit_dna_representation(c: &u8) -> Option<usize> {
    if *c == b'A' || *c == b'a' {
        // binary is 00
        Some(0)
    } else if *c == b'C' || *c == b'c' {
        // binary is 01
        Some(1)
    } else if *c == b'G' || *c == b'g' {
        // binary is 10
        Some(2)
    } else if *c == b'T' || *c == b't' {
        // binary is 11
        Some(3)
    } else {
        return None;
    }
}

pub fn compressed_representation(kmer: &[u8]) -> Option<usize> {
    let mut num = 0_usize;
    for c in kmer.iter() {
        num <<= 2;
        match two_bit_dna_representation(c) {
            None => {
                return None;
            }
            Some(n) => {
                num |= n;
            }
        }
    }
    Some(num)
}

pub fn get_range(kmer_length: usize, log_blocks: u32, block_index: usize) -> (usize, usize) {
    let n_blocks = 2_usize.pow(log_blocks);
    if block_index >= n_blocks {
        panic!(
            "Block index needs to be < {}. Block index {} was chosen.",
            n_blocks, block_index
        );
    }
    let block_size = 4_usize.pow(kmer_length as u32) / n_blocks;
    info!("{} blocks with size {} each", n_blocks, block_size);
    let (lowest_kmer, highest_kmer) = (block_size * block_index, block_size * (block_index + 1));
    info!(
        "accepting kmers in the range [{}, {})",
        lowest_kmer, highest_kmer
    );
    (lowest_kmer, highest_kmer)
}

// Creates a single bitmap containing k-mers from all files, if necessary
pub fn create_bitmap(files: Vec<PathBuf>, kmer_len: usize, canonical: bool) -> RoaringBitmap {
    let mut bitmap = RoaringBitmap::new();
    for file in files {
        let mut record_iter = get_fasta_iter_of_file(&file);
        while let Some(Ok(record)) = record_iter.next() {
            if record.seq().len() < kmer_len {
                continue;
            }
            for kmer in KmerIter::from(record.seq(), kmer_len, canonical) {
                bitmap.insert(kmer as u32);
            }
        }
    }
    bitmap
}

pub fn greedy_ordering(distances: &Vec<Vec<u32>>, start_index: usize) -> Vec<usize> {
    let mut connected_indices = HashSet::from([start_index]);
    let mut ordering = vec![start_index];
    let mut current_index = start_index;

    while ordering.len() < distances.len() {
        let mut next_index = 0_usize;
        let mut next_index_minimum = u32::MAX;

        let mut distance_iter = distances[current_index]
            .iter()
            .chain(
                distances[(current_index + 1)..]
                    .iter()
                    .map(|row| &row[current_index]),
            )
            .enumerate();

        while let Some((index, distance)) = distance_iter.next() {
            if *distance < next_index_minimum && !connected_indices.contains(&index) {
                next_index = index;
                next_index_minimum = *distance;
            }
        }

        ordering.push(next_index);
        connected_indices.insert(next_index);
        current_index = next_index;
    }

    ordering
}

pub fn ordering_statistics(ordering: &Vec<usize>, distances: &Vec<Vec<u32>>) -> (f64, u64) {
    let sum = ordering
        .windows(2)
        .map(|x| {
            if x[0] < x[1] {
                distances[x[1]][x[0]] as u64
            } else {
                distances[x[0]][x[1]] as u64
            }
        })
        .sum();
    (sum as f64 / (ordering.len() - 1) as f64, sum)
}
