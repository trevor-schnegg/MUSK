use bio::io::fasta;
use bio::io::fasta::Records;
use bio::utils::TextSlice;
use log::{debug, error, info, warn};
use roaring::RoaringBitmap;
use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use crate::kmer_iter::KmerIter;

pub const XOR_NUMBER: usize = 188_888_881;

pub fn get_fasta_files(reference_loc: &Path) -> Vec<String> {
    let dir_content =
        fs::read_dir(reference_loc).expect("could not read provided reference directory");
    dir_content
        .filter_map(|dir_entry| match dir_entry {
            Ok(entry) => {
                if entry.path().is_file() && entry.file_name().to_str().unwrap().ends_with(".fna") {
                    Some(entry.path().to_str().unwrap().to_string())
                } else {
                    warn!(
                        "Found directory entry '{:?}' that did not end with '.fna'. Skipping...",
                        entry
                    );
                    None
                }
            }
            Err(e) => {
                warn!(
                    "Error found while reading entry of reference directory {:?}",
                    reference_loc
                );
                error!("{}", e);
                warn!("Skipping the entry because of the above error");
                None
            }
        })
        .collect::<Vec<String>>()
}

pub fn convert_to_uppercase(sequence: TextSlice) -> String {
    match std::str::from_utf8(sequence) {
        Err(error) => panic!("{}", error),
        Ok(some_str) => some_str.to_uppercase(),
    }
}

pub fn get_fasta_iterator_of_file(file_path: &Path) -> Records<BufReader<File>> {
    match fasta::Reader::from_file(file_path) {
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

pub fn create_bitmap(
    files: &str,
    kmer_length: usize,
    lowest_kmer: usize,
    highest_kmer: usize,
    xor: bool,
    canonical: bool,
) -> RoaringBitmap {
    let mut bitmap = RoaringBitmap::new();
    for file in files.split(",") {
        let mut record_iter = get_fasta_iterator_of_file(Path::new(&file));
        while let Some(Ok(record)) = record_iter.next() {
            if record.seq().len() < kmer_length {
                continue;
            }
            for kmer in KmerIter::from(record.seq(), kmer_length, canonical) {
                let kmer = if xor { kmer ^ XOR_NUMBER } else { kmer };
                if lowest_kmer <= kmer && kmer < highest_kmer {
                    bitmap.insert(kmer as u32);
                }
            }
        }
    }
    bitmap
}

pub fn find_ordering(distances: &Vec<(Vec<u32>, String, u32)>, start_index: usize) -> Vec<usize> {
    let mut connected_indices = HashSet::from([start_index]);
    let mut ordering = vec![start_index];
    let mut current_index = start_index;
    while ordering.len() < distances.len() {
        let mut next_index = 0_usize;
        let mut next_index_minimum = u32::MAX;
        for (index, distance) in distances[current_index].0.iter().enumerate() {
            if *distance < next_index_minimum && !connected_indices.contains(&index) {
                next_index = index;
                next_index_minimum = *distance;
            }
        }
        ordering.push(next_index);
        connected_indices.insert(next_index);
        current_index = next_index;
        if ordering.len() % 2500 == 0 {
            debug!("found ordering for {} bitmaps", ordering.len());
        }
    }
    ordering
}

pub fn average_hamming_distance(
    ordering: &Vec<usize>,
    distances: &Vec<(Vec<u32>, String, u32)>,
) -> (f64, u64) {
    let sum = ordering
        .windows(2)
        .map(|x| distances[x[0]].0[x[1]] as u64)
        .sum();
    (sum as f64 / (ordering.len() - 1) as f64, sum)
}
