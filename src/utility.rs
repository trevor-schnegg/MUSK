use bio::io::fasta;
use bio::io::fasta::Records;
use bio::utils::TextSlice;
use log::{error, warn};
use statrs::function::beta;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

pub fn get_fasta_files(reference_loc: &Path) -> Vec<String> {
    let dir_content =
        fs::read_dir(reference_loc).expect("could not read provided reference directory");
    let fasta_files = dir_content
        .filter_map(|x| {
            match x {
                Ok(entry) => {
                    if entry.path().is_file() && entry.file_name().to_str().unwrap().ends_with(".fna") {
                        Some(entry.path().to_str().unwrap().to_string())
                    } else {
                        warn!("Found directory entry '{:?}' that did not end with '.fna'. Make sure this is intentional.", entry);
                        None
                    }
                },
                Err(e) => {
                    warn!("Error found while reading contents of directory {:?}", reference_loc);
                    warn!("{}", e);
                    None
                }
            }
        })
        .collect::<Vec<String>>();
    fasta_files
}

pub fn convert_to_uppercase(sequence: TextSlice) -> String {
    match std::str::from_utf8(sequence) {
        Err(error) => {
            error!(
                "Unable to convert record sequence to uppercase, \
                the following error was returned:\n\t{}",
                error
            );
            panic!()
        }
        Ok(some_str) => some_str.to_uppercase(),
    }
}

pub fn get_fasta_iterator_of_file(file_path: &Path) -> Records<BufReader<File>> {
    let returned_reader = match fasta::Reader::from_file(file_path) {
        Ok(reader) => reader,
        Err(error) => {
            panic!("{}", error)
        }
    };
    returned_reader.records()
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

pub fn sf(n: u64, p: f64, x: u64) -> f64 {
    if x >= n {
        1.0
    } else {
        let k = x;
        beta::beta_reg(k as f64 + 1.0, (n - k) as f64, p)
    }
}

pub fn vec_dna_bytes_to_u32(kmer: &[u8]) -> Option<u32> {
    let mut num = 0_u32;
    for (base_position, n) in kmer.iter().rev().enumerate() {
        if *n == b'A' || *n == b'a' {
            // binary is 00, don't need to change anything
            continue;
        } else if *n == b'C' || *n == b'c' {
            // binary is 01
            num |= 1_u32 << (base_position << 1)
        } else if *n == b'G' || *n == b'g' {
            // binary is 10
            num |= 2_u32 << (base_position << 1)
        } else if *n == b'T' || *n == b't' {
            // binary is 11
            num |= 3_u32 << (base_position << 1)
        } else {
            return None;
        }
    }
    Some(num)
}
