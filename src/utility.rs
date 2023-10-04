use bio::io::fasta;
use bio::io::fasta::Records;
use bio::utils::TextSlice;
use statrs::function::beta;
use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

pub enum Sequence {
    Single(String),
    Double(String, String),
}

pub fn get_fasta_files(reference_loc: &Path) -> Vec<String> {
    let dir_content =
        fs::read_dir(reference_loc).expect("could not read provided reference directory");
    let fasta_files = dir_content
        .filter_map(|x| {
            if let Ok(entry) = x {
                if entry.path().is_file() && entry.file_name().to_str().unwrap().ends_with(".fna") {
                    Some(entry.path().to_str().unwrap().to_string())
                } else {
                    None
                }
            } else {
                None
            }
        })
        .collect::<Vec<String>>();
    fasta_files
}

pub fn convert_to_uppercase(sequence: TextSlice) -> String {
    match std::str::from_utf8(sequence) {
        Err(error) => {
            panic!(
                "Unable to convert record sequence to uppercase, \
                the following error was returned:\n\t{}",
                error
            );
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

pub fn reverse_complement(string: &str) -> String {
    let mut return_string = String::new();
    for character in string.chars().rev() {
        if character == 'A' {
            return_string += "T"
        } else if character == 'T' {
            return_string += "A"
        } else if character == 'C' {
            return_string += "G"
        } else if character == 'G' {
            return_string += "C"
        } else {
            return_string += &*String::from(character.clone())
        }
    }
    return_string
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
    let mut num = 0;
    for (kmer_position, n) in kmer.iter().rev().enumerate() {
        if *n == b'A' {
            // binary is 00, don't need to change anything
            continue;
        } else if *n == b'C' {
            // binary is 01
            num |= 1_u32 << (kmer_position << 1)
        } else if *n == b'G' {
            // binary is 10
            num |= 2_u32 << (kmer_position << 1)
        } else if *n == b'T' {
            // binary is 11
            num |= 3_u32 << (kmer_position << 1)
        } else {
            return None;
        }
    }
    Some(num)
}

pub fn get_kmers_as_u32(sequence: Sequence, kmer_len: usize) -> HashSet<u32> {
    let mut kmer_set = HashSet::new();
    match sequence {
        Sequence::Single(sequence) => {
            for kmer in sequence
                .as_bytes()
                .windows(kmer_len)
                .map(|kmer_bytes| vec_dna_bytes_to_u32(kmer_bytes))
            {
                match kmer {
                    None => continue,
                    Some(int) => {
                        kmer_set.insert(int);
                    }
                }
            }
        }
        Sequence::Double(forward, reverse) => {
            for (kmer_1_int, kmer_2_int) in forward
                .as_bytes()
                .windows(kmer_len)
                .zip(reverse.as_bytes().windows(kmer_len))
                .map(|(kmer_1_bytes, kmer_2_bytes)| {
                    (
                        vec_dna_bytes_to_u32(kmer_1_bytes),
                        vec_dna_bytes_to_u32(kmer_2_bytes),
                    )
                })
            {
                match (kmer_1_int, kmer_2_int) {
                    (Some(int_1), Some(int_2)) => {
                        kmer_set.insert(int_1);
                        kmer_set.insert(int_2);
                    }
                    (Some(int), None) => {
                        kmer_set.insert(int);
                    }
                    (None, Some(int)) => {
                        kmer_set.insert(int);
                    }
                    (None, None) => continue,
                }
            }
        }
    }
    kmer_set
}
