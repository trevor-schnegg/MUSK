use bio::io::fasta;
use bio::io::fasta::Records;
use bio::utils::TextSlice;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

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

pub fn create_fasta_iterator_from_file(file_path: &Path) -> Records<BufReader<File>> {
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
