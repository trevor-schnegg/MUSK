use bio::io::fasta;
use bio::io::fasta::Records;
use bio::utils::TextSlice;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

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
