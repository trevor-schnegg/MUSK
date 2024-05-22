use log::{error, warn};
use serde::Deserialize;
use std::any::type_name;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::{io, vec};

pub fn split_string_to_taxid(line: String) -> (String, u32) {
    let split_line = line.split("\t").collect::<Vec<&str>>();
    let file_path = split_line[0].to_string();
    let taxid = split_line[1].parse::<u32>().unwrap();
    (file_path, taxid)
}

pub fn load_string2taxid(string2taxid: &Path) -> Vec<(String, u32)> {
    let open_file =
        File::open(string2taxid).expect(&*format!("could not read tsv at {:?}", string2taxid));
    let reader = BufReader::new(open_file).lines();
    let mut string2taxid = Vec::new();
    for (line_number, line) in reader.enumerate() {
        match line {
            Ok(line) => {
                string2taxid.push(split_string_to_taxid(line));
            }
            Err(error) => {
                warn!(
                    "line {} from {:?} gave the following error:",
                    line_number, string2taxid
                );
                error!("{:?}", error);
                warn!("skipping line {}...", line_number);
            }
        }
    }
    string2taxid
}

pub fn load_taxid2files(file2taxid: &Path) -> HashMap<u32, Vec<String>> {
    let open_file =
        File::open(file2taxid).expect(&*format!("could not read tsv at {:?}", file2taxid));
    let reader = BufReader::new(open_file).lines();
    let mut taxid2files = HashMap::new();
    for (line_number, line) in reader.enumerate() {
        match line {
            Ok(line) => {
                let (file_path, taxid) = split_string_to_taxid(line);
                match taxid2files.get_mut(&taxid) {
                    None => {
                        taxid2files.insert(taxid, vec![file_path]);
                    }
                    Some(vector) => {
                        vector.push(file_path);
                    }
                }
            }
            Err(error) => {
                warn!(
                    "line {} from {:?} gave the following error:",
                    line_number, taxid2files
                );
                error!("{:?}", error);
                warn!("skipping line {}...", line_number);
            }
        }
    }
    taxid2files
}

pub fn dump_data_to_file(data: Vec<u8>, file: &Path) -> io::Result<()> {
    let mut f = File::create(file).expect(&*format!("Could not create file {:?}", file));
    f.write_all(&*data)
}

pub fn load_data_from_file<T: for<'a> Deserialize<'a>>(path: &Path) -> T {
    let mut f = File::open(path).expect(&*format!("could not open file at {:?}", path));
    let mut buf: Vec<u8> = vec![];
    f.read_to_end(&mut buf).unwrap();
    bincode::deserialize(&*buf).expect(&*format!(
        "failed to deserialize data at {:?} into {}",
        path,
        type_name::<T>()
    ))
}
