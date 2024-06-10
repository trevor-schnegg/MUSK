use itertools::Itertools;
use serde::Deserialize;
use std::any::type_name;
use std::fs::File;
use std::io::Read;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::{io, vec};
use tracing::{error, warn};

pub fn split_string_to_taxid(line: String) -> Result<(String, usize), String> {
    let split_line = line.split("\t").collect::<Vec<&str>>();
    let file = split_line[0].to_string();

    // If there is no tab character, index 1 will not exist
    match split_line.get(1) {
        None => Err("line did not have a tab character".to_string()),
        Some(str) => {
            // Try to parse as a usize
            match str.parse::<usize>() {
                Ok(taxid) => Ok((file, taxid)),
                Err(e) => Err(e.to_string()),
            }
        }
    }
}

pub fn load_string2taxid(string2taxid: &Path) -> Vec<(String, usize)> {
    let file = File::open(string2taxid).expect(&*format!(
        "could not read string2taxid tsv at {:?}",
        string2taxid
    ));
    let reader = BufReader::new(file).lines();

    reader
        .enumerate()
        .filter_map(|(line_num, line)| {
            // Try to get the line from the file
            match line {
                Ok(line) => {
                    // Try to parse the line into a (file, taxid) tuple
                    match split_string_to_taxid(line) {
                        Ok(tuple) => Some(tuple),
                        Err(msg) => {
                            error!("{}", msg);
                            warn!(
                                "line {} from {:?} gave the previous error. skipping...",
                                line_num, string2taxid
                            );
                            None
                        }
                    }
                }
                Err(e) => {
                    error!("{:?}", e);
                    warn!(
                        "line {} from {:?} gave the previous error. skipping...",
                        line_num, string2taxid
                    );
                    None
                }
            }
        })
        .collect_vec()
}

pub fn dump_data_to_file(data: Vec<u8>, file: &Path) -> io::Result<()> {
    let mut f = File::create(file).expect(&*format!("could not create file {:?}", file));
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
