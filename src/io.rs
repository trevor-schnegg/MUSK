use log::{error, warn};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::{io, vec};

pub fn load_string2taxid(string2taxid: &Path) -> HashMap<String, u32> {
    let open_file =
        File::open(string2taxid).expect(&*format!("could not read tsv at {:?}", string2taxid));
    let reader = BufReader::new(open_file).lines();
    let mut string2taxid = HashMap::new();
    for (line_number, line) in reader.enumerate() {
        match line {
            Ok(line) => {
                let split_line = line.split("\t").collect::<Vec<&str>>();
                string2taxid.insert(
                    split_line[0].to_string(),
                    split_line[1].parse::<u32>().unwrap(),
                );
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

pub fn load_taxid2files(taxid2files: &Path) -> HashMap<u32, Vec<String>> {
    let open_file =
        File::open(taxid2files).expect(&*format!("could not read tsv at {:?}", taxid2files));
    let reader = BufReader::new(open_file).lines();
    let mut taxid2files = HashMap::new();
    for (line_number, line) in reader.enumerate() {
        match line {
            Ok(line) => {
                let split_line = line.split("\t").collect::<Vec<&str>>();
                let taxid = split_line[0].parse::<u32>().unwrap();
                match taxid2files.get_mut(&taxid) {
                    None => {
                        taxid2files.insert(taxid, vec![split_line[1].to_string()]);
                    }
                    Some(vector) => {
                        vector.push(split_line[1].to_string());
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
