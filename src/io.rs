use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

pub fn load_accession2taxid(accession2taxid_path: &Path) -> HashMap<String, u32> {
    let mut accession2taxid = HashMap::new();
    let file = File::open(accession2taxid_path).expect("could not read accession2taxid");
    for line in BufReader::new(file).lines() {
        if let Ok(line) = line {
            let split_line = line.split("\t").collect::<Vec<&str>>();
            let accession = split_line.get(0).unwrap().to_string();
            let taxid = split_line.get(1).unwrap().parse::<u32>().unwrap();
            accession2taxid.insert(accession, taxid);
        }
    }
    accession2taxid
}

pub fn dump_data_to_file(data: Vec<u8>, file: &Path) -> io::Result<()> {
    let mut f = File::create(file).expect(&*format!("Could not create file {:?}", file));
    f.write_all(&*data)
}
