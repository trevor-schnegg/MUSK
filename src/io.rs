use log::{error, warn};
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Lines, Write};
use std::path::Path;

struct TsvIter {
    tsv_iter: Lines<BufReader<File>>,
}

impl TsvIter {
    fn new(path: &Path) -> Self {
        let file = File::open(path).expect(&*format!("could not read tsv at {:?}", path));
        TsvIter {
            tsv_iter: BufReader::new(file).lines(),
        }
    }
}

impl Iterator for TsvIter {
    type Item = Vec<String>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.tsv_iter.next() {
            None => None,
            Some(Ok(line)) => Some(line.split("\t").map(|x| x.to_string()).collect()),
            Some(Err(x)) => {
                error!("{}", x);
                warn!("skipping the line that caused this error...");
                self.next()
            }
        }
    }
}

pub fn load_string2taxid(string2taxid: &Path) -> HashMap<String, u32> {
    let mut accession2taxid = HashMap::new();
    let file = File::open(string2taxid).expect("could not read accession2taxid");
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

pub fn load_accession2taxid(accession2taxid_path: &Path) -> HashMap<String, u32> {
    let mut tsv_iter = TsvIter::new(accession2taxid_path);
    let mut accession2taxid = HashMap::new();
    while let Some(vec) = tsv_iter.next() {
        let accession = vec.get(0).unwrap().clone();
        let taxid = vec.get(1).unwrap().parse::<u32>().unwrap();
        accession2taxid.insert(accession, taxid);
    }
    accession2taxid
}

pub fn load_taxid2files(file2taxid_path: &Path) -> HashMap<u32, Vec<String>> {
    let mut tsv_iter = TsvIter::new(file2taxid_path);
    let mut taxid2file = HashMap::new();
    while let Some(vec) = tsv_iter.next() {
        let file = vec.get(0).unwrap().clone();
        let taxid = vec.get(1).unwrap().parse::<u32>().unwrap();
        match taxid2file.get_mut(&taxid) {
            None => {
                taxid2file.insert(taxid, vec![file]);
            }
            Some(file_vec) => {
                file_vec.push(file);
            }
        }
    }
    taxid2file
}

pub fn dump_data_to_file(data: Vec<u8>, file: &Path) -> io::Result<()> {
    let mut f = File::create(file).expect(&*format!("Could not create file {:?}", file));
    f.write_all(&*data)
}
