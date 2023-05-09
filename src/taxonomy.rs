use crate::utility::create_fasta_iterator_from_file;
use bio::io::fasta::Records;
use log::info;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

pub fn get_needed_accessions(mut fasta_reader: Records<BufReader<File>>) -> HashSet<String> {
    let mut accession_set = HashSet::new();
    while let Some(Ok(record)) = fasta_reader.next() {
        accession_set.insert(String::from(record.id()));
    }
    accession_set
}

pub fn create_accession_to_tax_id_map<P: AsRef<Path>>(
    accession_set: HashSet<String>,
    accession_to_tax_id_dir: P,
) -> HashMap<String, u32> {
    let directory = accession_to_tax_id_dir.as_ref();
    let mut nucl_gb_reader =
        BufReader::new(File::open(directory.join("nucl_gb.accession2taxid")).unwrap()).lines();
    let mut nucl_wgs_reader =
        BufReader::new(File::open(directory.join("nucl_wgs.accession2taxid")).unwrap()).lines();
    let mut nucl_extra_reader =
        BufReader::new(File::open(directory.join("nucl_extra.accession2taxid")).unwrap()).lines();

    let mut accession_to_tax_id: HashMap<String, u32> = HashMap::new();

    while let Some(Ok(line)) = nucl_gb_reader.next() {
        let split_line: Vec<&str> = line.split("\t").collect();
        let accession = *split_line.get(1).unwrap();
        match accession_set.contains(accession) {
            true => accession_to_tax_id.insert(
                String::from(accession),
                split_line.get(2).unwrap().parse().unwrap(),
            ),
            false => continue,
        };
    }
    while let Some(Ok(line)) = nucl_wgs_reader.next() {
        let split_line: Vec<&str> = line.split("\t").collect();
        let accession = *split_line.get(1).unwrap();
        match accession_set.contains(accession) {
            true => {
                accession_to_tax_id.insert(
                    String::from(accession),
                    split_line.get(2).unwrap().parse().unwrap(),
                );
                ()
            }
            false => continue,
        }
    }
    while let Some(Ok(line)) = nucl_extra_reader.next() {
        let split_line: Vec<&str> = line.split("\t").collect();
        let accession = *split_line.get(1).unwrap();
        match accession_set.contains(accession) {
            true => {
                accession_to_tax_id.insert(
                    String::from(accession),
                    split_line.get(2).unwrap().parse().unwrap(),
                );
                ()
            }
            false => continue,
        }
    }
    accession_to_tax_id
}

pub fn dump_accession_to_tax_id<P: AsRef<Path>>(
    taxonomy_dir: P,
    accession_to_tax_id: &HashMap<String, u32>,
) {
    let taxonomy_dir = taxonomy_dir.as_ref();
    let mut f = File::create(taxonomy_dir.join("needed_accession2taxid"))
        .expect("Could not create accession2taxid file");
    let data_to_write = bincode::serialize(accession_to_tax_id)
        .expect("Could not serialize accession_to_tax_id map");
    f.write_all(&*data_to_write).unwrap();
}

pub fn load_accession_to_tax_id<P: AsRef<Path>>(taxonomy_dir: P) -> HashMap<String, u32> {
    let mut f = File::open(taxonomy_dir.as_ref().join("needed_accession2taxid")).unwrap();
    let mut buf: Vec<u8> = vec![];
    f.read_to_end(&mut buf).unwrap();
    bincode::deserialize(&*buf).unwrap()
}

pub fn get_accession_to_tax_id(taxonomy_dir: &Path, fasta_file: &Path) -> HashMap<String, u32> {
    match File::open(taxonomy_dir.join("needed_accession2taxid")) {
        Ok(_) => {
            info!("'needed_accession2taxid' file found, checking that all needed accessions are present...");
            let mut accession_to_tax_id = load_accession_to_tax_id(taxonomy_dir);
            let mut fasta_iter = create_fasta_iterator_from_file(fasta_file);
            while let Some(Ok(record)) = fasta_iter.next() {
                match accession_to_tax_id.get(record.id()) {
                    None => {
                        info!("not all accessions found, creating new 'needed_accession2taxid' file...");
                        accession_to_tax_id = create_accession_to_tax_id_map(
                            get_needed_accessions(create_fasta_iterator_from_file(fasta_file)),
                            taxonomy_dir,
                        );
                        dump_accession_to_tax_id(taxonomy_dir, &accession_to_tax_id);
                        info!("new 'needed_accession2taxid' file created!");
                        break;
                    }
                    Some(_) => continue,
                }
            }
            info!("all needed accessions are present!");
            accession_to_tax_id
        }
        Err(_) => {
            info!("no 'needed_accession2taxid' file found, creating...");
            let accession_to_tax_id = create_accession_to_tax_id_map(
                get_needed_accessions(create_fasta_iterator_from_file(fasta_file)),
                taxonomy_dir,
            );
            dump_accession_to_tax_id(taxonomy_dir, &accession_to_tax_id);
            info!("'needed_accession2taxid' file created!");
            accession_to_tax_id
        }
    }
}

