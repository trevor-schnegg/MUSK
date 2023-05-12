use std::collections::{HashMap, HashSet};
use std::path::Path;
use probabilitic_classifier::utility::{convert_to_uppercase, create_fasta_iterator_from_file, reverse_complement};
use clap::{Parser};
use log::info;
use num_traits::Float;
use f128::f128;
use statrs::distribution::DiscreteCDF;
use probabilitic_classifier::binomial::{Binomial, ulps_eq};
use probabilitic_classifier::taxonomy::get_accession_to_tax_id;
use probabilitic_classifier::database::Database;

/// Converts a fasta file to a database
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg()]
    /// Name of the file to create a reference from
    reference_file: String,

    #[arg()]
    /// Name of the file to create a reference from
    reads_file: String,

    #[arg()]
    /// Directory containing names.dmp, nodes.dmp, nucl_extra.accession2taxid,
    /// nucl_gb.accession2taxid, and nucl_wgs.accession2taxid
    taxonomy_dir: String,

    #[arg()]
    /// k-mers to test
    kmer_len: usize,
}

fn insert_kmers(set_ref: &mut HashSet<Vec<u8>>, seq: String, kmer_len: usize) {
    for index in 0..(seq.len() - kmer_len) {
        let kmer = &seq[index..(index + kmer_len)];
        if kmer.contains("N") {
            continue
        }
        set_ref.insert(Vec::from(kmer.as_bytes()));
    }
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let taxonomy_dir = Path::new(&args.taxonomy_dir);
    let reference_file = Path::new(&args.reference_file);
    let reads_file = Path::new(&args.reads_file);

    // Get accession2taxid and initialize the "database"
    let accession2taxid = get_accession_to_tax_id(taxonomy_dir, reference_file);
    let mut database = Database::new(args.kmer_len);

    // Insert points into the database
    let mut record_iter = create_fasta_iterator_from_file(reads_file);
    while let Some(Ok(record)) = record_iter.next() {
        let uppercase_record_seq = convert_to_uppercase(record.seq());
        let reverse_complement_seq = reverse_complement(&*uppercase_record_seq);
        database.insert_records(record.id().to_string(), vec![uppercase_record_seq, reverse_complement_seq])
    }

    database.init_probabilities();

    let needed_probability = f128::from(1.0e-100);
    info!("Beginning classification");
    let mut read_iter = create_fasta_iterator_from_file(reads_file);
    while let Some(Ok(read)) = read_iter.next() {
        let mut read_kmers = HashSet::new();
        let uppercase_record_seq = convert_to_uppercase(read.seq());
        insert_kmers(&mut read_kmers, uppercase_record_seq, args.kmer_len);

        let num_queries = read_kmers.len() as u64;
        let mut hit_counts = HashMap::new();
        for kmer in read_kmers {
            match database.query(&*kmer) {
                None => {continue}
                Some(vec) => {
                    for index in vec {
                        match hit_counts.get_mut(index) {
                            None => {hit_counts.insert(*index, 1_u64);}
                            Some(count) => {*count += 1}
                        }
                    }
                }
            }
        }

        let mut best_prob = f128::MAX;
        let mut best_prob_index = 0_u32;
        let mut zero_seen = false;
        for (index, num_hits) in hit_counts {
            let binomial = Binomial::new(f128::from(database.get_probability_of_index(index as usize)), num_queries).unwrap();
            let prob = binomial.sf(num_hits).abs();
            if ulps_eq(f128::ZERO, prob) {
                zero_seen = true;
                println!("{}\t{}", read.id(), accession2taxid.get(database.get_accession_of_index(index as usize)).unwrap());
            }
            // println!("{}", prob);
            // println!("taxid={}, binomial with p={}, n={}, x={}", tax_id, *probabilities.get(tax_id).unwrap(), num_queries, num_hits);
            if prob < best_prob {
                best_prob = prob;
                best_prob_index = accession2taxid.get(database.get_accession_of_index(index as usize)).unwrap().clone();
            }
        }
        if !zero_seen {
            if best_prob < needed_probability {
                println!("{}\t{}", read.id(), best_prob_index);
            } else {
                println!("{}\t0", read.id())
            }
        }
    } // end read iterator
    info!("Done!")

} // end main