use std::collections::{HashMap, HashSet};
use std::path::Path;
use bio::alphabets;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::fmindex::BackwardSearchResult::Complete;
use bio::data_structures::suffix_array::suffix_array;
use probabilitic_classifier::utility::{convert_to_uppercase, create_fasta_iterator_from_file, reverse_complement};
use clap::{Parser};
use log::info;
use num_traits::Float;
use f128::f128;
use statrs::distribution::DiscreteCDF;
use probabilitic_classifier::binomial::{Binomial, ulps_eq};
use probabilitic_classifier::taxonomy::get_accession_to_tax_id;
use probabilitic_classifier::database::{Database};

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

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let taxonomy_dir = Path::new(&args.taxonomy_dir);
    let reference_file = Path::new(&args.reference_file);
    let reads_file = Path::new(&args.reads_file);

    // Get accession2taxid
    info!("Loading accession2taxid from taxonomy directory: {}", args.taxonomy_dir);
    let accession2taxid = get_accession_to_tax_id(taxonomy_dir, reference_file);
    info!("accession2taxid loaded!");
    let mut database_stats = Database::new(args.kmer_len);

    info!("Creating database string and probabilities");
    // Create database string and get probabilities
    let mut accession2probabilities = HashMap::new();
    let mut database_string = String::new();
    let mut record_iter = create_fasta_iterator_from_file(reference_file);
    while let Some(Ok(record)) = record_iter.next() {
        let database_string_start = database_string.len();
        let uppercase_record_seq = convert_to_uppercase(record.seq());
        let reverse_complement_seq = reverse_complement(&*uppercase_record_seq);
        let prob = database_stats.calculate_probability((&*uppercase_record_seq, &*reverse_complement_seq));
        accession2probabilities.insert(record.id().to_string(), prob);
        database_string += &*(uppercase_record_seq + "$");
        database_string += &*(reverse_complement_seq + "$");
        database_stats.push_interval((database_string_start, database_string.len() -1, record.id().to_string()));
    }
    info!("Database string and probabilities created!");

    // Set probabilities
    database_stats.set_probabilities(accession2probabilities);

    info!("Creating FM index");
    // Create fm index
    let alphabet = alphabets::dna::n_alphabet();
    let sa = suffix_array(database_string.as_bytes());
    let bwt = bwt(database_string.as_bytes(), &sa);
    let less = less(&bwt, &alphabet);
    let occ = Occ::new(&bwt, 3, &alphabet);
    let fmindex = FMIndex::new(&bwt, &less, &occ);
    info!("FM index created!");

    let needed_probability = f128::from(1.0e-100);
    info!("Beginning classification");
    let mut read_iter = create_fasta_iterator_from_file(reads_file);
    while let Some(Ok(read)) = read_iter.next() {
        // Find all distinct kmers in the read
        let mut read_kmers = HashSet::new();
        let uppercase_record_seq = convert_to_uppercase(read.seq());
        for index in 0..(uppercase_record_seq.len() - args.kmer_len) {
            let kmer = &uppercase_record_seq[index..(index + args.kmer_len)];
            if !kmer.contains("N") {
                read_kmers.insert(kmer);
            }
        }
        // Query all distinct kmers in the FM Index
        let num_queries = read_kmers.len() as u64;
        let mut hit_counts = HashMap::new();
        for kmer in read_kmers {
            match fmindex.backward_search(kmer.as_bytes().iter()) {
                Complete(interval) => {
                    let mut accessions = HashSet::new();
                    for accession in interval.occ(&sa).iter().map(|x| database_stats.get_accession_of_index(*x)) {
                        accessions.insert(accession);
                    }
                    for accession in accessions {
                        match hit_counts.get_mut(&accession) {
                            None => {hit_counts.insert(accession, 1);},
                            Some(count) => {*count += 1},
                        }
                    }
                }
                _ => continue,
            }
        }

        let mut best_prob = f128::MAX;
        let mut best_prob_index = 0_u32;
        let mut zero_seen = false;
        for (accession, num_hits) in hit_counts {
            let binomial = Binomial::new(f128::from(database_stats.get_probability(&*accession)), num_queries).unwrap();
            let prob = binomial.sf(num_hits).abs();
            if ulps_eq(f128::ZERO, prob) {
                zero_seen = true;
                println!("{}\t{}", read.id(), accession2taxid.get(&*accession).unwrap());
            }
            // println!("{}", prob);
            // println!("taxid={}, binomial with p={}, n={}, x={}", tax_id, *probabilities.get(tax_id).unwrap(), num_queries, num_hits);
            if prob < best_prob {
                best_prob = prob;
                best_prob_index = accession2taxid.get(&*accession).unwrap().clone();
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