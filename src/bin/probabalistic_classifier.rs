use std::collections::{HashMap, HashSet};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::path::Path;
use probabilitic_classifier::utility::{convert_to_uppercase, create_fasta_iterator_from_file};
use clap::{Parser};
use log::info;
use statrs::distribution::{Binomial, DiscreteCDF};
use probabilitic_classifier::taxonomy::get_accession_to_tax_id;

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
    k_mer_size: usize,
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
    let mut database = HashMap::new();

    // Insert points into the database
    let mut ref_fasta_iter = create_fasta_iterator_from_file(reference_file);
    while let Some(Ok(record)) = ref_fasta_iter.next() {
        let accession = String::from(record.id());
        let accession_set = match database.get_mut(&*accession) {
            None => {
                database.insert(accession.clone(), HashSet::new());
                database.get_mut(&*accession).unwrap()
            },
            Some(set) => set
        };
        let uppercase_record_seq = convert_to_uppercase(record.seq());
        for index in 0..(uppercase_record_seq.len() - args.k_mer_size) {
            let kmer = &uppercase_record_seq[index..(index + args.k_mer_size)];
            if kmer.contains("N") {
                continue
            }
            let mut hasher = DefaultHasher::new();
            kmer.hash(&mut hasher);
            accession_set.insert(hasher.finish());
        }
    }

    let mut probabilities = HashMap::new();
    for (accession, kmer_set) in &database {
        let probability = kmer_set.len() as f64 / 4_usize.pow(args.k_mer_size as u32) as f64;
        probabilities.insert(accession.clone(), probability);
    }

    info!("Beginning classificaiton");
    let mut reads_fasta_iter = create_fasta_iterator_from_file(reads_file);
    while let Some(Ok(read)) = reads_fasta_iter.next() {
        let mut read_kmers = HashSet::new();
        let uppercase_record_seq = convert_to_uppercase(read.seq());
        for index in 0..(uppercase_record_seq.len() - args.k_mer_size) {
            let kmer = &uppercase_record_seq[index..(index + args.k_mer_size)];
            if kmer.contains("N") {
                continue
            }
            let mut hasher = DefaultHasher::new();
            kmer.hash(&mut hasher);
            read_kmers.insert( hasher.finish());
        }

        let num_queries = read_kmers.len();
        let mut hit_counts = HashMap::new();
        for kmer in read_kmers {
            for (accession, kmer_set) in &database {
                if kmer_set.contains(&kmer) {
                    match hit_counts.get_mut(&*accession) {
                        None => {
                            hit_counts.insert(accession.clone(), 0);
                            *hit_counts.get_mut(&*accession).unwrap() += 1
                        },
                        Some(count) => *count += 1,
                    };
                }
            }
        }

        let mut best_prob = f64::MAX;
        let mut best_prob_tax_id = 0_u32;
        let mut zero_seen = false;
        for (accession, num_hits) in hit_counts {
            let binomial = Binomial::new(*probabilities.get(&*accession).unwrap(), num_queries as u64).unwrap();
            let prob = binomial.sf(num_hits).abs();
            if prob == 0.0 {
                zero_seen = true;
                println!("{}\t{}", read.id(), accession2taxid.get(&*accession).unwrap());
            }
            // println!("{}", prob);
            // println!("taxid={}, binomial with p={}, n={}, x={}", tax_id, *probabilities.get(tax_id).unwrap(), num_queries, num_hits);
            if prob < best_prob {
                best_prob = prob;
                best_prob_tax_id = accession2taxid.get(&*accession).unwrap().clone();
            }
        }
        if !zero_seen {
            if best_prob < 0.000000000000001 {
                println!("{}\t{}", read.id(), best_prob_tax_id);
            } else {
                println!("{}\t0", read.id())
            }
        }
    } // end read iterator
    info!("Done!")

} // end main