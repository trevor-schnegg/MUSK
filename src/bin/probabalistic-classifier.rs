use clap::Parser;
use log::info;
use probabilitic_classifier::binomial::Binomial;
use probabilitic_classifier::taxonomy::get_accession_to_tax_id;
use probabilitic_classifier::utility::{
    convert_to_uppercase, create_fasta_iterator_from_file, reverse_complement,
};
use rug::Float;
use std::collections::hash_map::DefaultHasher;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::path::Path;

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

fn insert_kmers(set_ref: &mut HashSet<u64>, seq: String, kmer_len: usize) {
    for index in 0..(seq.len() - kmer_len) {
        let kmer = &seq[index..(index + kmer_len)];
        if kmer.contains("N") {
            continue;
        }
        let mut hasher = DefaultHasher::new();
        kmer.hash(&mut hasher);
        set_ref.insert(hasher.finish());
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
    let mut database = HashMap::new();

    // Insert points into the database
    for record in create_fasta_iterator_from_file(reference_file) {
        let record = record.unwrap();
        let accession_set = match database.get_mut(record.id()) {
            None => {
                database.insert(String::from(record.id()), HashSet::new());
                database.get_mut(record.id()).unwrap()
            }
            Some(set) => set,
        };
        let uppercase_record_seq = convert_to_uppercase(record.seq());
        let reverse_complement_seq = reverse_complement(&*uppercase_record_seq);
        insert_kmers(accession_set, uppercase_record_seq, args.kmer_len);
        insert_kmers(accession_set, reverse_complement_seq, args.kmer_len);
    }

    let mut probabilities = HashMap::new();
    for (accession, kmer_set) in &database {
        let probability = kmer_set.len() as f64 / 4_usize.pow(args.kmer_len as u32) as f64;
        probabilities.insert(accession.clone(), probability);
    }

    let needed_probability = Float::with_val(256, 1.0e-100);
    info!("Beginning classification");
    let mut read_iter = create_fasta_iterator_from_file(reads_file);
    while let Some(Ok(read)) = read_iter.next() {
        let mut read_kmers = HashSet::new();
        let uppercase_record_seq = convert_to_uppercase(read.seq());
        insert_kmers(&mut read_kmers, uppercase_record_seq, args.kmer_len);

        let num_queries = read_kmers.len();
        let mut hit_counts = HashMap::new();
        for kmer in read_kmers {
            for (accession, kmer_set) in &database {
                if kmer_set.contains(&kmer) {
                    match hit_counts.get_mut(&*accession) {
                        None => {
                            hit_counts.insert(accession.clone(), 0);
                            *hit_counts.get_mut(&*accession).unwrap() += 1
                        }
                        Some(count) => *count += 1,
                    };
                }
            }
        }

        let mut best_prob = Float::with_val(256, 1.0);
        let mut best_prob_tax_id = 0_u32;
        for (accession, num_hits) in hit_counts {
            let binomial = Binomial::new(
                Float::with_val(256, *probabilities.get(&*accession).unwrap()),
                num_queries as u64,
            )
            .unwrap();
            let prob = binomial.sf(num_hits).abs();
            if prob < best_prob {
                best_prob = prob;
                best_prob_tax_id = accession2taxid.get(&*accession).unwrap().clone();
            }
        }
        if best_prob < needed_probability {
            println!("{}\t{}", read.id(), best_prob_tax_id);
        } else {
            println!("{}\t0", read.id())
        }
    } // end read iterator
    info!("Done!")
} // end main
