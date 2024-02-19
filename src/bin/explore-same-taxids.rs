use std::collections::HashSet;
use clap::Parser;
use musk::io::load_taxid2files;
use musk::utility::get_fasta_iterator_of_file;
use std::path::Path;
use log::{debug, info};
use musk::kmer_iter::KmerIter;

/// Explores similarities between files with the same species tax id
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 15)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg()]
    /// the file2taxid file
    file2taxid: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let file2taxid_path = Path::new(&args.file2taxid);

    info!("loading file2taxid at {}", args.file2taxid);
    let file2taxid = load_taxid2files(file2taxid_path);
    info!("file2taxid loaded! finding the average Hamming distance between the same tax ids");
    println!("taxid\tavg_set_size\tavg_hamming_distance");
    for (taxid, files) in file2taxid {
        if files.len() == 1 {
            continue;
        }
        let mut kmer_sets = vec![];
        for file in files {
            let mut record_iter = get_fasta_iterator_of_file(Path::new(&file));
            let mut kmer_set = HashSet::new();
            while let Some(Ok(record)) = record_iter.next() {
                if record.seq().len() < args.kmer_length {
                    continue;
                }
                let kmers = KmerIter::from(record.seq(), args.kmer_length);
                for kmer in kmers {
                    kmer_set.insert(kmer);
                }
            }
            kmer_sets.push(kmer_set);
        }
        let avg_set_size = {
            let count = kmer_sets.len();
            let sum: f64 = kmer_sets.iter().map(|x| x.len() as f64).sum();
            sum / count as f64
        };
        let mut hamming_distances = vec![];
        for (i1, kmer_set_1) in kmer_sets.iter().enumerate() {
            for (i2, kmer_set_2) in kmer_sets.iter().enumerate() {
                if i2 <= i1 {
                    continue;
                }
                let intersect_size = kmer_set_1.intersection(&kmer_set_2).map(|x| *x).collect::<Vec<usize>>().len();
                debug!("{}\t{}\t{}", kmer_set_1.len(), kmer_set_2.len(), intersect_size);
                hamming_distances.push((kmer_set_1.len() - intersect_size) + (kmer_set_2.len() - intersect_size))
            }
        }
        let avg_hamming_dist = {
            let count = hamming_distances.len();
            let sum: usize = hamming_distances.iter().sum();
            sum as f64 / count as f64
        };
        println!("{}\t{}\t{}", taxid, avg_set_size, avg_hamming_dist);
    }
}
