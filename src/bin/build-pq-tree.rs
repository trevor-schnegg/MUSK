use clap::Parser;
use log::{debug, info};
use musk::io::load_string2taxid;
use musk::kmer_iter::KmerIter;
use musk::utility::get_fasta_iterator_of_file;
use pq_tree::PQTree;
use std::collections::HashSet;
use std::path::Path;
use std::sync::mpsc;
use threadpool::ThreadPool;

fn create_bitmap(files: &str, kmer_length: usize) -> HashSet<u32> {
    let mut bitset = HashSet::new();
    for file in files.split(",") {
        let mut record_iter = get_fasta_iterator_of_file(Path::new(&file));
        while let Some(Ok(record)) = record_iter.next() {
            if record.seq().len() < kmer_length {
                continue;
            }
            for kmer in KmerIter::from(record.seq(), kmer_length) {
                bitset.insert(kmer as u32);
            }
        }
    }
    bitset
}

/// Creates a file to tax id mapping where files with the same tax id are grouped
/// together if their k-mer spectra are similar enough.
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 14)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg(short, long, default_value_t = 12)]
    /// Length of k-mer to use in the database
    thread_number: usize,

    #[arg()]
    /// the file2taxid file
    file2taxid: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let file2taxid_path = Path::new(&args.file2taxid);
    let mut pq_tree = PQTree::from_leaves(&[0, 1, 2, 3, 4, 5, 6, 7]).unwrap();
    let mut matrix_by_row = vec![vec![]; 4_usize.pow(args.kmer_length as u32)];

    info!("loading files2taxid at {}", args.file2taxid);
    let file2taxid = load_string2taxid(file2taxid_path);
    info!("creating bitmaps for each group...");
    let (sender, receiver) = mpsc::channel();
    let pool = ThreadPool::new(args.thread_number);
    for (files, taxid) in file2taxid {
        let sender_clone = sender.clone();
        pool.execute(move || {
            let bitmap = create_bitmap(&*files, args.kmer_length);
            sender_clone.send((bitmap, files, taxid)).unwrap();
        })
    }
    drop(sender);
    for (index, (kmer_set, _file, _taxid)) in receiver.iter().enumerate() {
        for kmer in kmer_set {
            matrix_by_row[kmer as usize].push(index);
        }
    }

    let mut constraints_inserted = 0_usize;
    for (kmer, sequence_vector) in matrix_by_row.iter().enumerate() {
        if !sequence_vector.is_empty() {
            match pq_tree.reduction(sequence_vector) {
                Ok(tree_node) => {
                    constraints_inserted += 1;
                    pq_tree = tree_node;
                },
                Err(x) => {
                    debug!("failed on kmer {}", kmer);
                    debug!("{} actual constraints inserted", constraints_inserted);
                    panic!("{:?}", x);
                }
            }
        }
    }
}
