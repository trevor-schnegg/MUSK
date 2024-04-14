use clap::Parser;
use log::{debug, info};
use musk::explore::connected_components;
use musk::io::load_taxid2files;
use musk::kmer_iter::KmerIter;
use musk::utility::get_fasta_iterator_of_file;
use std::collections::HashSet;
use std::path::Path;
use std::sync::mpsc;
use threadpool::ThreadPool;

fn create_bitmaps(files: &Vec<String>, kmer_length: usize, thread_number: usize) -> Vec<Vec<u32>> {
    let mut sorted_vector_sets = vec![];
    let (sender, receiver) = mpsc::channel();
    let pool = ThreadPool::new(thread_number);
    for (index, file) in files.iter().enumerate() {
        let sender_clone = sender.clone();
        let file = file.clone();
        pool.execute(move || {
            let mut kmer_set = HashSet::new();
            let mut record_iter = get_fasta_iterator_of_file(Path::new(&file));
            while let Some(Ok(record)) = record_iter.next() {
                if record.seq().len() < kmer_length {
                    continue;
                }
                for kmer in KmerIter::from(record.seq(), kmer_length) {
                    kmer_set.insert(kmer as u32);
                }
            }
            let mut sorted_vector = kmer_set.into_iter().collect::<Vec<u32>>();
            sorted_vector.sort();
            sender_clone.send((index, sorted_vector)).unwrap();
        });
    }
    drop(sender);
    for tuple in receiver {
        sorted_vector_sets.push(tuple);
    }
    sorted_vector_sets.sort_by_key(|(i, _)| *i);
    sorted_vector_sets
        .into_iter()
        .map(|(_, sorted_vector_set)| sorted_vector_set)
        .collect::<Vec<Vec<u32>>>()
}

/// Creates a matrix of (hamming) distances between bitmaps
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

    info!("loading file2taxid at {}", args.file2taxid);
    let taxid2files = load_taxid2files(file2taxid_path);
    info!("file2taxid loaded! exploring files with the same tax id");
    for (taxid, files) in taxid2files {
        if files.len() == 1 {
            println!("{}\t{}", files[0], taxid);
            continue;
        }
        debug!(
            "creating sets for taxid '{}' with {} files...",
            taxid,
            files.len()
        );
        let sorted_vectors = create_bitmaps(&files, args.kmer_length, args.thread_number);
        debug!("hashsets created! performing comparisons...");
        let connected_components = connected_components(sorted_vectors, 0.9, args.thread_number);
        for component in connected_components {
            let mut files_string = String::new();
            for file_index in component {
                if files_string.is_empty() {
                    files_string += &*files[file_index];
                } else {
                    files_string += &*(String::from(",") + &*files[file_index])
                }
            }
            println!("{}\t{}", files_string, taxid);
        }
    }
}
