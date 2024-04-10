use clap::Parser;
use log::info;
use musk::intersect::IntersectIterator;
use musk::io::{dump_data_to_file, load_string2taxid};
use musk::kmer_iter::KmerIter;
use musk::utility::get_fasta_iterator_of_file;
use std::collections::HashSet;
use std::path::Path;
use std::sync::mpsc;
use threadpool::ThreadPool;

fn create_bit_vector(files: &str, kmer_length: usize) -> Vec<u32> {
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
    let mut bit_vector = bitset.into_iter().collect::<Vec<u32>>();
    bit_vector.sort();
    bit_vector
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
    /// Location to output the serialzed distances
    output_file: String,

    #[arg()]
    /// the file2taxid file
    file2taxid: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let file2taxid_path = Path::new(&args.file2taxid);
    let output_file_path = Path::new(&args.output_file);

    info!("loading files2taxid at {}", args.file2taxid);
    let file2taxid = load_string2taxid(file2taxid_path);
    info!("creating bitmaps for each group...");
    let mut bit_vectors = vec![];
    let (sender, receiver) = mpsc::channel();
    let pool = ThreadPool::new(args.thread_number);
    for (files, taxid) in file2taxid {
        let sender_clone = sender.clone();
        pool.execute(move || {
            let bitmap = create_bit_vector(&*files, args.kmer_length);
            sender_clone.send((bitmap, files, taxid)).unwrap();
        })
    }
    drop(sender);
    for triple in receiver {
        bit_vectors.push(triple);
    }

    info!("bitmaps computed, computing pairwise distances...");

    let mut all_distances = vec![vec![]; bit_vectors.len()];
    for index_1 in 0..bit_vectors.len() {
        for index_2 in 0..bit_vectors.len() {
            if index_2 <= index_1 {
                continue;
            }
            let distance =
                IntersectIterator::from(&bit_vectors[index_1].0, &bit_vectors[index_2].0).count();
            all_distances[index_1].push(distance);
        }
    }

    let data_dump = (
        bit_vectors
            .into_iter()
            .map(|(_, file, taxid)| (file, taxid))
            .collect::<Vec<(String, u32)>>(),
        all_distances,
    );

    dump_data_to_file(
        bincode::serialize(&data_dump).unwrap(),
        output_file_path,
    )
    .unwrap();
}
