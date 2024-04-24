use clap::Parser;
use itertools::Itertools;
use log::{debug, info};
use musk::io::{dump_data_to_file, load_string2taxid};
use musk::kmer_iter::KmerIter;
use musk::utility::get_fasta_iterator_of_file;
use roaring::RoaringBitmap;
use std::collections::HashSet;
use std::path::Path;
use std::sync::{mpsc, Arc};
use threadpool::ThreadPool;

fn create_bitmap(files: String, kmer_length: usize, taxid: u32) -> (RoaringBitmap, String, u32) {
    let mut kmer_set = HashSet::new();
    for file in files.split(",") {
        let mut record_iter = get_fasta_iterator_of_file(Path::new(&file));
        while let Some(Ok(record)) = record_iter.next() {
            if record.seq().len() < kmer_length {
                continue;
            }
            for kmer in KmerIter::from(record.seq(), kmer_length) {
                kmer_set.insert(kmer as u32);
            }
        }
    }
    let mut kmers = Vec::from_iter(kmer_set.into_iter());
    kmers.sort();
    (RoaringBitmap::from_sorted_iter(kmers).unwrap(), files, taxid)
}

fn distance(size_1: u64, size_2: u64, intersection_size: u64) -> u32 {
    // |A| + |B| - (2 * |A & B|)
    (size_1 + size_2 - (2 * intersection_size)) as u32
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
    info!("creating sorted kmer vectors for each group...");
    let mut sequences = vec![];
    let (sender, receiver) = mpsc::channel();
    let pool = ThreadPool::new(args.thread_number);
    for (files, taxid) in file2taxid {
        let sender_clone = sender.clone();
        pool.execute(move || {
            let sequence = create_bitmap(files, args.kmer_length, taxid);
            sender_clone.send(sequence).unwrap();
        })
    }
    drop(sender);
    for triple in receiver {
        sequences.push(triple);
    }

    info!("roaring bitmaps computed, computing pairwise distances...");

    let mut all_distances = sequences
        .iter()
        .map(|(_bitmap, file, taxid)| (vec![0_u32; sequences.len()], file.clone(), *taxid))
        .collect_vec();
    let (sender, receiver) = mpsc::channel();
    let sequences_arc = Arc::new(sequences);
    for sequence_index_1 in 0..sequences_arc.len() {
        let sender_clone = sender.clone();
        let sequences_arc_clone = sequences_arc.clone();
        pool.execute(move || {
            for sequence_index_2 in 0..sequences_arc_clone.len() {
                if sequence_index_2 <= sequence_index_1 {
                    continue;
                }
                let (sequence_1, sequence_2) = (&sequences_arc_clone[sequence_index_1].0, &sequences_arc_clone[sequence_index_2].0);
                let intersection_size = sequence_1.intersection_len(sequence_2);
                let distance = distance(sequence_1.len(), sequence_2.len(), intersection_size);
                sender_clone.send((sequence_index_1, sequence_index_2, distance)).unwrap();
            }
        });
    }
    drop(sender);
    let mut maximum_distance_computations = 0;
    for (distance_computation, (index_1, index_2, distance)) in receiver.into_iter().enumerate() {
        all_distances[index_1].0[index_2] = distance;
        all_distances[index_2].0[index_1] = distance;
        if distance_computation % 1000000 == 0 {
            debug!("done with {} distance computations", distance_computation);
        }
        maximum_distance_computations += 1
    }
    info!(
        "completed {} distance computations",
        maximum_distance_computations
    );

    let data = bincode::serialize(&all_distances).unwrap();
    dump_data_to_file(data, output_file_path).unwrap();
}
