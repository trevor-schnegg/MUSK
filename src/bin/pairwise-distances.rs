use clap::Parser;
use log::{debug, info};
use musk::io::{dump_data_to_file, load_string2taxid};
use musk::kmer_iter::KmerIter;
use musk::sorted_vector_sets::IntersectIterator;
use musk::utility::get_fasta_iterator_of_file;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::{mpsc, Arc};
use std::time::{Duration, Instant};
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
    let (sender, receiver) = mpsc::channel();
    let bit_vectors_arc = Arc::new(bit_vectors);
    for index_1 in 0..bit_vectors_arc.len() {
        let sender_clone = sender.clone();
        let bit_vectors_arc_clone = bit_vectors_arc.clone();
        pool.execute(move || {
            let mut total_duration = Duration::new(0, 0);
            let mut total_intersection_count = 0;
            let mut distances = vec![];
            for index_2 in 0..bit_vectors_arc_clone.len() {
                if index_2 <= index_1 {
                    continue;
                }
                let (bit_vector_1, bit_vector_2) = (
                    &bit_vectors_arc_clone[index_1].0,
                    &bit_vectors_arc_clone[index_2].0,
                );
                let start = Instant::now();
                let intersection_size = IntersectIterator::from(bit_vector_1, bit_vector_2).count();
                total_duration += start.elapsed();
                total_intersection_count += 1;
                // |A| + |B| - 2 * |A and B|
                distances.push(
                    (bit_vector_1.len() + bit_vector_2.len() - (2 * intersection_size)) as u32,
                );
            }
            sender_clone
                .send((
                    index_1,
                    (total_duration, total_intersection_count),
                    distances,
                    bit_vectors_arc_clone[index_1].1.clone(),
                    bit_vectors_arc_clone[index_1].2,
                ))
                .unwrap();
        })
    }
    drop(sender);
    let mut map = HashMap::new();
    let mut total_duration = Duration::new(0, 0);
    let mut total_computations = 0;
    for (index, duration_information, distances, file, taxid) in receiver {
        total_duration += duration_information.0;
        total_computations += duration_information.1;
        all_distances[index] = distances;
        map.insert(index, (file, taxid));
        if map.len() % 1000 == 0 {
            debug!("done with {} sequences", map.len());
        }
    }
    debug!(
        "average time to compute intersection was {:?}",
        total_duration / total_computations
    );

    let data_dump = all_distances
        .into_iter()
        .enumerate()
        .map(|(index, distances)| {
            (
                distances,
                map.get(&index).unwrap().0.clone(),
                map.get(&index).unwrap().1,
            )
        })
        .collect::<Vec<(Vec<u32>, String, u32)>>();

    dump_data_to_file(bincode::serialize(&data_dump).unwrap(), output_file_path).unwrap();
}
