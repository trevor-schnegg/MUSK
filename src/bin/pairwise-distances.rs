use clap::Parser;
use log::{debug, info};
use musk::io::{dump_data_to_file, load_string2taxid};
use musk::kmer_iter::KmerIter;
use musk::utility::get_fasta_iterator_of_file;
use std::cmp::{min, Reverse};
use std::collections::HashSet;
use std::path::Path;
use std::sync::{mpsc, Arc};
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

fn geometric_search(sorted_vector: &[u32], value: u32) -> Result<usize, usize> {
    let mut end = 0;
    while sorted_vector[end] < value && end < sorted_vector.len() - 1 {
        end = min(if end == 0 { 1 } else { end << 1 }, sorted_vector.len() - 1);
    }
    if sorted_vector[end] == value {
        Ok(end)
    } else {
        match sorted_vector[end / 2..=end].binary_search(&value) {
            Ok(index) => Ok((end / 2) + index),
            Err(index) => Err((end / 2) + index),
        }
    }
}

fn intersect(vector_1: &Vec<u32>, vector_2: &Vec<u32>) -> u32 {
    let mut intersect_size = 0;
    let (mut index_1, mut index_2) = (0_usize, 0_usize);
    while index_1 < vector_1.len() && index_2 < vector_2.len() {
        let (current_1, current_2) = (vector_1[index_1], vector_2[index_2]);
        if current_1 < current_2 {
            match geometric_search(&vector_1[index_1..], current_2) {
                Ok(index) => {
                    intersect_size += 1;
                    index_1 += index + 1;
                }
                Err(index) => index_1 += index,
            }
        } else if current_2 < current_1 {
            match geometric_search(&vector_2[index_2..], current_1) {
                Ok(index) => {
                    intersect_size += 1;
                    index_2 += index + 1;
                }
                Err(index) => index_2 += index,
            }
        } else {
            panic!("impossible case reached");
        }
    }
    intersect_size
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

    let mut all_distances = vec![];
    let bitsets_arc = Arc::new(bit_vectors);
    let (sender, receiver) = mpsc::channel();
    for index_1 in 0..bitsets_arc.len() {
        let bitmaps_arc_clone = bitsets_arc.clone();
        let sender_clone = sender.clone();
        pool.execute(move || {
            let mut distances = vec![];
            for index_2 in 0..bitmaps_arc_clone.len() {
                if index_2 <= index_1 {
                    continue;
                }
                let (bitset_1, bitset_2) =
                    (&bitmaps_arc_clone[index_1].0, &bitmaps_arc_clone[index_2].0);
                let intersection_size = intersect(bitset_1, bitset_2);
                // |A| + |B| - 2 * |A and B|
                distances
                    .push(bitset_1.len() as u32 + bitset_2.len() as u32 - (2 * intersection_size));
            }
            sender_clone
                .send((
                    distances,
                    bitmaps_arc_clone[index_1].1.clone(),
                    bitmaps_arc_clone[index_1].2,
                ))
                .unwrap();
        })
    }
    drop(sender);
    drop(bitsets_arc);
    for distances in receiver {
        all_distances.push(distances);

        if all_distances.len() % 1000 == 0 {
            debug!("completed {} sequences", all_distances.len());
        }
    }

    info!("all distances computed! sorting...");
    all_distances.sort_by_key(|x| Reverse(x.0.len()));
    info!("sorting completed! outputting to file...");

    dump_data_to_file(
        bincode::serialize(&all_distances).unwrap(),
        output_file_path,
    )
    .unwrap();
}
