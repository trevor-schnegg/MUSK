use clap::Parser;
use itertools::Itertools;
use log::{debug, info};
use musk::io::{dump_data_to_file, load_string2taxid};
use musk::kmer_iter::KmerIter;
use musk::utility::get_fasta_iterator_of_file;
use roaring::{MultiOps, RoaringBitmap};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::{mpsc, Arc};
use threadpool::ThreadPool;

enum Bitmap {
    One(RoaringBitmap, String, u32),
    Many(RoaringBitmap, Vec<(RoaringBitmap, String)>, u32),
}

fn create_bitmaps(files: &str, kmer_length: usize, taxid: u32) -> Bitmap {
    let mut sorted_kmer_vectors = vec![];
    for file in files.split(",") {
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
        let mut sorted_kmer_vector = kmer_set.into_iter().collect::<Vec<u32>>();
        sorted_kmer_vector.sort();
        sorted_kmer_vectors.push((sorted_kmer_vector, file.to_string()));
    }
    if sorted_kmer_vectors.len() == 1 {
        let sequence = sorted_kmer_vectors.into_iter().nth(0).unwrap();
        Bitmap::One(
            RoaringBitmap::from_sorted_iter(sequence.0).unwrap(),
            sequence.1,
            taxid,
        )
    } else {
        let bitmaps = sorted_kmer_vectors
            .into_iter()
            .map(|(vector, files)| (RoaringBitmap::from_sorted_iter(vector).unwrap(), files))
            .collect_vec();
        let union = bitmaps.iter().map(|x| &x.0).union();
        let difference_vectors = bitmaps
            .into_iter()
            .map(|(bitmap, file)| (&union - bitmap, file))
            .collect_vec();
        Bitmap::Many(union, difference_vectors, taxid)
    }
}

fn create_file_to_index_map(
    sequences: &Vec<Bitmap>,
) -> (Vec<(String, u32)>, HashMap<String, usize>) {
    let files = sequences
        .into_iter()
        .map(|sequence| match sequence {
            Bitmap::One(_set, file, taxid) => vec![(file.clone(), *taxid)],
            Bitmap::Many(_union, files, taxid) => files
                .into_iter()
                .map(|(_difference, file)| (file.clone(), *taxid))
                .collect_vec(),
        })
        .flatten()
        .collect_vec();
    let map = HashMap::from_iter(
        files
            .iter()
            .enumerate()
            .map(|(index, (file, _taxid))| (file.clone(), index)),
    );
    (files, map)
}

fn distance(size_1: u64, size_2: u64, intersection_size: u64) -> u32 {
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
            let sequence = create_bitmaps(&*files, args.kmer_length, taxid);
            sender_clone.send(sequence).unwrap();
        })
    }
    drop(sender);
    for triple in receiver {
        sequences.push(triple);
    }

    let (index_to_file_and_taxid, file_to_index) = create_file_to_index_map(&sequences);

    info!("sorted kmer vectors computed, computing pairwise distances...");

    let mut all_distances = index_to_file_and_taxid
        .into_iter()
        .map(|(file, taxid)| (vec![0_u32; file_to_index.len()], file, taxid))
        .collect_vec();
    let (sender, receiver) = mpsc::channel();
    let sequences_arc = Arc::new(sequences);
    let file_to_index_arc = Arc::new(file_to_index);
    for sequence_index_1 in 0..sequences_arc.len() {
        let sender_clone = sender.clone();
        let sequences_arc_clone = sequences_arc.clone();
        let file_to_index_arc_clone = file_to_index_arc.clone();
        pool.execute(move || {
            for sequence_index_2 in 0..sequences_arc_clone.len() {
                if sequence_index_2 < sequence_index_1 {
                    continue;
                } else if sequence_index_2 == sequence_index_1 {
                    if let Bitmap::Many(union, differences, _) =
                        &sequences_arc_clone[sequence_index_1]
                    {
                        for difference_index_1 in 0..differences.len() {
                            for difference_index_2 in 0..differences.len() {
                                if difference_index_2 <= difference_index_1 {
                                    continue;
                                }
                                let (difference_1, difference_2) = (
                                    &differences[difference_index_1],
                                    &differences[difference_index_2],
                                );
                                let intersection_size = union.difference_len(
                                    &vec![&difference_1.0, &difference_2.0].union(),
                                );
                                let size_1 = union.len() - difference_1.0.len();
                                let size_2 = union.len() - difference_2.0.len();
                                let distance = distance(size_1, size_2, intersection_size);
                                let (send_index_1, send_index_2) = (
                                    *file_to_index_arc_clone.get(&difference_1.1).unwrap(),
                                    *file_to_index_arc_clone.get(&difference_2.1).unwrap(),
                                );
                                sender_clone
                                    .send((send_index_1, send_index_2, distance))
                                    .unwrap();
                            }
                        }
                    }
                } else {
                    match (
                        &sequences_arc_clone[sequence_index_1],
                        &sequences_arc_clone[sequence_index_2],
                    ) {
                        (
                            Bitmap::Many(union_1, differences_1, _taxid_1),
                            Bitmap::Many(union_2, differences_2, _taxid_2),
                        ) => {
                            let union_intersection = vec![union_1, union_2].intersection();
                            for difference_1 in differences_1 {
                                for difference_2 in differences_2 {
                                    let intersection_size = union_intersection.difference_len(
                                        &vec![&difference_1.0, &difference_2.0].union(),
                                    );
                                    let size_1 = union_1.len() - difference_1.0.len();
                                    let size_2 = union_2.len() - difference_2.0.len();
                                    let distance = distance(size_1, size_2, intersection_size);
                                    let (send_index_1, send_index_2) = (
                                        *file_to_index_arc_clone.get(&difference_1.1).unwrap(),
                                        *file_to_index_arc_clone.get(&difference_2.1).unwrap(),
                                    );
                                    sender_clone
                                        .send((send_index_1, send_index_2, distance))
                                        .unwrap();
                                }
                            }
                        }
                        (
                            Bitmap::Many(union, differences, _taxid_1),
                            Bitmap::One(set, file, _taxid_2),
                        ) => {
                            let intersection = vec![union, set].intersection();
                            let size_1 = set.len();
                            for difference in differences {
                                let intersection_size = intersection.difference_len(&difference.0);
                                let size_2 = union.len() - difference.0.len();
                                let distance = distance(size_1, size_2, intersection_size);
                                let (send_index_1, send_index_2) = (
                                    *file_to_index_arc_clone.get(&difference.1).unwrap(),
                                    *file_to_index_arc_clone.get(file).unwrap(),
                                );
                                sender_clone
                                    .send((send_index_1, send_index_2, distance))
                                    .unwrap();
                            }
                        }
                        (
                            Bitmap::One(set, file, _taxid_1),
                            Bitmap::Many(union, differences, _taxid_2),
                        ) => {
                            let intersection = vec![union, set].intersection();
                            let size_1 = set.len();
                            for difference in differences {
                                let intersection_size = intersection.difference_len(&difference.0);
                                let size_2 = union.len() - difference.0.len();
                                let distance = distance(size_1, size_2, intersection_size);
                                let (send_index_1, send_index_2) = (
                                    *file_to_index_arc_clone.get(&difference.1).unwrap(),
                                    *file_to_index_arc_clone.get(file).unwrap(),
                                );
                                sender_clone
                                    .send((send_index_1, send_index_2, distance))
                                    .unwrap();
                            }
                        }
                        (
                            Bitmap::One(set_1, file_1, _taxid_1),
                            Bitmap::One(set_2, file_2, _taxid_2),
                        ) => {
                            let intersection_size = set_1.intersection_len(set_2);
                            let distance = distance(set_1.len(), set_2.len(), intersection_size);
                            let (send_index_1, send_index_2) = (
                                *file_to_index_arc_clone.get(file_1).unwrap(),
                                *file_to_index_arc_clone.get(file_2).unwrap(),
                            );
                            sender_clone
                                .send((send_index_1, send_index_2, distance))
                                .unwrap();
                        }
                    }
                }
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
