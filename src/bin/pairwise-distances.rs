use clap::Parser;
use log::{debug, info};
use musk::io::{dump_data_to_file, load_string2taxid};
use musk::kmer_iter::KmerIter;
use musk::sorted_vector_utilities::{DifferenceIterator, IntersectIterator, UnionIterator};
use musk::utility::get_fasta_iterator_of_file;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::mpsc::Sender;
use std::sync::{mpsc, Arc};
use threadpool::ThreadPool;

enum Sequence {
    One(Vec<u32>, String, u32),
    Many(Vec<u32>, Vec<(Vec<u32>, String)>, u32),
}

fn distance(length_1: usize, length_2: usize, intersection_size: usize) -> u32 {
    (length_1 + length_2 - (2 * intersection_size)) as u32
}

fn create_bit_vectors(files: &str, kmer_length: usize, taxid: u32) -> Sequence {
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
        Sequence::One(sequence.0, sequence.1, taxid)
    } else {
        let union = UnionIterator::from(
            sorted_kmer_vectors
                .iter()
                .map(|vector| &vector.0)
                .collect::<Vec<&Vec<u32>>>(),
        )
        .map(|kmer| *kmer)
        .collect::<Vec<u32>>();
        let difference_vectors = sorted_kmer_vectors
            .into_iter()
            .map(|(sorted_kmer_vector, file)| {
                (
                    DifferenceIterator::from(&union, vec![&sorted_kmer_vector])
                        .map(|kmer| *kmer)
                        .collect::<Vec<u32>>(),
                    file,
                )
            })
            .collect::<Vec<(Vec<u32>, String)>>();
        Sequence::Many(union, difference_vectors, taxid)
    }
}

fn create_file_to_index_map(
    sequences: &Vec<Sequence>,
) -> (Vec<(String, u32)>, HashMap<String, usize>) {
    let files = sequences
        .into_iter()
        .map(|sequence| match sequence {
            Sequence::One(_set, file, taxid) => vec![(file.clone(), *taxid)],
            Sequence::Many(_union, files, taxid) => files
                .into_iter()
                .map(|(_difference, file)| (file.clone(), *taxid))
                .collect::<Vec<(String, u32)>>(),
        })
        .flatten()
        .collect::<Vec<(String, u32)>>();
    let map = HashMap::from_iter(
        files
            .iter()
            .enumerate()
            .map(|(index, (file, _taxid))| (file.clone(), index)),
    );
    (files, map)
}

fn self_matrix(
    many_sequences: (&Vec<u32>, &Vec<(Vec<u32>, String)>),
    sender: &Sender<(usize, usize, u32)>,
    file_to_index: &Arc<HashMap<String, usize>>,
) -> () {
    let (union, difference_vectors) = many_sequences;
    for index_1 in 0..difference_vectors.len() {
        for index_2 in 0..difference_vectors.len() {
            if index_2 <= index_1 {
                continue;
            }
            let (difference_1, difference_2) =
                (&difference_vectors[index_1], &difference_vectors[index_2]);
            let intersection_size =
                DifferenceIterator::from(union, vec![&difference_1.0, &difference_2.0]).count();
            let distance = distance(
                union.len() - difference_1.0.len(),
                union.len() - difference_2.0.len(),
                intersection_size,
            );
            let (sequence_index_1, sequence_index_2) = (
                *file_to_index.get(&difference_1.1).unwrap(),
                *file_to_index.get(&difference_2.1).unwrap(),
            );
            sender
                .send((sequence_index_1, sequence_index_2, distance))
                .unwrap();
        }
    }
}

fn one_to_one(
    sequence_1: (&Vec<u32>, &String),
    sequence_2: (&Vec<u32>, &String),
    sender: &Sender<(usize, usize, u32)>,
    file_to_index: &Arc<HashMap<String, usize>>,
) -> () {
    let intersection_size = IntersectIterator::from(&sequence_1.0, &sequence_2.0).count();
    let distance = distance(sequence_1.0.len(), sequence_2.0.len(), intersection_size);
    let (sequence_index_1, sequence_index_2) = (
        *file_to_index.get(sequence_1.1).unwrap(),
        *file_to_index.get(sequence_2.1).unwrap(),
    );
    sender
        .send((sequence_index_1, sequence_index_2, distance))
        .unwrap();
}

fn many_to_one(
    many_sequences: (&Vec<u32>, &Vec<(Vec<u32>, String)>),
    one: (&Vec<u32>, &String),
    sender: &Sender<(usize, usize, u32)>,
    file_to_index: &Arc<HashMap<String, usize>>,
) -> () {
    let (union, differences) = many_sequences;
    let intersection = IntersectIterator::from(union, one.0)
        .map(|kmer| *kmer)
        .collect::<Vec<u32>>();
    for difference in differences {
        let intersection_size =
            DifferenceIterator::from(&intersection, vec![&difference.0]).count();
        let distance = distance(
            union.len() - difference.0.len(),
            one.0.len(),
            intersection_size,
        );
        let (sequence_index_1, sequence_index_2) = (
            *file_to_index.get(one.1).unwrap(),
            *file_to_index.get(&difference.1).unwrap(),
        );
        sender
            .send((sequence_index_1, sequence_index_2, distance))
            .unwrap();
    }
}

fn many_to_many(
    many_sequences_1: (&Vec<u32>, &Vec<(Vec<u32>, String)>),
    many_sequences_2: (&Vec<u32>, &Vec<(Vec<u32>, String)>),
    sender: &Sender<(usize, usize, u32)>,
    file_to_index: &Arc<HashMap<String, usize>>,
) -> () {
    let (union_1, differences_1) = many_sequences_1;
    let (union_2, differences_2) = many_sequences_2;
    let intersection = IntersectIterator::from(union_1, union_2)
        .map(|kmer| *kmer)
        .collect::<Vec<u32>>();
    for difference_1 in differences_1 {
        for difference_2 in differences_2 {
            let intersection_size =
                DifferenceIterator::from(&intersection, vec![&difference_1.0, &difference_2.0])
                    .count();
            let distance = distance(
                union_1.len() - difference_1.0.len(),
                union_2.len() - difference_2.0.len(),
                intersection_size,
            );
            let (sequence_index_1, sequence_index_2) = (
                *file_to_index.get(&difference_1.1).unwrap(),
                *file_to_index.get(&difference_2.1).unwrap(),
            );
            sender
                .send((sequence_index_1, sequence_index_2, distance))
                .unwrap();
        }
    }
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
            let sequence = create_bit_vectors(&*files, args.kmer_length, taxid);
            sender_clone.send(sequence).unwrap();
        })
    }
    drop(sender);
    for triple in receiver {
        sequences.push(triple);
    }

    info!("sorted kmer vectors computed, computing pairwise distances...");

    let (index_to_file_and_taxid, file_to_index) = create_file_to_index_map(&sequences);

    let mut all_distances = index_to_file_and_taxid
        .into_iter()
        .map(|(file, taxid)| (vec![0_u32; file_to_index.len()], file, taxid))
        .collect::<Vec<(Vec<u32>, String, u32)>>();
    let (sender, receiver) = mpsc::channel();
    let sequences_arc = Arc::new(sequences);
    let pool_arc = Arc::new(pool);
    let file_to_index_arc = Arc::new(file_to_index);
    for index_1 in 0..sequences_arc.len() {
        let sender_clone = sender.clone();
        let sequences_arc_clone = sequences_arc.clone();
        let file_to_index_arc_clone = file_to_index_arc.clone();
        pool_arc.execute(move || {
            for index_2 in 0..sequences_arc_clone.len() {
                if index_2 < index_1 {
                    continue;
                } else if index_2 == index_1 {
                    if let Sequence::Many(union, differences, _) = &sequences_arc_clone[index_1] {
                        self_matrix(
                            (union, differences),
                            &sender_clone,
                            &file_to_index_arc_clone,
                        )
                    }
                } else {
                    match (&sequences_arc_clone[index_1], &sequences_arc_clone[index_2]) {
                        (
                            Sequence::Many(union_1, differences_1, _taxid_1),
                            Sequence::Many(union_2, differences_2, _taxid_2),
                        ) => {
                            many_to_many(
                                (union_1, differences_1),
                                (union_2, differences_2),
                                &sender_clone,
                                &file_to_index_arc_clone,
                            );
                        }
                        (
                            Sequence::Many(union, differences, _taxid_1),
                            Sequence::One(set, file, _taxid_2),
                        ) => {
                            many_to_one(
                                (union, differences),
                                (set, file),
                                &sender_clone,
                                &file_to_index_arc_clone,
                            );
                        }
                        (
                            Sequence::One(set, file, _taxid_1),
                            Sequence::Many(union, differences, _taxid_2),
                        ) => {
                            many_to_one(
                                (union, differences),
                                (set, file),
                                &sender_clone,
                                &file_to_index_arc_clone,
                            );
                        }
                        (
                            Sequence::One(set_1, file_1, _taxid_1),
                            Sequence::One(set_2, file_2, _taxid_2),
                        ) => {
                            one_to_one(
                                (set_1, file_1),
                                (set_2, file_2),
                                &sender_clone,
                                &file_to_index_arc_clone,
                            );
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
        if distance_computation % 100000 == 0 {
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
