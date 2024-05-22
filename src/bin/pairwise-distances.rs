use clap::Parser;
use itertools::Itertools;
use log::{debug, info};
use musk::io::{dump_data_to_file, load_string2taxid};
use musk::kmer_iter::KmerIter;
use musk::utility::get_fasta_iterator_of_file;
use roaring::RoaringBitmap;
use std::path::Path;
use std::sync::{mpsc, Arc};
use threadpool::ThreadPool;
use rayon::prelude::*;

const XOR_NUMBER: usize = 188_888_881;

fn create_bitmap(files: String, kmer_length: usize, taxid: u32, low: usize, high: usize) -> (RoaringBitmap, String, u32) {
    let mut bitmap = RoaringBitmap::new();
    for file in files.split(",") {
        let mut record_iter = get_fasta_iterator_of_file(Path::new(&file));
        while let Some(Ok(record)) = record_iter.next() {
            if record.seq().len() < kmer_length {
                continue;
            }
            for kmer in KmerIter::from(record.seq(), kmer_length) {
                let kmer = kmer ^ XOR_NUMBER;
                if low <= kmer && kmer < high {
                    bitmap.insert(kmer as u32);
                }
            }
        }
    }
    (
        bitmap,
        files,
        taxid,
    )
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

    #[arg(short, long, default_value_t = 0)]
    /// 2^{log_blocks} partitions
    log_blocks: u32,

    #[arg(short, long, default_value_t = 0)]
    /// The index of the block to use
    block_i: usize,

    #[arg(short, long)]
    /// The directory prefix of the fasta files
    old_directory_prefix: Option<String>,

    #[arg(short, long)]
    /// The directory prefix of the fasta files
    new_directory_prefix: Option<String>,

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

    let (lowest_kmer, highest_kmer) = {
        let n_blocks = 2_usize.pow(args.log_blocks);
        if args.block_i >= n_blocks {
            panic!("Block index needs to be < {}. Block index {} was chosen.", n_blocks, args.block_i);
        }
        let block_size = 4_usize.pow(args.kmer_length as u32) / n_blocks;
        debug!("{} blocks with size {}", n_blocks, block_size);
        (block_size * args.block_i, block_size * (args.block_i + 1))
    };
    info!("accepting kmers in the range [{}, {})", lowest_kmer, highest_kmer);

    info!("loading files2taxid at {}", args.file2taxid);
    let mut file2taxid = load_string2taxid(file2taxid_path);
    if let (Some(old_prefix), Some(new_prefix))= (args.old_directory_prefix, args.new_directory_prefix) {
        file2taxid = file2taxid.into_iter().map(|(files, taxid)| (files.replace(&*old_prefix, &*new_prefix), taxid)).collect_vec();
    }
    info!("{} groups total", file2taxid.len());
    info!("creating roaring bitmaps for each group...");
    let sequences = file2taxid.into_par_iter().map(|(files, taxid)| {
        create_bitmap(files, args.kmer_length, taxid, lowest_kmer, highest_kmer)

    }).collect::<Vec<(RoaringBitmap, String, u32)>>();
    info!("roaring bitmaps computed, allocating matrix...");

    let mut all_distances = sequences
        .iter()
        .map(|(_bitmap, file, taxid)| (vec![0_u32; sequences.len()], file.clone(), *taxid))
        .collect_vec();
    info!("matrix allocated, computing pairwise distances...");
    let pool = ThreadPool::new(args.thread_number);
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
                let (sequence_1, sequence_2) = (
                    &sequences_arc_clone[sequence_index_1].0,
                    &sequences_arc_clone[sequence_index_2].0,
                );
                let intersection_size = sequence_1.intersection_len(sequence_2);
                // |A| + |B| - (2 * |A & B|)
                let distance = (sequence_1.len() + sequence_2.len() - (2 * intersection_size)) as u32;
                sender_clone
                    .send((sequence_index_1, sequence_index_2, distance))
                    .unwrap();
            }
        });
    }
    drop(sender);
    let mut total_distance_computations = 0;
    for (distance_computation, (index_1, index_2, distance)) in receiver.into_iter().enumerate() {
        all_distances[index_1].0[index_2] = distance;
        all_distances[index_2].0[index_1] = distance;
        if distance_computation % 1000000 == 0 {
            debug!("done with {} distance computations", distance_computation);
        }
        total_distance_computations += 1
    }
    info!(
        "completed all {} distance computations",
        total_distance_computations
    );

    let data = bincode::serialize(&all_distances).unwrap();
    dump_data_to_file(data, output_file_path).unwrap();
}
