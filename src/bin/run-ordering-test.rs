use clap::Parser;
use indicatif::{ParallelProgressIterator, ProgressIterator};
use itertools::Itertools;
use musk::{
    io::{load_data_from_file, load_string2taxid},
    kmer_iter::KmerIter,
    rle::{BuildRunLengthEncoding, RunLengthEncoding},
    utility::{get_fasta_iterator_of_file, get_range, greedy_ordering, XOR_NUMBER},
};
use rand::{
    distributions::{Distribution, Uniform},
    thread_rng,
};
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::{
    collections::{HashMap, HashSet},
    path::Path,
};
use tracing::info;

fn create_bitmaps(
    files: &str,
    kmer_length: usize,
    blocks: &Vec<(usize, usize, usize)>,
) -> Vec<RoaringBitmap> {
    let mut bitmaps = vec![RoaringBitmap::new(); blocks.len()];
    for file in files.split(",") {
        let mut record_iter = get_fasta_iterator_of_file(Path::new(&file));
        while let Some(Ok(record)) = record_iter.next() {
            if record.seq().len() < kmer_length {
                continue;
            }
            for kmer in KmerIter::from(record.seq(), kmer_length, false) {
                let kmer = kmer ^ XOR_NUMBER;
                for (index, block) in blocks.iter().enumerate() {
                    if block.1 <= kmer && kmer < block.2 {
                        bitmaps[index].insert(kmer as u32);
                        break;
                    }
                }
            }
        }
    }
    bitmaps
}

/// Creates a sample of k-mers from the matrix
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 14)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg(short, long, default_value_t = 0)]
    /// 2^{log_blocks} partitions
    log_blocks: u32,

    #[arg()]
    /// The old directory prefix of the fasta files
    old_directory_prefix: String,

    #[arg()]
    /// The old directory prefix of the fasta files
    new_directory_prefix: String,

    #[arg()]
    /// The ordering of the sequences for the full matrix
    matrix_ordering_file: String,

    #[arg()]
    /// the file2taxid file
    file2taxid: String,
}

fn main() {
    tracing_subscriber::fmt::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let file2taxid_path = Path::new(&args.file2taxid);
    let matrix_ordering_file_path = Path::new(&args.matrix_ordering_file);

    let num_blocks_to_test = 64;
    let max_block_number = 2_usize.pow(args.log_blocks);

    info!(
        "randomly generating {} blocks to test...",
        num_blocks_to_test
    );

    let mut test_block_indices = HashSet::new();
    let mut rng = thread_rng();
    let distribution = Uniform::new(0_usize, max_block_number);
    while test_block_indices.len() < num_blocks_to_test {
        let sample = distribution.sample(&mut rng);
        test_block_indices.insert(sample);
    }

    let blocks = test_block_indices
        .into_iter()
        .map(|block_index| {
            let (low, high) = get_range(args.kmer_length, args.log_blocks, block_index);
            (block_index, low, high)
        })
        .collect_vec();

    info!(
        "blocks generated! loading files2taxid at {} and full matrix ordering at {}...",
        args.file2taxid, args.matrix_ordering_file
    );

    let file2taxid = load_string2taxid(file2taxid_path);
    let matrix_ordering = load_data_from_file::<Vec<(String, u32)>>(matrix_ordering_file_path)
        .into_iter()
        .map(|(files, taxid)| {
            (
                files.replace(&*args.old_directory_prefix, &*args.new_directory_prefix),
                taxid,
            )
        })
        .collect_vec();

    info!(
        "files2taxid and full matrix ordering loaded! creating roaring bitmaps for each group..."
    );

    let all_file_bitmaps = file2taxid
        .into_par_iter()
        .progress()
        .map(|(files, taxid)| {
            let bitmaps = create_bitmaps(&*files, args.kmer_length, &blocks);
            bitmaps
                .into_iter()
                .map(|bitmap| (bitmap, files.clone(), taxid))
                .collect_vec()
        })
        .collect::<Vec<Vec<(RoaringBitmap, String, u32)>>>();
    let mut block_to_bitmaps = HashMap::new();
    for block in blocks.iter() {
        block_to_bitmaps.insert(*block, vec![]);
    }
    for file_bitmaps in all_file_bitmaps {
        for (block, triple) in blocks.iter().zip(file_bitmaps.into_iter()) {
            block_to_bitmaps.get_mut(block).unwrap().push(triple);
        }
    }

    info!("done creating and formatting bitmaps! performing test for each block...");

    for (block, bitmaps) in block_to_bitmaps {
        info!("testing block {}, computing distances...", block.0);

        let mut distances = bitmaps
            .par_iter()
            .progress()
            .enumerate()
            .map(|(index_1, (bitmap_1, files_1, taxid_1))| {
                let inner_distances = bitmaps
                    .par_iter()
                    .enumerate()
                    .filter_map(|(index_2, (bitmap_2, _files_2, _taxid_2))| {
                        if index_2 <= index_1 {
                            None
                        } else {
                            let intersection_size = bitmap_1.intersection_len(bitmap_2);
                            // |A| + |B| - (2 * |A & B|)
                            let distance =
                                (bitmap_1.len() + bitmap_2.len() - (2 * intersection_size)) as u32;
                            Some(distance)
                        }
                    })
                    .collect::<Vec<u32>>();
                (inner_distances, files_1.clone(), *taxid_1)
            })
            .collect::<Vec<(Vec<u32>, String, u32)>>();

        info!(
            "done computing distances for block {}! filling out matrix...",
            block.0
        );

        let mut all_distances: Vec<(Vec<u32>, String, u32)> = vec![];
        for index in 0..distances.len() {
            let mut full_distance_vector = vec![];
            for prior_index in 0..index {
                full_distance_vector.push(all_distances[prior_index].0[index])
            }
            full_distance_vector.push(0);
            full_distance_vector.append(&mut distances[index].0);
            all_distances.push((
                full_distance_vector,
                distances[index].1.clone(),
                distances[index].2.clone(),
            ));
        }
        drop(distances);

        info!(
            "done filling out matrix for block {}, finding ordering...",
            block.0
        );

        let block_ordering = greedy_ordering(&all_distances, 0)
            .into_iter()
            .map(|x| (all_distances[x].1.clone(), all_distances[x].2))
            .collect::<Vec<(String, u32)>>();
        drop(all_distances);

        info!("ordering found for block {}! creating databases to compare to full matrix ordering for block {}...", block.0, block.0);

        let bitmaps: HashMap<String, RoaringBitmap> =
            HashMap::from_iter(bitmaps.into_iter().map(|triple| (triple.1, triple.0)));

        let mut block_ordering_database =
            vec![BuildRunLengthEncoding::new(); 4_usize.pow(args.kmer_length as u32)];
        for (index, (files, _taxid)) in block_ordering.into_iter().progress().enumerate() {
            let bitmap = bitmaps.get(&files).unwrap();
            for kmer in bitmap {
                block_ordering_database[kmer as usize].push(index);
            }
        }
        let naive_block_ordering_runs = block_ordering_database
            .iter()
            .map(|build_rle| build_rle.get_vector().len())
            .sum::<usize>();
        let block_ordering_compressed_database = block_ordering_database
            .into_par_iter()
            .map(|build_rle| build_rle.to_rle())
            .collect::<Vec<RunLengthEncoding>>();
        let block_ordering_compressed_runs = block_ordering_compressed_database
            .iter()
            .map(|rle| rle.get_vector().len())
            .sum::<usize>();
        drop(block_ordering_compressed_database);

        let mut matrix_ordering_database =
            vec![BuildRunLengthEncoding::new(); 4_usize.pow(args.kmer_length as u32)];
        for (index, (files, _taxid)) in matrix_ordering.iter().progress().enumerate() {
            let bitmap = bitmaps.get(files).unwrap();
            for kmer in bitmap {
                matrix_ordering_database[kmer as usize].push(index);
            }
        }
        let naive_matrix_ordering_runs = matrix_ordering_database
            .iter()
            .map(|build_rle| build_rle.get_vector().len())
            .sum::<usize>();
        let matrix_ordering_compressed_database = matrix_ordering_database
            .into_par_iter()
            .map(|build_rle| build_rle.to_rle())
            .collect::<Vec<RunLengthEncoding>>();
        let matrix_ordering_compressed_runs = matrix_ordering_compressed_database
            .iter()
            .map(|rle| rle.get_vector().len())
            .sum::<usize>();
        drop(matrix_ordering_compressed_database);

        println!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            args.log_blocks,
            block.0,
            block_ordering_compressed_runs,
            matrix_ordering_compressed_runs,
            naive_block_ordering_runs,
            naive_matrix_ordering_runs
        );

        info!("done!");
    }
}
