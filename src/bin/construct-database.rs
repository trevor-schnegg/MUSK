use clap::Parser;
use itertools::Itertools;
use log::{debug, info};
use musk::io::{dump_data_to_file, load_data_from_file};
use musk::kmer_iter::KmerIter;
use musk::rle::{BuildRunLengthEncoding, RunLengthEncoding};
use musk::utility::get_fasta_iterator_of_file;
use roaring::RoaringBitmap;
use std::path::Path;
use rayon::prelude::*;

// 999,999,937
const XOR_NUMBER: usize = 0b_00111011100110101100100111000001;

fn create_bitmap(files: &str, kmer_length: usize, low: usize, high: usize) -> RoaringBitmap {
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
    bitmap
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
    /// the file2taxid file
    ordering_file: String,

    #[arg()]
    /// Location to output the serialzed distances
    output_file: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let ordering_file_path = Path::new(&args.ordering_file);
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

    info!("loading ordering at {}", args.ordering_file);
    let mut ordering = load_data_from_file::<Vec<(String, u32)>>(ordering_file_path);
    if let (Some(old_prefix), Some(new_prefix))= (args.old_directory_prefix, args.new_directory_prefix) {
        ordering = ordering.into_iter().map(|(files, taxid)| (files.replace(&*old_prefix, &*new_prefix), taxid)).collect_vec();
    }
    info!("creating sorted kmer vectors for each group...");
    let bitmaps = ordering.par_iter().map(|(files, _taxid)| {
        create_bitmap(files, args.kmer_length, lowest_kmer, highest_kmer)

    }).collect::<Vec<RoaringBitmap>>();
    info!("kmer vectors computed, creating database...");

    let mut database = vec![BuildRunLengthEncoding::new(); 4_usize.pow(args.kmer_length as u32)];
    for (index, bitmap) in bitmaps.into_iter().enumerate() {
        for kmer in bitmap {
            database[kmer as usize].push(index);
        }
        if index % 1000 == 0 && index != 0 {
            debug!("done inserting {} bitmaps into the database", index);
        }
    }
    let naive_runs = database
        .iter()
        .map(|build_rle| build_rle.get_vector().len())
        .sum::<usize>();
    info!("Total naive runs for the ordering {}", naive_runs);
    let compressed_database = database
        .into_par_iter()
        .map(|build_rle| build_rle.to_rle())
        .collect::<Vec<RunLengthEncoding>>();
    let compressed_runs = compressed_database
        .iter()
        .map(|rle| rle.get_vector().len())
        .sum::<usize>();
    info!("Total compressed runs for the ordering {}", compressed_runs);
    info!("Saving the compressed runs to the output file...");
    dump_data_to_file(
        bincode::serialize(&compressed_database).unwrap(),
        output_file_path,
    )
    .unwrap();
}
