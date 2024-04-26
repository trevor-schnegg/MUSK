use clap::Parser;
use itertools::Itertools;
use log::{debug, info};
use musk::io::{dump_data_to_file, load_data_from_file};
use musk::kmer_iter::KmerIter;
use musk::rle::BuildRunLengthEncoding;
use musk::utility::get_fasta_iterator_of_file;
use roaring::RoaringBitmap;
use std::collections::HashMap;
use std::path::Path;
use std::sync::mpsc;
use threadpool::ThreadPool;

fn create_bitmap(files: &str, kmer_length: usize) -> RoaringBitmap {
    let mut bitset = RoaringBitmap::new();
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
    bitset
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

    info!("loading ordering at {}", args.ordering_file);
    let ordering = load_data_from_file::<Vec<(String, u32)>>(ordering_file_path);
    info!("creating sorted kmer vectors for each group...");
    let mut bitmaps = HashMap::new();
    let (sender, receiver) = mpsc::channel();
    let pool = ThreadPool::new(args.thread_number);
    for (files, _) in ordering.clone() {
        let sender_clone = sender.clone();
        pool.execute(move || {
            let bitmap = create_bitmap(&*files, args.kmer_length);
            sender_clone.send((files, bitmap)).unwrap();
        })
    }
    drop(sender);
    for (files, sorted_kmer_vector) in receiver {
        bitmaps.insert(files, sorted_kmer_vector);
    }
    info!("kmer vectors computed, creating database...");

    let mut database = vec![BuildRunLengthEncoding::new(); 4_usize.pow(args.kmer_length as u32)];
    for (index, (files, _taxid)) in ordering.into_iter().enumerate() {
        for kmer in bitmaps.get(&files).unwrap() {
            database[kmer as usize].push(index);
        }
        if index % 1000 == 0 && index != 0 {
            debug!("done inserting {} bitmaps into the database", index);
        }
    }
    let naive_runs = database.iter().map(|build_rle| build_rle.get_vector().len()).sum::<usize>();
    info!("Total naive runs for the ordering {}", naive_runs);
    let compressed_database = database.into_iter().map(|build_rle| build_rle.to_rle()).collect_vec();
    let compressed_runs = compressed_database.iter().map(|rle| rle.get_vector().len()).sum::<usize>();
    info!("Total compressed runs for the ordering {}", compressed_runs);
    info!("Saving the compressed runs to the output file...");
    dump_data_to_file(bincode::serialize(&compressed_database).unwrap(), output_file_path).unwrap();
}
