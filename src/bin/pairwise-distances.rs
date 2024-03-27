use clap::Parser;
use log::info;
use musk::io::{dump_data_to_file, load_string2taxid};
use musk::kmer_iter::KmerIter;
use musk::utility::get_fasta_iterator_of_file;
use roaring::RoaringBitmap;
use std::cmp::Reverse;
use std::path::Path;
use std::sync::{mpsc, Arc};
use threadpool::ThreadPool;

fn create_bitmap(files: &str, kmer_length: usize) -> RoaringBitmap {
    let mut bitmap = RoaringBitmap::new();
    for file in files.split(",") {
        let mut record_iter = get_fasta_iterator_of_file(Path::new(&file));
        while let Some(Ok(record)) = record_iter.next() {
            if record.seq().len() < kmer_length {
                continue;
            }
            for kmer in KmerIter::from(record.seq(), kmer_length) {
                bitmap.insert(kmer as u32);
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
    let mut bitmaps = vec![];
    let (sender, receiver) = mpsc::channel();
    let pool = ThreadPool::new(args.thread_number);
    for (files, taxid) in file2taxid {
        let sender_clone = sender.clone();
        pool.execute(move || {
            let bitmap = create_bitmap(&*files, args.kmer_length);
            sender_clone.send((bitmap, files, taxid)).unwrap();
        })
    }
    drop(sender);
    for triple in receiver {
        bitmaps.push(triple);
    }

    info!("bitmaps computed, computing pairwise distances...");

    let mut all_distances = vec![];
    let bitmaps_arc = Arc::new(bitmaps);
    let (sender, receiver) = mpsc::channel();
    for index_1 in 0..bitmaps_arc.len() {
        let bitmaps_arc_clone = bitmaps_arc.clone();
        let sender_clone = sender.clone();
        pool.execute(move || {
            let mut distances = vec![];
            for index_2 in 0..bitmaps_arc_clone.len() {
                if index_2 <= index_1 {
                    continue;
                }
                let (bitmap_1, bitmap_2) =
                    (&bitmaps_arc_clone[index_1].0, &bitmaps_arc_clone[index_2].0);
                let intersection_size = bitmap_1.intersection_len(bitmap_2);
                // |A| + |B| - 2 * |A and B|
                distances.push(bitmap_1.len() + bitmap_2.len() - (2 * intersection_size));
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
    drop(bitmaps_arc);
    for distances in receiver {
        all_distances.push(distances);
    }
    all_distances.sort_by_key(|x| Reverse(x.0.len()));

    dump_data_to_file(
        bincode::serialize(&all_distances).unwrap(),
        output_file_path,
    )
    .unwrap();
}
