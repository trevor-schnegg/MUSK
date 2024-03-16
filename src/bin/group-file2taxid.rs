use clap::Parser;
use log::{debug, info};
use musk::explore::connected_components;
use musk::io::load_taxid2files;
use musk::kmer_iter::KmerIter;
use musk::utility::get_fasta_iterator_of_file;
use roaring::RoaringBitmap;
use std::path::Path;
use std::sync::mpsc;
use threadpool::ThreadPool;

fn create_bitmaps(
    files: &Vec<String>,
    kmer_length: usize,
    thread_number: usize,
) -> Vec<RoaringBitmap> {
    let mut bitmaps = vec![];
    let (sender, receiver) = mpsc::channel();
    let pool = ThreadPool::new(thread_number);
    for file in files.clone() {
        let sender_clone = sender.clone();
        pool.execute(move || {
            let mut bitmap = RoaringBitmap::new();
            let mut record_iter = get_fasta_iterator_of_file(Path::new(&file));
            while let Some(Ok(record)) = record_iter.next() {
                if record.seq().len() < kmer_length {
                    continue;
                }
                for kmer in KmerIter::from(record.seq(), kmer_length) {
                    bitmap.insert(kmer as u32);
                }
            }
            sender_clone.send(bitmap).unwrap();
        });
    }
    drop(sender);
    for bitmap in receiver {
        bitmaps.push(bitmap);
    }
    bitmaps
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
    taxid2file: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let file2taxid_path = Path::new(&args.taxid2file);

    info!("loading file2taxid at {}", args.taxid2file);
    let file2taxid = load_taxid2files(file2taxid_path);
    info!("file2taxid loaded! exploring files with the same tax id");
    for (taxid, files) in file2taxid {
        if files.len() == 1 {
            println!("{}\t{}", files[0], taxid);
            continue;
        }
        debug!(
            "creating sets for taxid '{}' with {} files...",
            taxid,
            files.len()
        );
        let bitmaps = create_bitmaps(&files, args.kmer_length, args.thread_number);
        debug!("hashsets created! performing comparisons...");
        let connected_components = connected_components(bitmaps, 0.8, args.thread_number);
        for component in connected_components {
            let mut files_string = String::new();
            for file_index in component {
                if files_string.is_empty() {
                    files_string += &*files[file_index];
                } else {
                    files_string += &*(String::from(",") + &*files[file_index])
                }
            }
            println!("{}\t{}", files_string, taxid);
        }
    }
}
