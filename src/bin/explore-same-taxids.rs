use clap::Parser;
use log::{debug, info};
use musk::io::load_taxid2files;
use musk::kmer_iter::KmerIter;
use musk::utility::get_fasta_iterator_of_file;
use std::collections::HashSet;
use std::path::Path;
use musk::explore::connected_components;

/// Explores similarities between files with the same species tax id
#[derive(Parser)]
#[clap(version, about)]
#[clap(author = "Trevor S. <trevor.schneggenburger@gmail.com>")]
struct Args {
    #[arg(short, long, default_value_t = 15)]
    /// Length of k-mer to use in the database
    kmer_length: usize,

    #[arg()]
    /// the file2taxid file
    file2taxid: String,
}

fn main() {
    env_logger::init();

    // Parse arguments from the command line
    let args = Args::parse();
    let file2taxid_path = Path::new(&args.file2taxid);

    info!("loading file2taxid at {}", args.file2taxid);
    let file2taxid = load_taxid2files(file2taxid_path);
    info!("file2taxid loaded! exploring files with the same tax id");
    for (taxid, files) in file2taxid {
        if files.len() == 1 {
            println!("{}\t{}", taxid, files[0]);
            continue;
        }
        debug!(
            "creating hashsets for taxid '{}' with {} files...",
            taxid,
            files.len()
        );
        let mut sets = vec![];
        for file in &files {
            let mut record_iter = get_fasta_iterator_of_file(Path::new(file));
            let mut total_kmer_set = HashSet::new();
            while let Some(Ok(record)) = record_iter.next() {
                if record.seq().len() < args.kmer_length {
                    continue;
                }
                let kmers: HashSet<usize> = HashSet::from_iter(KmerIter::from(record.seq(), args.kmer_length));
                total_kmer_set = HashSet::from_iter(total_kmer_set.union(&kmers).map(|x| *x));
            }
            sets.push(total_kmer_set);
        }
        debug!("hashsets created! performing comparisons...");
        let connected_components = connected_components(sets, 0.8);
        for component in connected_components {
            let mut files_string = String::new();
            for file_index in component {
                if files_string.is_empty() {
                    files_string += &*files[file_index];
                } else {
                    files_string += &*(String::from(",") + &*files[file_index])
                }
            }
            println!("{}\t{}", taxid, files_string);
        }
    }
}
