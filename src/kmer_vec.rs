use crate::kmer_vec::BinarySearchResult::{ExistingVectorAt, IndexToInsert};

enum BinarySearchResult {
    IndexToInsert(usize),
    ExistingVectorAt(usize),
}

pub(crate) struct KmerVec {
    data: Vec<(u32, Vec<usize>)>,
}

impl KmerVec {
    pub fn new() -> Self {
        KmerVec { data: Vec::new() }
    }

    pub fn insert(&mut self, kmer: u32, accession_index: usize) -> bool {
        if self.data.len() == 0 {
            self.data.push((kmer, Vec::from(vec![accession_index])));
            return true;
        }
        match self.binary_search(kmer) {
            IndexToInsert(index) => {
                self.data
                    .insert(index, (kmer, Vec::from(vec![accession_index])));
            }
            ExistingVectorAt(index) => {
                let existing_vec = &mut self.data.get_mut(index).unwrap().1;
                if existing_vec.contains(&accession_index) {
                    return false;
                } else {
                    existing_vec.push(accession_index)
                }
            }
        }
        true
    }

    pub fn query(&self, kmer: u32) -> Option<&Vec<usize>> {
        match self.binary_search(kmer) {
            IndexToInsert(_) => None,
            ExistingVectorAt(index) => Some(&self.data.get(index).unwrap().1),
        }
    }

    fn binary_search(&self, kmer: u32) -> BinarySearchResult {
        self.binary_search_helper(kmer, 0, self.data.len() - 1)
    }

    fn binary_search_helper(
        &self,
        search_kmer: u32,
        low: usize,
        high: usize,
    ) -> BinarySearchResult {
        let mid = (low + high) / 2;
        let (kmer, _) = self.data.get(mid).unwrap();
        // If the value was found
        if search_kmer == *kmer {
            ExistingVectorAt(mid)
        }
        // If we need to insert at either end of the data vector
        else if mid == 0 && search_kmer < *kmer {
            IndexToInsert(0)
        } else if mid == (self.data.len() - 1) && search_kmer > *kmer {
            IndexToInsert(self.data.len())
        }
        // If we need to insert somewhere in the middle
        else if search_kmer < *kmer && search_kmer > self.data.get(mid - 1).unwrap().0 {
            IndexToInsert(mid)
        } else if search_kmer > *kmer && search_kmer < self.data.get(mid + 1).unwrap().0 {
            IndexToInsert(mid + 1)
        }
        // Otherwise, recurse appropriately
        else if search_kmer < *kmer {
            self.binary_search_helper(search_kmer, low, mid - 1)
        } else {
            self.binary_search_helper(search_kmer, mid + 1, high)
        }
    }
}
