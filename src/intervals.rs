pub struct Intervals {
    vec: Vec<(usize, usize, String)>
}

impl Intervals {
    pub fn new() -> Self {
        Intervals {
            vec: Vec::new()
        }
    }

    pub fn push(&mut self, val: (usize, usize, String)) -> () {
        self.vec.push(val)
    }

    pub fn get_accession_of_index(&self, index: usize) -> String {
        self.binary_search(index, 0, self.vec.len() - 1)
    }

    fn binary_search(&self, index: usize, low: usize, high: usize) -> String {
        let mid = (low + high) / 2;
        let (start, end, accession) = self.vec.get(mid).unwrap();
        if index >= *start && index <= *end {
            accession.clone()
        } else if index < *start {
            self.binary_search(index, low, mid - 1)
        } else {
            self.binary_search(index, mid + 1, high)
        }
    }
}