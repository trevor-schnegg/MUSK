#[derive(Clone)]
pub struct HitCounter {
    last_was_hit: bool,
    match2match: usize,
    match2nonmatch: usize,
    nonmatch2match: usize,
    nonmatch2nonmatch: usize,
}

impl HitCounter {
    pub fn new() -> Self {
        HitCounter {
            last_was_hit: false,
            match2match: 0,
            match2nonmatch: 0,
            nonmatch2match: 0,
            nonmatch2nonmatch: 0,
        }
    }

    pub fn hit(&mut self) -> () {
        if self.last_was_hit {
            self.match2match += 1;
        } else {
            self.nonmatch2match += 1;
            self.last_was_hit = true;
        }
    }

    pub fn miss(&mut self) -> () {
        if self.last_was_hit {
            self.match2nonmatch += 1;
            self.last_was_hit = false;
        } else {
            self.nonmatch2nonmatch += 1;
        }
    }

    pub fn get_counts(&self) -> (usize, usize, usize, usize) {
        (self.match2match, self.match2nonmatch, self.nonmatch2match, self.nonmatch2nonmatch)
    }
}