use std::slice::Iter;

pub struct IntersectIterator<'a> {
    value_1: u32,
    value_2: u32,
    iterator_1: Iter<'a, u32>,
    iterator_2: Iter<'a, u32>,
    initialized: bool,
}

impl<'a> IntersectIterator<'a> {
    pub fn from(sorted_vector_1: &'a [u32], sorted_vector_2: &'a [u32]) -> Self {
        IntersectIterator {
            value_1: 0,
            value_2: 0,
            iterator_1: sorted_vector_1.iter(),
            iterator_2: sorted_vector_2.iter(),
            initialized: false,
        }
    }

    fn search(&mut self) -> Option<u32> {
        loop {
            if self.value_1 < self.value_2 {
                // Advance iterator 1 until the value is >= value_2
                loop {
                    match self.iterator_1.next() {
                        None => return None,
                        Some(next_1) => {
                            if *next_1 >= self.value_2 {
                                self.value_1 = *next_1;
                                break;
                            }
                        }
                    }
                }
            } else if self.value_2 < self.value_1 {
                // Advance iterator_2 until the value is >= value_1
                loop {
                    match self.iterator_2.next() {
                        None => return None,
                        Some(next_2) => {
                            if *next_2 >= self.value_1 {
                                self.value_2 = *next_2;
                                break;
                            }
                        }
                    }
                }
            } else {
                // the values must be equal
                return Some(self.value_1);
            }
        }
    }
}

impl<'a> Iterator for IntersectIterator<'a> {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.initialized {
            self.value_1 = match self.iterator_1.next() {
                None => return None,
                Some(x) => *x,
            };
            self.value_2 = match self.iterator_2.next() {
                None => return None,
                Some(x) => *x,
            };
            self.initialized = true;
            if self.value_1 == self.value_2 {
                return Some(self.value_1);
            }
            self.search()
        } else {
            self.value_1 = match self.iterator_1.next() {
                None => return None,
                Some(x) => *x,
            };
            self.search()
        }
    }
}
