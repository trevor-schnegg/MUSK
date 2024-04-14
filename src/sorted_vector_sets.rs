use std::slice::Iter;

pub struct IntersectIterator<'a> {
    value_1: u32,
    value_2: u32,
    iterator_1: Iter<'a, u32>,
    iterator_2: Iter<'a, u32>,
}

impl<'a> IntersectIterator<'a> {
    pub fn from(sorted_vector_1: &'a [u32], sorted_vector_2: &'a [u32]) -> Self {
        IntersectIterator {
            value_1: 0,
            value_2: 0,
            iterator_1: sorted_vector_1.iter(),
            iterator_2: sorted_vector_2.iter(),
        }
    }
}

impl<'a> Iterator for IntersectIterator<'a> {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        if let (Some(value_1), Some(value_2)) = (self.iterator_1.next(), self.iterator_2.next()) {
            self.value_1 = *value_1;
            self.value_2 = *value_2;
        } else {
            return None;
        }
        loop {
            if self.value_1 < self.value_2 {
                let value_1_at_start = self.value_1;
                while let Some(next_1) = self.iterator_1.next() {
                    if *next_1 >= self.value_2 {
                        self.value_1 = *next_1;
                        if *next_1 == self.value_2 {
                            return Some(*next_1);
                        }
                        break;
                    }
                }
                if value_1_at_start == self.value_1 {
                    return None;
                }
            } else if self.value_2 < self.value_1 {
                let value_2_at_start = self.value_2;
                while let Some(next_2) = self.iterator_2.next() {
                    if *next_2 >= self.value_1 {
                        self.value_2 = *next_2;
                        if *next_2 == self.value_1 {
                            return Some(*next_2);
                        }
                        break;
                    }
                }
                if value_2_at_start == self.value_2 {
                    return None;
                }
            } else {
                return Some(self.value_1);
            }
        }
    }
}
