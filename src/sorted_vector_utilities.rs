use std::{iter::Peekable, slice::Iter};
use itertools::{Itertools, KMerge};

pub enum Sequence {
    Single(Vec<u32>),
    Many(Vec<u32>, Vec<Vec<u32>>)
}

pub struct DifferenceIterator<'a> {
    left_iterator: Iter<'a , u32>,
    right_iterator: UnionIterator<'a>,
    next_skip_value: Option<&'a u32>,
}

impl<'a> DifferenceIterator<'a> {
    pub fn from(left_sorted_vector: &'a [u32], right_sorted_vectors: Vec<&'a [u32]>) -> Self {
        let mut right_iterator = UnionIterator::from(right_sorted_vectors);
        let first_skip_value = right_iterator.next();
        DifferenceIterator {
            left_iterator: left_sorted_vector.iter(),
            right_iterator,
            next_skip_value: first_skip_value,
        }
    }
}

impl<'a> Iterator for DifferenceIterator<'a> {
    type Item = &'a u32;

    fn next(&mut self) -> Option<Self::Item> {
        'left_value: while let Some(left_value) = self.left_iterator.next() {
            match self.next_skip_value {
                None => return Some(left_value),
                Some(current_skip_value) => {
                    if *left_value < *current_skip_value {
                        return Some(left_value);
                    } else if *left_value == *current_skip_value {
                        self.next_skip_value = self.right_iterator.next();
                        continue;
                    }
                }
            }
            // At this point, *left_value > *self.next_skip_value.unwrap()
            while let Some(next_skip_value) = self.right_iterator.next() {
                if *next_skip_value > *left_value {
                    self.next_skip_value = Some(next_skip_value);
                    return Some(left_value);
                } else if *next_skip_value == *left_value {
                    self.next_skip_value = self.right_iterator.next();
                    continue 'left_value;
                }
            }
            self.next_skip_value = None;
            return Some(left_value);
        } 
        None
    }
}

pub struct UnionIterator<'a> {
    iterator: Peekable<KMerge<Iter<'a, u32>>>,
}

impl<'a> UnionIterator<'a> {
    pub fn from(sorted_vectors: Vec<&'a [u32]>) -> Self {
        UnionIterator {
            iterator: sorted_vectors.into_iter().kmerge().peekable(),
        }
    }
}

impl<'a> Iterator for UnionIterator<'a> {
    type Item = &'a u32;

    fn next(&mut self) -> Option<Self::Item> {
        let next = match self.iterator.next() {
            None => return None,
            Some(next) => next,
        };
        while let Some(new_next) = self.iterator.peek() {
            if **new_next == *next {
                self.iterator.next();
            } else {
                break;
            }
        }
        Some(next)
    }
}

pub struct IntersectIterator<'a> {
    value_1: &'a u32,
    value_2: &'a u32,
    iterator_1: Iter<'a, u32>,
    iterator_2: Iter<'a, u32>,
}

impl<'a> IntersectIterator<'a> {
    pub fn from(sorted_vector_1: &'a [u32], sorted_vector_2: &'a [u32]) -> Self {
        IntersectIterator {
            value_1: &0,
            value_2: &0,
            iterator_1: sorted_vector_1.iter(),
            iterator_2: sorted_vector_2.iter(),
        }
    }
}

impl<'a> Iterator for IntersectIterator<'a> {
    type Item = &'a u32;

    fn next(&mut self) -> Option<Self::Item> {
        if let (Some(value_1), Some(value_2)) = (self.iterator_1.next(), self.iterator_2.next()) {
            self.value_1 = value_1;
            self.value_2 = value_2;
        } else {
            return None;
        }
        'outer: loop {
            if self.value_1 < self.value_2 {
                while let Some(next_1) = self.iterator_1.next() {
                    if next_1 == self.value_2 {
                        return Some(next_1);
                    } else if next_1 > self.value_2 {
                        self.value_1 = next_1;
                        continue 'outer;
                    }
                }
                return None;
            } else if self.value_2 < self.value_1 {
                while let Some(next_2) = self.iterator_2.next() {
                    if next_2 == self.value_1 {
                        return Some(next_2);
                    } else if next_2 > self.value_2 {
                        self.value_2 = next_2;
                        continue 'outer;
                    }
                }
                return None;
            } else {
                return Some(self.value_1);
            }
        }
    }
}
