fn geometric_search(sorted_vector: &[u32], value: u32) -> Result<usize, usize> {
    let mut current_index = 0;
    let mut last_index = 0;
    while sorted_vector[current_index] < value {
        last_index = current_index;
        current_index = if current_index == 0 {
            1
        } else {
            current_index << 1
        };
        if current_index >= sorted_vector.len() {
            current_index = sorted_vector.len();
            break;
        }
    }
    let result = if current_index == sorted_vector.len() {
        sorted_vector[last_index..].binary_search(&value)
    } else if sorted_vector[current_index] == value {
        return Ok(current_index);
    } else {
        sorted_vector[last_index..current_index].binary_search(&value)
    };
    match result {
        Ok(index) => Ok(last_index + index),
        Err(index) => Err(last_index + index),
    }
}

pub struct IntersectIterator<'a> {
    position_1: usize,
    position_2: usize,
    value_1: u32,
    value_2: u32,
    sorted_vector_1: &'a [u32],
    sorted_vector_2: &'a [u32],
    initialized: bool,
    at_end: bool,
}

impl<'a> IntersectIterator<'a> {
    pub fn from(sorted_vector_1: &'a [u32], sorted_vector_2: &'a [u32]) -> Self {
        IntersectIterator {
            position_1: 0,
            position_2: 0,
            value_1: 0,
            value_2: 0,
            sorted_vector_1,
            sorted_vector_2,
            initialized: false,
            at_end: false,
        }
    }

    fn search(&mut self) -> Option<u32> {
        loop {
            if self.value_1 < self.value_2 {
                match geometric_search(&self.sorted_vector_1[self.position_1..], self.value_2) {
                    Ok(index) => {
                        self.position_1 += index;
                        self.value_1 = self.sorted_vector_1[self.position_1];
                        assert_eq!(self.value_1, self.value_2);
                        return Some(self.value_1);
                    }
                    Err(index) => {
                        self.position_1 += index;
                        if self.position_1 >= self.sorted_vector_1.len() {
                            self.at_end = true;
                            return None;
                        }
                        self.value_1 = self.sorted_vector_1[self.position_1];
                    }
                }
            } else if self.value_2 < self.value_1 {
                match geometric_search(&self.sorted_vector_2[self.position_2..], self.value_1) {
                    Ok(index) => {
                        self.position_2 += index;
                        self.value_2 = self.sorted_vector_2[self.position_2];
                        assert_eq!(self.value_1, self.value_2);
                        return Some(self.value_2);
                    }
                    Err(index) => {
                        self.position_2 += index;
                        if self.position_2 >= self.sorted_vector_2.len() {
                            self.at_end = true;
                            return None;
                        }
                        self.value_2 = self.sorted_vector_2[self.position_2];
                    }
                }
            } else {
                panic!("impossible case")
            }
        }
    }
}

impl<'a> Iterator for IntersectIterator<'a> {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        if self.at_end {
            return None;
        }

        if !self.initialized {
            self.value_1 = self.sorted_vector_1[self.position_1];
            self.value_2 = self.sorted_vector_2[self.position_2];
            self.initialized = true;
            if self.value_1 == self.value_2 {
                return Some(self.value_1);
            } else {
                self.search()
            }
        } else {
            self.position_1 += 1;
            if self.position_1 >= self.sorted_vector_1.len() {
                self.at_end = true;
                return None;
            } else {
                self.value_1 = self.sorted_vector_1[self.position_1];
            }
            self.search()
        }
    }
}
