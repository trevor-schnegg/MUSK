use std::sync::Arc;
use std::sync::mpsc::Sender;
use std::collections::HashMap;
use crate::sorted_vector_utilities::{DifferenceIterator, IntersectIterator};

pub enum Sequence {
    One(Vec<u32>, String, u32),
    Many(Vec<u32>, Vec<(Vec<u32>, String)>, u32),
}

fn distance(length_1: usize, length_2: usize, intersection_size: usize) -> u32 {
    (length_1 + length_2 - (2 * intersection_size)) as u32
}

pub fn self_matrix(
    many_sequences: (&Vec<u32>, &Vec<(Vec<u32>, String)>),
    sender: &Sender<(usize, usize, u32)>,
    file_to_index: &Arc<HashMap<String, usize>>,
) -> () {
    let (union, difference_vectors) = many_sequences;
    for index_1 in 0..difference_vectors.len() {
        for index_2 in 0..difference_vectors.len() {
            if index_2 <= index_1 {
                continue;
            }
            let (difference_1, difference_2) =
                (&difference_vectors[index_1], &difference_vectors[index_2]);
            let intersection_size =
                DifferenceIterator::from(union, vec![&difference_1.0, &difference_2.0]).count();
            let distance = distance(
                union.len() - difference_1.0.len(),
                union.len() - difference_2.0.len(),
                intersection_size,
            );
            let (sequence_index_1, sequence_index_2) = (
                *file_to_index.get(&difference_1.1).unwrap(),
                *file_to_index.get(&difference_2.1).unwrap(),
            );
            sender
                .send((sequence_index_1, sequence_index_2, distance))
                .unwrap();
        }
    }
}

pub fn one_to_one(
    sequence_1: (&Vec<u32>, &String),
    sequence_2: (&Vec<u32>, &String),
    sender: &Sender<(usize, usize, u32)>,
    file_to_index: &Arc<HashMap<String, usize>>,
) -> () {
    let intersection_size = IntersectIterator::from(&sequence_1.0, &sequence_2.0).count();
    let distance = distance(sequence_1.0.len(), sequence_2.0.len(), intersection_size);
    let (sequence_index_1, sequence_index_2) = (
        *file_to_index.get(sequence_1.1).unwrap(),
        *file_to_index.get(sequence_2.1).unwrap(),
    );
    sender
        .send((sequence_index_1, sequence_index_2, distance))
        .unwrap();
}

pub fn many_to_one(
    many_sequences: (&Vec<u32>, &Vec<(Vec<u32>, String)>),
    one: (&Vec<u32>, &String),
    sender: &Sender<(usize, usize, u32)>,
    file_to_index: &Arc<HashMap<String, usize>>,
) -> () {
    let (union, differences) = many_sequences;
    let intersection = IntersectIterator::from(union, one.0)
        .map(|kmer| *kmer)
        .collect::<Vec<u32>>();
    for difference in differences {
        let intersection_size =
            DifferenceIterator::from(&intersection, vec![&difference.0]).count();
        let distance = distance(
            union.len() - difference.0.len(),
            one.0.len(),
            intersection_size,
        );
        let (sequence_index_1, sequence_index_2) = (
            *file_to_index.get(one.1).unwrap(),
            *file_to_index.get(&difference.1).unwrap(),
        );
        sender
            .send((sequence_index_1, sequence_index_2, distance))
            .unwrap();
    }
}

pub fn many_to_many(
    many_sequences_1: (&Vec<u32>, &Vec<(Vec<u32>, String)>),
    many_sequences_2: (&Vec<u32>, &Vec<(Vec<u32>, String)>),
    sender: &Sender<(usize, usize, u32)>,
    file_to_index: &Arc<HashMap<String, usize>>,
) -> () {
    let (union_1, differences_1) = many_sequences_1;
    let (union_2, differences_2) = many_sequences_2;
    let intersection = IntersectIterator::from(union_1, union_2)
        .map(|kmer| *kmer)
        .collect::<Vec<u32>>();
    for difference_1 in differences_1 {
        for difference_2 in differences_2 {
            let intersection_size =
                DifferenceIterator::from(&intersection, vec![&difference_1.0, &difference_2.0])
                    .count();
            let distance = distance(
                union_1.len() - difference_1.0.len(),
                union_2.len() - difference_2.0.len(),
                intersection_size,
            );
            let (sequence_index_1, sequence_index_2) = (
                *file_to_index.get(&difference_1.1).unwrap(),
                *file_to_index.get(&difference_2.1).unwrap(),
            );
            sender
                .send((sequence_index_1, sequence_index_2, distance))
                .unwrap();
        }
    }
}
