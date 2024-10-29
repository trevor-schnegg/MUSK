use std::collections::HashSet;

pub fn greedy_ordering(distances: &Vec<Vec<u32>>, start_index: usize) -> Vec<usize> {
    let mut connected_indices = HashSet::from([start_index]);
    let mut ordering = vec![start_index];
    let mut current_index = start_index;

    while ordering.len() < distances.len() {
        let mut next_index = 0_usize;
        let mut next_index_minimum = u32::MAX;

        let mut distance_iter = distances[current_index]
            .iter()
            .chain(
                distances[(current_index + 1)..]
                    .iter()
                    .map(|row| &row[current_index]),
            )
            .enumerate();

        while let Some((index, distance)) = distance_iter.next() {
            if *distance < next_index_minimum && !connected_indices.contains(&index) {
                next_index = index;
                next_index_minimum = *distance;
            }
        }

        ordering.push(next_index);
        connected_indices.insert(next_index);
        current_index = next_index;
    }

    ordering
}

pub fn ordering_statistics(ordering: &Vec<usize>, distances: &Vec<Vec<u32>>) -> (f64, u64) {
    let sum = ordering
        .windows(2)
        .map(|x| {
            if x[0] < x[1] {
                distances[x[1]][x[0]] as u64
            } else {
                distances[x[0]][x[1]] as u64
            }
        })
        .sum();
    (sum as f64 / (ordering.len() - 1) as f64, sum)
}
