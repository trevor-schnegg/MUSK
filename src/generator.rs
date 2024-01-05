use rand::distributions::Uniform;
use rand::{thread_rng, Rng};

fn u8_to_base(n: u8) -> char {
    if n == 0 {
        'A'
    } else if n == 1 {
        'C'
    } else if n == 2 {
        'G'
    } else if n == 3 {
        'T'
    } else {
        panic!("Impossible case")
    }
}

pub fn create_random_read(len: usize) -> String {
    let rng = thread_rng();
    let range = Uniform::from(0..4);
    rng.sample_iter(&range)
        .take(len)
        .map(|n| u8_to_base(n))
        .collect()
}
