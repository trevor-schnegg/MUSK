use musk::my_float::{integer_decode_f32, MyFloat};

fn main() {
    // let f1: f32 = 1.0e-12;
    // println!("{:032b}", f1.to_bits());
    // println!("{:032b}", f1.to_bits() << 1);
    // // println!("{:032b}", f1.to_bits() << 3);
    // println!("{:08b}", (f1.to_bits() << 1).to_be_bytes()[0]);
    // println!("{}", get_exp_f32(f1));
    // println!("{:032b}", set_exp_to_zero(f1).to_bits());
    // println!("{}", get_exp_f32(set_exp_to_zero(f1)));

    // let f2: f32 = 1.0e-12;
    // let f3: f32 = 1.0e20;
    // let mut my_f2 = MyFloat::new(f2);
    // my_f2.multiply(f3);
    // println!("{:?}", my_f2)

    let f4: f32 = 0.5;
    let f5: f32 = 1.0 * 2_usize.pow(24) as f32;
    println!("{}", f4);
    println!("{:?}", integer_decode_f32(f4));
    println!("{}", f5);
    println!("{}", f4 + f5);
    println!("{:?}", integer_decode_f32(f4 / f5));
    println!("{:?}", MyFloat::new(f4) * MyFloat::new(f5));
}