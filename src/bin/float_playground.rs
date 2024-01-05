use musk::my_float::MyFloat;

fn main() {
    let f1 = 0.19873465345e-11_f64;
    let my_f1 = MyFloat::from_f64(f1);
    println!("{}, {}", f1.sqrt(), my_f1.sqrt().as_f64());
    println!("{}, {}", f1.powi(2), my_f1.square().as_f64());

    let f2= 0.000000000000000571238572_f64;
    println!("{}", f2.exp() as f32);
    let my_f2 = MyFloat::from_f64(f2);
    println!("{}, {}", f2.exp(), my_f2.exp().as_f64());
}