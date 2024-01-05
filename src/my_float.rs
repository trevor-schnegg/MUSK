use std::ops::{Add, Div, Mul, MulAssign, Neg, Sub};
use num_traits::{One, Zero};
use crate::decode::{integer_decode_f32, integer_decode_f64};

#[derive(Debug, Clone, Copy, PartialOrd)]
pub struct MyFloat {
    exp: i32,
    float: f32,
}

impl MyFloat {
    pub fn from_f32(f: f32) -> Self {
        let (zeroed_exp_f, exp) = integer_decode_f32(f);
        MyFloat {
            float: zeroed_exp_f,
            exp
        }
    }

    pub fn from_f64(f: f64) -> Self {
        let (zeroed_exp_f, exp) = integer_decode_f64(f);
        MyFloat {
            float: zeroed_exp_f,
            exp
        }
    }

    pub fn ln(&self) -> Self {
        let (zeroed_exp_f, exp) = integer_decode_f32(self.float.ln() + (self.exp as f32 * 2.0_f32.ln()));
        MyFloat {
            float: zeroed_exp_f,
            exp
        }
    }

    pub fn exp(&self) -> Self {
        let base = MyFloat::from_f32(self.float.exp());
        if self.exp.is_positive() {
            let mut acc = base;
            for _ in 0..self.exp {
                acc = acc.square()
            }
            acc
        } else if self.exp.is_negative() {
            let mut acc = base;
            for _ in 0..self.exp.neg() {
                acc = acc.sqrt()
            }
            acc
        } else {
            base
        }
    }

    pub fn sqrt(&self) -> Self {
        if self.exp % 2 == 0 {
            let (zeroed_exp_f, exp) = integer_decode_f32(self.float.sqrt());
            let exp = (self.exp / 2) + exp;
            MyFloat {
                float: zeroed_exp_f,
                exp
            }
        } else {
            let (zeroed_exp_f, exp) = integer_decode_f32(self.float.sqrt() * 2.0_f32.sqrt());
            let exp = ((self.exp - 1) / 2) + exp;
            MyFloat {
                float: zeroed_exp_f,
                exp
            }
        }
    }

    pub fn square(&self) -> Self {
        *self * *self
    }

    pub fn as_f64(&self) -> f64 {
        self.float as f64 * 2.0_f64.powi(self.exp)
    }
}

impl Mul for MyFloat {
    type Output = MyFloat;

    fn mul(self, rhs: Self) -> Self::Output {
        let res = self.float * rhs.float;
        let (zeroed_exp_f, exp) = integer_decode_f32(res);
        MyFloat {
            float: zeroed_exp_f,
            exp: self.exp + rhs.exp + exp,
        }
    }
}

impl MulAssign for MyFloat {
    fn mul_assign(&mut self, rhs: Self) {
        let res = self.float * rhs.float;
        let (zeroed_exp_f, exp) = integer_decode_f32(res);
        self.float = zeroed_exp_f;
        self.exp += rhs.exp + exp;
    }
}

impl Div for MyFloat {
    type Output = MyFloat;

    fn div(self, rhs: Self) -> Self::Output {
        let res = self.float / rhs.float;
        let (zeroed_exp_f, exp) = integer_decode_f32(res);
        MyFloat {
            float: zeroed_exp_f,
            exp: self.exp - rhs.exp + exp,
        }
    }
}

impl Add for MyFloat {
    type Output = MyFloat;

    fn add(self, rhs: Self) -> Self::Output {
        if self.exp == rhs.exp {
            let res = self.float + rhs.float;
            let (zeroed_exp_f, exp) = integer_decode_f32(res);
            MyFloat {
                float: zeroed_exp_f,
                exp: self.exp + exp,
            }
        } else {
            // Need to normalize the exponents so that the floats are comparable
            let left_minus_right = self.exp - rhs.exp;
            if left_minus_right > 0 {
                let res = self.float + (rhs.float * 2.0_f32.powi(left_minus_right.neg()));
                let (zeroed_exp_f, exp) = integer_decode_f32(res);
                MyFloat {
                    float: zeroed_exp_f,
                    exp: self.exp + exp,
                }
            } else { // If left_minus_right < 0
                let res = (self.float * 2.0_f32.powi(left_minus_right)) + rhs.float;
                let (zeroed_exp_f, exp) = integer_decode_f32(res);
                MyFloat {
                    float: zeroed_exp_f,
                    exp: rhs.exp + exp,
                }

            }
        }
    }
}

impl Sub for MyFloat {
    type Output = MyFloat;

    fn sub(self, rhs: Self) -> Self::Output {
        if self.exp == rhs.exp {
            let res = self.float - rhs.float;
            let (zeroed_exp_f, exp) = integer_decode_f32(res);
            MyFloat {
                float: zeroed_exp_f,
                exp: self.exp + exp,
            }
        } else {
            // Need to normalize the exponents so that the floats are comparable
            let left_minus_right = self.exp - rhs.exp;
            if left_minus_right > 0 {
                let res = self.float - (rhs.float * 2.0_f32.powi(left_minus_right.neg()));
                let (zeroed_exp_f, exp) = integer_decode_f32(res);
                MyFloat {
                    float: zeroed_exp_f,
                    exp: self.exp + exp,
                }
            } else { // If left_minus_right < 0
                let res = (self.float * 2.0_f32.powi(left_minus_right)) - rhs.float;
                let (zeroed_exp_f, exp) = integer_decode_f32(res);
                MyFloat {
                    float: zeroed_exp_f,
                    exp: rhs.exp + exp,
                }
            }
        }
    }
}

impl Neg for MyFloat {
    type Output = MyFloat;

    fn neg(self) -> Self::Output {
        MyFloat {
            exp: self.exp,
            float: self.float.neg(),
        }
    }
}

impl One for MyFloat {
    fn one() -> Self {
        MyFloat::from_f32(1.0)
    }
}

impl Zero for MyFloat {
    fn zero() -> Self {
        MyFloat::from_f32(0.0)
    }

    fn is_zero(&self) -> bool {
        if self.exp == 0 && self.float == 0.0 {
            true
        } else {
            false
        }
    }
}

impl PartialEq<Self> for MyFloat {
    fn eq(&self, other: &Self) -> bool {
        if self.float == other.float && self.exp == other.exp {
            true
        } else {
            false
        }
    }
}
