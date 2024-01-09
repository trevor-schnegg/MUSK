use crate::big_exp_float::BigExpFloat;
use crate::consts::Consts;
use approx::ulps_eq;
use num_traits::{FloatConst, One, Zero};
use statrs::StatsError;

pub fn sf(p: f64, n: u64, x: u64, consts: &Consts) -> BigExpFloat {
    if x >= n {
        BigExpFloat::zero()
    } else {
        let k = x;
        checked_beta_reg(k as f64 + 1.0, (n - k) as f64, p, consts).unwrap()
    }
}

fn checked_beta_reg(a: f64, b: f64, x: f64, consts: &Consts) -> Result<BigExpFloat, StatsError> {
    if a <= 0.0 {
        Err(StatsError::ArgMustBePositive("a"))
    } else if b <= 0.0 {
        Err(StatsError::ArgMustBePositive("b"))
    } else if !(0.0..=1.0).contains(&x) {
        Err(StatsError::ArgIntervalIncl("x", 0.0, 1.0))
    } else {
        let bt = if x.is_zero() || ulps_eq!(x, 1.0) {
            BigExpFloat::zero()
        } else {
            (ln_gamma(a + b, consts) - ln_gamma(a, consts) - ln_gamma(b, consts)
                + BigExpFloat::from_f64(a * x.ln())
                + BigExpFloat::from_f64(b * (1.0 - x).ln()))
            .exp()
        };
        let symm_transform = x >= (a + 1.0) / (a + b + 2.0);

        let mut a = BigExpFloat::from_f64(a);
        let mut b = BigExpFloat::from_f64(b);
        let mut x = BigExpFloat::from_f64(x);
        if symm_transform {
            let swap = a;
            x = BigExpFloat::one() - x;
            a = b;
            b = swap;
        }

        let qab = a + b;
        let qap = a + BigExpFloat::one();
        let qam = a - BigExpFloat::one();
        let mut c = BigExpFloat::one();
        let mut d = BigExpFloat::one() - qab * x / qap;

        d = BigExpFloat::one() / d;
        let mut h = d;

        for m in 1..141 {
            let m = BigExpFloat::from_f32(m as f32);
            let m2 = m * BigExpFloat::from_f64(2.0);
            let mut aa = m * (b - m) * x / ((qam + m2) * (a + m2));
            d = BigExpFloat::one() + aa * d;

            c = BigExpFloat::one() + aa / c;

            d = BigExpFloat::one() / d;
            h = h * d * c;
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
            d = BigExpFloat::one() + aa * d;

            c = BigExpFloat::one() + aa / c;

            d = BigExpFloat::one() / d;
            let del = d * c;
            h *= del;
        }

        if symm_transform {
            Ok(BigExpFloat::one() - bt * h / a)
        } else {
            Ok(bt * h / a)
        }
    }
}

fn ln_gamma(x: f64, consts: &Consts) -> BigExpFloat {
    if x < 0.5 {
        let s = consts
            .gamma_dk
            .iter()
            .enumerate()
            .skip(1)
            .fold(consts.gamma_dk[0], |s, t| {
                s + *t.1 / (BigExpFloat::from_f32(t.0 as f32) - BigExpFloat::from_f64(x))
            });

        consts.ln_pi
            - BigExpFloat::from_f64((f64::PI() * x).sin().ln())
            - s.ln()
            - consts.ln_2_sqrt_e_over_pi
            - BigExpFloat::from_f64((0.5 - x) * ((0.5 - x + consts.gamma_r) / f64::E()).ln())
    } else {
        let s = consts
            .gamma_dk
            .iter()
            .enumerate()
            .skip(1)
            .fold(consts.gamma_dk[0], |s, t| {
                s + *t.1
                    / (BigExpFloat::from_f64(x) + BigExpFloat::from_f32(t.0 as f32)
                        - BigExpFloat::one())
            });

        s.ln()
            + consts.ln_2_sqrt_e_over_pi
            + BigExpFloat::from_f64((x - 0.5) * ((x - 0.5 + consts.gamma_r) / f64::E()).ln())
    }
}
