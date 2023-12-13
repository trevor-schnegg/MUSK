use num_traits::{FloatConst, One, Zero};
use statrs::StatsError;
use crate::my_float::MyFloat;
use approx::ulps_eq;

struct Consts {
    gamma_r: f64,
    gamma_dk: Vec<MyFloat>,
    ln_pi: MyFloat,
    ln_2_sqrt_e_over_pi: MyFloat,
}

impl Consts {
    fn new() -> Self {
        Consts {
            gamma_r: 10.900511,
            gamma_dk: vec![MyFloat::from_f64(2.48574089138753565546e-5),
                           MyFloat::from_f64(1.05142378581721974210),
                           MyFloat::from_f64(-3.45687097222016235469),
                           MyFloat::from_f64(4.51227709466894823700),
                           MyFloat::from_f64(-2.98285225323576655721),
                           MyFloat::from_f64(1.05639711577126713077),
                           MyFloat::from_f64(-1.95428773191645869583e-1),
                           MyFloat::from_f64(1.70970543404441224307e-2),
                           MyFloat::from_f64(-5.71926117404305781283e-4),
                           MyFloat::from_f64(4.63399473359905636708e-6),
                           MyFloat::from_f64(-2.71994908488607703910e-9),
            ],
            ln_pi: MyFloat::from_f64(1.1447298858494001741434273513530587116472948129153),
            ln_2_sqrt_e_over_pi: (MyFloat::from_f64(0.6207822376352452223455184457816472122518527279025978)),
        }
    }
}

pub fn sf(p: f64, n: u64, x: u64) -> MyFloat {
    if x >= n {
        MyFloat::zero()
    } else {
        let k = x;
        checked_beta_reg(k as f64 + 1.0, (n - k) as f64, p).unwrap()
    }
}

fn checked_beta_reg(a: f64, b: f64, x: f64) -> Result<MyFloat, StatsError> {
    if a <= 0.0 {
        Err(StatsError::ArgMustBePositive("a"))
    } else if b <= 0.0 {
        Err(StatsError::ArgMustBePositive("b"))
    } else if !(0.0..=1.0).contains(&x) {
        Err(StatsError::ArgIntervalIncl("x", 0.0, 1.0))
    } else {
        let bt = if x.is_zero() || ulps_eq!(x, 1.0) {
            MyFloat::zero()
        } else {
            (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b)
                + MyFloat::from_f64(a * x.ln())
                + MyFloat::from_f64(b * (1.0 - x).ln()))
                .exp()
        };
        let symm_transform = x >= (a + 1.0) / (a + b + 2.0);

        let mut a = MyFloat::from_f64(a);
        let mut b = MyFloat::from_f64(b);
        let mut x = MyFloat::from_f64(x);
        if symm_transform {
            let swap = a;
            x = MyFloat::one() - x;
            a = b;
            b = swap;
        }

        let qab = a + b;
        let qap = a + MyFloat::one();
        let qam = a - MyFloat::one();
        let mut c = MyFloat::one();
        let mut d = MyFloat::one() - qab * x / qap;

        d = MyFloat::one() / d;
        let mut h = d;

        for m in 1..141 {
            let m = MyFloat::from_f64(m as f64);
            let m2 = m * MyFloat::from_f64(2.0);
            let mut aa = m * (b - m) * x / ((qam + m2) * (a + m2));
            d = MyFloat::one() + aa * d;

            c = MyFloat::one() + aa / c;

            d = MyFloat::one() / d;
            h = h * d * c;
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
            d = MyFloat::one() + aa * d;

            c = MyFloat::one() + aa / c;

            d = MyFloat::one() / d;
            let del = d * c;
            h *= del;
        }

        if symm_transform {
            Ok(MyFloat::one() - bt * h / a)
        } else {
            Ok(bt * h / a)
        }
    }
}

fn ln_gamma(x: f64) -> MyFloat {
    let consts = Consts::new();
    if x < 0.5 {
        let s = consts.gamma_dk
            .iter()
            .enumerate()
            .skip(1)
            .fold(consts.gamma_dk[0], |s, t| s + *t.1 / (MyFloat::from_f32(t.0 as f32) - MyFloat::from_f64(x)));

        consts.ln_pi
            - MyFloat::from_f64((f64::PI() * x).sin().ln())
            - s.ln()
            - consts.ln_2_sqrt_e_over_pi
            - MyFloat::from_f64((0.5 - x) * ((0.5 - x + consts.gamma_r) / f64::E()).ln())
    } else {
        let s = consts.gamma_dk
            .iter()
            .enumerate()
            .skip(1)
            .fold(consts.gamma_dk[0], |s, t| s + *t.1 / (MyFloat::from_f64(x) + MyFloat::from_f32(t.0 as f32) - MyFloat::one()));

        s.ln()
            + consts.ln_2_sqrt_e_over_pi
            + MyFloat::from_f64((x - 0.5) * ((x - 0.5 + consts.gamma_r) / f64::E()).ln())
    }
}