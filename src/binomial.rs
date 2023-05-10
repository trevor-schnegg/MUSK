use statrs::distribution::DiscreteCDF;
use statrs::{StatsError};
use f128::f128;
use statrs::statistics::{Max, Min};
use num_traits::float::Float;
use num_traits::Zero;

// const F128_PREC: f128 = f128::EPSILON;

pub fn ulps_eq(a: f128, b: f128) -> bool {
    if (a - b).abs() < f128::MIN_POSITIVE_NORMAL {
        true
    } else {
        false
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Binomial {
    p: f128,
    n: u64,
    ln_pi: f128,
    ln_2_sqrt_e_over_pi: f128,
    gamma_r: f128,
    gamma_dk: Vec<f128>,
}

impl Binomial {
    pub fn new(p: f128, n: u64) -> Result<Binomial, StatsError> {
        if p.is_nan() || p < f128::from(0.0) || p > f128::from(1.0) {
            Err(StatsError::BadParams)
        } else {
            Ok(Binomial { p, n,
            ln_pi: f128::from(1.1447298858494001741434273513530587116472948129153),
            ln_2_sqrt_e_over_pi: f128::from(0.6207822376352452223455184457816472122518527279025978),
            gamma_r: f128::from(10.900511),
            gamma_dk: Vec::from([
                f128::from(2.48574089138753565546e-5),
                f128::from(1.05142378581721974210),
                f128::from(-3.45687097222016235469),
                f128::from(4.51227709466894823700),
                f128::from(-2.98285225323576655721),
                f128::from(1.05639711577126713077),
                f128::from(-1.95428773191645869583e-1),
                f128::from(1.70970543404441224307e-2),
                f128::from(-5.71926117404305781283e-4),
                f128::from(4.63399473359905636708e-6),
                f128::from(-2.71994908488607703910e-9)])
            })
        }
    }

    pub fn p(&self) -> f128 {
        self.p
    }

    pub fn n(&self) -> u64 {
        self.n
    }

    fn beta_reg(&self, a: f128, b: f128, x: f128) -> f128 {
        self.checked_beta_reg(a, b, x).unwrap()
    }

    fn checked_beta_reg(&self, a: f128, b: f128, x: f128) -> Result<f128, StatsError> {
        if a <= f128::ZERO {
            Err(StatsError::ArgMustBePositive("a"))
        } else if b <= f128::ZERO {
            Err(StatsError::ArgMustBePositive("b"))
        } else if !(f128::ZERO..=f128::ONE).contains(&x) {
            Err(StatsError::ArgIntervalIncl("x", 0.0, 1.0))
        } else {
            let bt = if x.is_zero() || ulps_eq(x, f128::ONE) {
                f128::ZERO
            } else {
                (self.ln_gamma(a + b) - self.ln_gamma(a) - self.ln_gamma(b)
                    + a * x.ln()
                    + b * (f128::ONE - x).ln())
                    .exp()
            };
            let symm_transform = x >= (a + f128::ONE) / (a + b + f128::TWO);
            let eps = f128::EPSILON;
            let fpmin = f128::MIN_POSITIVE_NORMAL;

            let mut a = a;
            let mut b = b;
            let mut x = x;
            if symm_transform {
                let swap = a;
                x = f128::ONE - x;
                a = b;
                b = swap;
            }

            let qab = a + b;
            let qap = a + f128::ONE;
            let qam = a - f128::ONE;
            let mut c = f128::ONE;
            let mut d = f128::ONE - qab * x / qap;
            if d.abs() < fpmin {
                d = fpmin;
            }
            d = f128::ONE / d;
            let mut h = d;

            for m in 1..141 {
                let m = f128::from(m);
                let m2 = m * f128::TWO;
                let mut aa = m * (b - m) * x / ((qam + m2) * (a + m2));
                d = f128::ONE + aa * d;

                if d.abs() < fpmin {
                    d = fpmin;
                }

                c = f128::ONE + aa / c;
                if c.abs() < fpmin {
                    c = fpmin;
                }

                d = f128::ONE / d;
                h = h * d * c;
                aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
                d = f128::ONE + aa * d;

                if d.abs() < fpmin {
                    d = fpmin;
                }

                c = f128::ONE + aa / c;

                if c.abs() < fpmin {
                    c = fpmin;
                }

                d = f128::ONE / d;
                let del = d * c;
                h *= del;

                if (del - f128::ONE).abs() <= eps {
                    return if symm_transform {
                        Ok(f128::ONE - bt * h / a)
                    } else {
                        Ok(bt * h / a)
                    };
                }
            }

            if symm_transform {
                Ok(f128::ONE - bt * h / a)
            } else {
                Ok(bt * h / a)
            }
        }
    }

    fn ln_gamma(&self, x: f128) -> f128 {
        if x < f128::from(0.5) {
            let s = self.gamma_dk
                .iter()
                .enumerate()
                .skip(1)
                .fold(self.gamma_dk[0], |s, t| s + t.1 / (f128::from(t.0) - x));
            self.ln_pi
                - (f128::PI * x).sin().ln()
                - s.ln()
                - self.ln_2_sqrt_e_over_pi
                - (f128::from(0.5) - x) * ((f128::from(0.5) - x + self.gamma_r) / f128::E).ln()
        } else {
            let s = self.gamma_dk
                .iter()
                .enumerate()
                .skip(1)
                .fold(self.gamma_dk[0], |s, t| s + t.1 / (x + f128::from(t.0) - f128::ONE));

            s.ln()
                + self.ln_2_sqrt_e_over_pi
                + (x - f128::from(0.5)) * ((x - f128::from(0.5) + self.gamma_r) / f128::E).ln()
        }
    }
}

impl Min<u64> for Binomial {
    fn min(&self) -> u64 {
        0
    }
}

impl Max<u64> for Binomial {
    fn max(&self) -> u64 {
        self.n
    }
}

impl DiscreteCDF<u64, f128> for Binomial {
    fn cdf(&self, x: u64) -> f128 {
        if x >= self.n {
            f128::from(1.0)
        } else {
            let k = x;
            let one = f128::from(1.0);
            self.beta_reg(f128::from(self.n - k), f128::from(k) + one, one - self.p)
        }
    }

    fn sf(&self, x: u64) -> f128 {
        if x >= self.n {
            f128::from(0.0)
        } else {
            let k = x;
            let one = f128::from(1.0);
            self.beta_reg(f128::from(k) + one, f128::from(self.n - k), self.p)
        }
    }
}
