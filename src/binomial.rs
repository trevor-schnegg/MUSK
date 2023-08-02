use crate::constants::Constants;
use rug::ops::CompleteRound;
use rug::Assign;
use rug::Float;
use statrs::StatsError;

#[derive(Debug)]
pub struct Binomial {
    p: Float,
    n: u64,
    constants: Constants,
}

impl Binomial {
    pub fn new(p: Float, n: u64) -> Result<Binomial, StatsError> {
        let zero = Float::with_val(256, 0.0);
        let one = Float::with_val(256, 1.0);
        if p.is_nan() || p < zero || p > one {
            Err(StatsError::BadParams)
        } else {
            let constants = Constants::new();
            Ok(Binomial { p, n, constants })
        }
    }

    pub fn p(&self) -> &Float {
        &self.p
    }

    pub fn n(&self) -> u64 {
        self.n
    }

    pub fn sf(&self, x: u64) -> Float {
        if x >= self.n {
            self.constants.zero.clone()
        } else {
            let k = x;
            self.beta_reg(
                Float::with_val(256, k) + &self.constants.one,
                Float::with_val(256, self.n - k),
                self.p.clone(),
            )
        }
    }

    fn beta_reg(&self, a: Float, b: Float, x: Float) -> Float {
        self.checked_beta_reg(a, b, x).unwrap()
    }

    fn checked_beta_reg(&self, a: Float, b: Float, x: Float) -> Result<Float, StatsError> {
        if a <= self.constants.zero {
            Err(StatsError::ArgMustBePositive("a"))
        } else if b <= self.constants.zero {
            Err(StatsError::ArgMustBePositive("b"))
        } else if !(self.constants.zero.clone()..=self.constants.one.clone()).contains(&x) {
            Err(StatsError::ArgIntervalIncl("x", 0.0, 1.0))
        } else {
            let bt = if x.is_zero() || x == self.constants.one {
                self.constants.zero.clone()
            } else {
                (self.ln_gamma((&a + &b).complete(256))
                    - self.ln_gamma(a.clone())
                    - self.ln_gamma(b.clone())
                    + &a * x.clone().ln()
                    + &b * (&self.constants.one - &x).complete(256).ln())
                .exp()
            };
            let symm_transform = x
                >= (&a + &self.constants.one).complete(256)
                    / ((&a + &b).complete(256) + &self.constants.two);
            let eps = self.constants.eps.clone();
            // let fpmin = constants::FPMIN;

            let mut a = a;
            let mut b = b;
            let mut x = x;
            if symm_transform {
                let swap = a;
                x = &self.constants.one - x;
                a = b;
                b = swap;
            }

            let qab = (&a + &b).complete(256);
            let qap = (&a + &self.constants.one).complete(256);
            let qam = (&a - &self.constants.one).complete(256);
            let mut c = self.constants.one.clone();
            let mut d =
                Float::with_val(256, &self.constants.one - (&qab * &x).complete(256) / &qap);
            // if d.abs() < fpmin {
            //     d = fpmin;
            // }
            d = (&self.constants.one / &d).complete(256);
            let mut h = d.clone();

            for m in 1..141 {
                let m = Float::with_val(256, m);
                let m2 = (&m * &self.constants.two).complete(256);
                let mut aa = Float::new(256);
                aa.assign(
                    &m * (&b - &m).complete(256) * &x
                        / ((&qam + &m2).complete(256) * (&a + &m2).complete(256)),
                );
                d = (&self.constants.one + &aa * &d).complete(256);

                // if d.abs() < fpmin {
                //     d = fpmin;
                // }

                c.assign(&self.constants.one + (&aa / &c).complete(256));
                // if c.abs() < fpmin {
                //     c = fpmin;
                // }

                d = (&self.constants.one / &d).complete(256);
                h.assign((&h * &d).complete(256) * &c);
                aa.assign(
                    -(&a + &m).complete(256) * (&qab + &m).complete(256) * &x
                        / ((&a + &m2).complete(256) * (&qap + &m2).complete(256)),
                );
                d = (&self.constants.one + &aa * &d).complete(256);

                // if d.abs() < fpmin {
                //     d = fpmin;
                // }

                c.assign(&self.constants.one + (&aa / &c).complete(256));

                // if c.abs() < fpmin {
                //     c = fpmin;
                // }

                d = (&self.constants.one / &d).complete(256);
                let del = (&d * &c).complete(256);
                h *= &del;

                if (&del - &self.constants.one).complete(256).abs() <= eps {
                    return if symm_transform {
                        Ok(&self.constants.one - bt * h / a)
                    } else {
                        Ok(bt * h / a)
                    };
                }
            }

            if symm_transform {
                Ok(&self.constants.one - bt * h / a)
            } else {
                Ok(bt * h / a)
            }
        }
    }

    fn ln_gamma(&self, x: Float) -> Float {
        if x < Float::with_val(256, 0.5) {
            let s = self
                .constants
                .gamma_dk
                .iter()
                .enumerate()
                .skip(1)
                .fold(self.constants.gamma_dk[0].clone(), |s, t| {
                    s + t.1 / (Float::with_val(256, t.0) - &x)
                });
            &self.constants.ln_pi
                - (&self.constants.pi * &x).complete(256).sin().ln()
                - s.ln()
                - &self.constants.ln_2_sqrt_e_over_pi
                - (Float::with_val(256, 0.5) - &x)
                    * ((Float::with_val(256, 0.5) - x + &self.constants.gamma_r)
                        / &self.constants.e)
                        .ln()
        } else {
            let s = self
                .constants
                .gamma_dk
                .iter()
                .enumerate()
                .skip(1)
                .fold(self.constants.gamma_dk[0].clone(), |s, t| {
                    s + t.1 / (&x + Float::with_val(256, t.0) - &self.constants.one)
                });

            s.ln()
                + &self.constants.ln_2_sqrt_e_over_pi
                + (&x - Float::with_val(256, 0.5))
                    * ((x - Float::with_val(256, 0.5) + &self.constants.gamma_r)
                        / &self.constants.e)
                        .ln()
        }
    }
}
