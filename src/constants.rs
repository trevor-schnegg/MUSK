use rug::float::{Constant, Special};
use rug::Float;

#[derive(Debug)]
pub struct Constants {
    pub(crate) zero: Float,
    pub(crate) one: Float,
    pub(crate) two: Float,
    pub(crate) eps: Float,
    pub(crate) pi: Float,
    pub(crate) e: Float,
    pub(crate) ln_pi: Float,
    pub(crate) ln_2_sqrt_e_over_pi: Float,
    pub(crate) gamma_r: Float,
    pub(crate) gamma_dk: Vec<Float>,
}

impl Constants {
    pub fn new() -> Self {
        Constants {
            zero: Float::with_val(256, Special::Zero),
            one: Float::with_val(256, 1.0),
            two: Float::with_val(256, 2.0),
            eps: Float::with_val(256, f64::EPSILON),
            pi: Float::with_val(256, Constant::Pi),
            e: Float::with_val(256, std::f64::consts::E),
            ln_pi: Float::with_val(256, 1.1447298858494001741434273513530587116472948129153),
            ln_2_sqrt_e_over_pi: Float::with_val(
                256,
                0.6207822376352452223455184457816472122518527279025978,
            ),
            gamma_r: Float::with_val(256, 10.900511),
            gamma_dk: Vec::from([
                Float::with_val(256, 2.48574089138753565546e-5),
                Float::with_val(256, 1.05142378581721974210),
                Float::with_val(256, -3.45687097222016235469),
                Float::with_val(256, 4.51227709466894823700),
                Float::with_val(256, -2.98285225323576655721),
                Float::with_val(256, 1.05639711577126713077),
                Float::with_val(256, -1.95428773191645869583e-1),
                Float::with_val(256, 1.70970543404441224307e-2),
                Float::with_val(256, -5.71926117404305781283e-4),
                Float::with_val(256, 4.63399473359905636708e-6),
                Float::with_val(256, -2.71994908488607703910e-9),
            ]),
        }
    }
}
