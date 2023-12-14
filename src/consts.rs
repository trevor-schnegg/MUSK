use crate::my_float::MyFloat;

pub struct Consts {
    pub(crate) gamma_r: f64,
    pub(crate) gamma_dk: Vec<MyFloat>,
    pub(crate) ln_pi: MyFloat,
    pub(crate) ln_2_sqrt_e_over_pi: MyFloat,
}

impl Consts {
    pub fn new() -> Self {
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
            ln_2_sqrt_e_over_pi: MyFloat::from_f64(0.6207822376352452223455184457816472122518527279025978),
        }
    }
}
