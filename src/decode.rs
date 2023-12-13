pub fn integer_decode_f32(f: f32) -> (f32, i32) {
    let bits = f.to_bits();
    let mut exponent = ((bits >> 23) & 0xff) as i32;
    let mantissa = f32::from_bits((bits & 0x807fffff) | 0x3f800000);
    // Exponent bias
    exponent -= 127;
    (mantissa, exponent)
}

pub fn integer_decode_f64(f: f64) -> (f32, i32) {
    let bits = f.to_bits();
    let mut exponent = ((bits >> 52) & 0x7ff) as i32;
    let mantissa = f64::from_bits((bits & 0x800fffffffffffff) | 0x3ff0000000000000) as f32;
    // Exponent bias + mantissa shift
    exponent -= 1023;
    (mantissa, exponent)
}