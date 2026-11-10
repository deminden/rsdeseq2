//! Positive Gamma normalization from Fullerton's public-domain SLATEC series.
//! Tables: https://netlib.org/slatec/fnlib/dgamma.f . The SLATEC guide states
//! public-domain distribution: https://netlib.org/slatec/guide (section 4).
//! This standalone Rust implementation is used only by posterior components;
//! it neither calls R nor changes the existing GLM probability kernels.
#![allow(clippy::excessive_precision)]
const COEFFICIENTS: [f64; 22] = [
    0.00857119559098933,
    0.004415381324841007,
    0.05685043681599363,
    -0.00421983539641856,
    0.0013268081812124603,
    -0.00018930245297988805,
    3.606925327441245e-05,
    -6.056761904460864e-06,
    1.0558295463022833e-06,
    -1.811967365542384e-07,
    3.117724964715322e-08,
    -5.354219639019687e-09,
    9.193275519859589e-10,
    -1.5779412802883398e-10,
    2.7079806229349544e-11,
    -4.64681865382573e-12,
    7.97335019200742e-13,
    -1.368078209830916e-13,
    2.3473194865638007e-14,
    -4.027432614949067e-15,
    6.910051747372101e-16,
    -1.185584500221993e-16,
];
fn chebyshev(x: f64) -> f64 {
    let mut current = 0.0;
    let mut previous = 0.0;
    let mut before = 0.0;
    for &c in COEFFICIENTS.iter().rev() {
        before = previous;
        previous = current;
        current = (2.0 * x).mul_add(previous, -before) + c;
    }
    (current - before) / 2.0
}
pub(super) fn gamma_positive(x: f64) -> f64 {
    if x.is_nan() || x <= 0.0 {
        return f64::NAN;
    }
    if x <= 10.0 {
        let integer = x.trunc() as usize;
        let fraction = x - integer as f64;
        let mut value = chebyshev(2.0 * fraction - 1.0) + 0.9375;
        if integer == 0 {
            return value / x;
        }
        for i in 1..integer {
            value *= fraction + i as f64;
        }
        return value;
    }
    if x <= 50.0 && x.fract() == 0.0 {
        return (2..x as usize).fold(1.0, |v, i| v * i as f64);
    }
    let inverse = 1.0 / x;
    let square = inverse * inverse;
    let correction = inverse
        * (1.0 / 12.0
            + square
                * (-1.0 / 360.0
                    + square
                        * (1.0 / 1260.0
                            + square
                                * (-1.0 / 1680.0
                                    + square * (1.0 / 1188.0 - square * (691.0 / 360360.0))))));
    ((x - 0.5).mul_add(x.ln(), -x) + 0.9189385332046728 + correction).exp()
}
