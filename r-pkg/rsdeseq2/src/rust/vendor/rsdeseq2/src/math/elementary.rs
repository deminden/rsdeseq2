//! Elementary expression semantics shared by formula and posterior calculations.

// The expression-level exponential rounds its first-order expansion near zero.
// Internal distribution kernels retain their ordinary mathematical exp calls.
pub(crate) fn language_exp(x: f64) -> f64 {
    if x.abs() <= f64::EPSILON.sqrt() {
        1.0 + x
    } else {
        x.exp()
    }
}

// Vector expression evaluation handles exact zero before the small-argument
// branch. Keep these semantics out of distribution and optimizer kernels.
pub(crate) fn language_log1p(x: f64) -> f64 {
    if x == 0.0 {
        0.0
    } else if x.abs() <= f64::EPSILON {
        x
    } else {
        x.ln_1p()
    }
}
pub(crate) fn language_expm1(x: f64) -> f64 {
    if x == 0.0 {
        0.0
    } else if x.abs() <= f64::EPSILON {
        x
    } else {
        x.exp_m1()
    }
}
