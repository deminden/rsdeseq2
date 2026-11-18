//! Standalone probability calculations shared by fitted and primitive inference.

use crate::numeric::math::ProfileMath;
/// Compensate argument reduction against sqrt(2), including the low part of the
/// irrational constant. Without it, rounding q/sqrt(2) alone loses several ULPs
/// in small tails even when erfc itself is accurate.
pub fn normal_two_sided_tail(q: f64) -> f64 {
    let q = q.abs();
    let x = q / std::f64::consts::SQRT_2;
    let tail = libm::erfc(x);
    if !x.is_finite() || tail == 0.0 {
        return tail;
    }
    const SQRT_2_LOW: f64 = -9.667293313452913e-17;
    let remainder = (-x).mul_add(std::f64::consts::SQRT_2, q) - x * SQRT_2_LOW;
    let correction = std::f64::consts::FRAC_2_SQRT_PI
        * (-x * x).exp_p()
        * (remainder / std::f64::consts::SQRT_2);
    (tail - correction).clamp(0.0, 1.0)
}

/// `P(Z > q)` for a standard normal `Z`.
pub fn normal_upper_tail(q: f64) -> f64 {
    let tail = 0.5 * normal_two_sided_tail(q);
    if q < 0.0 { 1.0 - tail } else { tail }
}

/// Shared inverse upper-normal tail, refined against the compensated tail.
pub fn normal_upper_quantile(probability: f64) -> f64 {
    if !(0.0..=1.0).contains(&probability) {
        return f64::NAN;
    }
    if probability == 0.0 {
        return f64::INFINITY;
    }
    if probability == 1.0 {
        return f64::NEG_INFINITY;
    }
    let tail = probability.min(1.0 - probability);
    let mut q = -crate::numeric::special::standard_normal_inverse_cdf(tail);
    for _ in 0..2 {
        let density = 0.3989422804014327 * (-0.5 * q * q).exp_p();
        let correction = (normal_upper_tail(q) - tail) / density;
        if !correction.is_finite() {
            break;
        }
        q += correction;
    }
    if probability > 0.5 { -q } else { q }
}
