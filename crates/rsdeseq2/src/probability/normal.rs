//! Standalone probability calculations shared by fitted and primitive inference.

// Compensate argument reduction against sqrt(2), including the low part of the
// irrational constant. Without it, rounding q/sqrt(2) alone loses several ULPs
// in small tails even when erfc itself is accurate.
pub(crate) fn normal_two_sided_tail(q: f64) -> f64 {
    let q = q.abs();
    let x = q / std::f64::consts::SQRT_2;
    let tail = libm::erfc(x);
    if !x.is_finite() || tail == 0.0 {
        return tail;
    }
    const SQRT_2_LOW: f64 = -9.667293313452913e-17;
    let remainder = (-x).mul_add(std::f64::consts::SQRT_2, q) - x * SQRT_2_LOW;
    let correction =
        std::f64::consts::FRAC_2_SQRT_PI * (-x * x).exp() * (remainder / std::f64::consts::SQRT_2);
    (tail - correction).clamp(0.0, 1.0)
}

pub(crate) fn normal_upper_tail(q: f64) -> f64 {
    let tail = 0.5 * normal_two_sided_tail(q);
    if q < 0.0 { 1.0 - tail } else { tail }
}

pub(crate) fn normal_log_upper_tail(q: f64) -> f64 {
    if q < 0.0 {
        return (-normal_upper_tail(-q)).ln_1p();
    }
    if q <= 8.0 || q.is_nan() {
        return normal_upper_tail(q).ln();
    }
    // Laplace's continued fraction for the normal Mills ratio. In this
    // range it converges rapidly without first forming an underflowed tail.
    // See NIST DLMF 7.9.1, with z = q / sqrt(2).
    let mut fraction = 0.0;
    for numerator in (1..=64).rev() {
        fraction = numerator as f64 / (q + fraction);
    }
    let mills_ratio = 1.0 / (q + fraction);
    -0.5 * q * q - 0.918_938_533_204_672_8 + mills_ratio.ln()
}
