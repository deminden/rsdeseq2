//! Standard truncated-normal moments with finite-arithmetic boundary bounds.
use crate::{math::elementary::language_exp, math::normal_log_cdf};

pub(super) fn log_phi(x: f64) -> f64 {
    // Match the pinned platform's original dnorm rounding. A one-ulp density
    // change is amplified when far-tail posterior variances subtract moments.
    #[cfg(all(target_os = "linux", target_arch = "x86_64"))]
    {
        -(0.918_938_533_204_672_8 + 0.5 * x * x)
    }
    #[cfg(not(all(target_os = "linux", target_arch = "x86_64")))]
    {
        -(0.5 * x).mul_add(x, 0.918_938_533_204_672_8)
    }
}
// Finite-arithmetic log difference, including coincident/zero endpoints.
fn log_difference(larger: f64, smaller: f64) -> f64 {
    if larger.is_nan() || smaller.is_nan() {
        return f64::NAN;
    }
    if larger <= smaller {
        return f64::NEG_INFINITY;
    }
    if larger == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    larger + (1.0 - language_exp(smaller - larger)).ln()
}
pub(super) fn truncated_mean(mut a: f64, mut b: f64) -> f64 {
    let flipped = (a > 0.0 && b > 0.0) || b > a.abs();
    if flipped {
        (a, b) = (-b, -a);
    }
    let denominator = log_difference(normal_log_cdf(b), normal_log_cdf(a));
    let mut mean = -language_exp(log_difference(log_phi(b), log_phi(a)) - denominator);
    if denominator.is_infinite() {
        mean = (a + b) / 2.0;
    }
    if b < 0.0 {
        let bound = (b + 1.0 / b).max((a + b) / 2.0);
        if mean < bound || mean > b {
            mean = bound;
        }
    }
    if flipped { -mean } else { mean }
}
pub(super) fn truncated_second(mut a: f64, mut b: f64) -> f64 {
    if a > 0.0 && b > 0.0 {
        (a, b) = (-b, -a);
    }
    let denominator = log_difference(normal_log_cdf(b), normal_log_cdf(a));
    let mut moment = 1.0;
    if a.is_finite() {
        moment = 1.0 + a * language_exp(log_phi(a) - denominator);
    }
    if b.is_finite() {
        moment -= b * language_exp(log_phi(b) - denominator);
    }
    if denominator.is_infinite() {
        moment = (a + b).powi(2) / 4.0;
    }
    if b < 0.0 {
        let endpoint = b * b;
        let bound = (endpoint + 2.0 * (1.0 + 1.0 / endpoint)).min((a * a + a * b + b * b) / 3.0);
        if moment < endpoint || moment > bound {
            moment = bound;
        }
    }
    moment
}

pub(super) fn scaled_mean(a: f64, b: f64, mean: f64, sd: f64) -> f64 {
    if sd == 0.0 {
        mean.clamp(a, b)
    } else {
        mean + sd * truncated_mean((a - mean) / sd, (b - mean) / sd)
    }
}
pub(super) fn scaled_sd(a: f64, b: f64, mean: f64, sd: f64) -> f64 {
    let alpha = (a - mean) / sd;
    let beta = (b - mean) / sd;
    let m = truncated_mean(alpha, beta);
    let mut v = truncated_second(alpha, beta) - m * m;
    let nearest = alpha.abs().min(beta.abs());
    if nearest.is_finite() && nearest > 30.0 {
        v = (1.0 / (nearest * nearest)).min((beta - alpha).powi(2) / 12.0);
    }
    if v.is_nan() {
        v = 0.0;
    }
    (sd * sd * v).sqrt()
}
