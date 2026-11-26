//! Standalone probability calculations shared by fitted and primitive inference.

use crate::numeric::math::ProfileMath;
pub use crate::numeric::normal_tail::{normal_two_sided_tail, normal_upper_tail};

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
