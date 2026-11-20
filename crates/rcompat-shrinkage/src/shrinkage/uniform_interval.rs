//! Interval integration shared by portable location likelihoods.
use crate::numeric::math::ProfileMath;

/// Simpson integration in the interval's original units.
pub(super) struct NarrowInterval {
    log_density: f64,
    left_weight: f64,
    right_weight: f64,
    total: f64,
}

impl NarrowInterval {
    pub(super) fn from_log_ratios(center: f64, left: f64, right: f64) -> Option<Self> {
        if !center.is_finite() {
            return None;
        }
        let left_weight = left.exp_p();
        let right_weight = right.exp_p();
        let total = left_weight + 4.0 + right_weight;
        Some(Self {
            log_density: center + (total / 6.0).ln_p(),
            left_weight,
            right_weight,
            total,
        })
    }

    pub(super) fn log_density(&self, se: f64) -> f64 {
        self.log_density - se.ln_p()
    }

    fn centered_moments(&self, half: f64) -> (f64, f64) {
        let centered_mean = half * ((self.right_weight - self.left_weight) / self.total);
        let centered_second = half * (half * ((self.right_weight + self.left_weight) / self.total));
        (centered_mean, centered_second)
    }

    pub(super) fn moments(&self, a: f64, b: f64) -> (f64, f64) {
        let center = a.midpoint(b);
        let (centered_mean, centered_second) = self.centered_moments((b - a) / 2.0);
        (
            center + centered_mean,
            center * (center + 2.0 * centered_mean) + centered_second,
        )
    }

    pub(super) fn mean_variance(&self, a: f64, b: f64) -> (f64, f64) {
        let (mean_delta, second_delta) = self.centered_moments((b - a) / 2.0);
        (
            a.midpoint(b) + mean_delta,
            second_delta - mean_delta * mean_delta,
        )
    }
}

pub(super) fn symmetric_log_mass(mut a: f64, mut b: f64, log_cdf: impl Fn(f64) -> f64) -> f64 {
    if a >= 0.0 {
        (a, b) = (-b, -a);
    }
    let upper = log_cdf(b);
    let lower = log_cdf(a);
    if lower == f64::NEG_INFINITY {
        return upper;
    }
    upper + (-(lower - upper).exp_m1_p()).ln_p()
}
