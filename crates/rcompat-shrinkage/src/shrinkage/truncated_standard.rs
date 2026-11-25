//! Standard truncated-normal moments with finite-arithmetic boundary bounds.
use super::{
    log_f_quadrature::weight_moments,
    uniform_interval::{NarrowInterval, symmetric_log_mass},
};
use crate::numeric::math::ProfileMath;
use crate::{math::elementary::language_exp, math::normal_log_cdf};

pub(super) fn log_phi(x: f64) -> f64 {
    // A one-ulp density change is amplified when far-tail posterior
    // variances subtract moments, so the rounding is fixed here.
    -(0.5 * x).mul_add(x, 0.918_938_533_204_672_8)
}

pub(super) fn log_interval_mass(a: f64, b: f64) -> f64 {
    symmetric_log_mass(a, b, normal_log_cdf)
}

pub(super) fn far_tail_cdf(a: f64, t: f64, b: f64, estimate: f64, se: f64) -> Option<f64> {
    let estimate_left = estimate < a;
    let (distance, near_width, far_width) = if estimate_left {
        ((a - estimate) / se, t - a, b - t)
    } else {
        ((estimate - b) / se, b - t, t - a)
    };
    if !distance.is_finite() || distance < 65_536.0 {
        return None;
    }
    // Beyond this distance the leading Mills ratio 1/d gives the relative
    // tail probabilities to f64 precision. Difference the squared distances
    // analytically using the original interval widths.
    let log_tail_ratio = |d: f64, step: f64| -(d + 0.5 * step) * step - (step / d).ln_1p_p();
    let near_step = near_width / se;
    let near_tail = log_tail_ratio(distance, near_step);
    let full_mass = -log_tail_ratio(distance, (b - a) / se).exp_m1_p();
    if full_mass == 0.0 {
        return None;
    }
    let mass = if estimate_left {
        -near_tail.exp_m1_p()
    } else if near_tail == f64::NEG_INFINITY {
        0.0
    } else {
        near_tail.exp_p() * -log_tail_ratio(distance + near_step, far_width / se).exp_m1_p()
    };
    Some((mass / full_mass).clamp(0.0, 1.0))
}

pub(super) fn narrow_interval(a: f64, b: f64, estimate: f64, se: f64) -> Option<NarrowInterval> {
    let half = ((b - a) / se) / 2.0;
    let center = (a.midpoint(b) - estimate) / se;
    let slope = center * half;
    let quadratic = 0.5 * half * half;
    // Bound the log-density variation across the interval. Form the
    // endpoint differences algebraically to preserve tiny widths.
    if slope.abs() + quadratic > 1e-3 {
        return None;
    }
    NarrowInterval::from_log_ratios(log_phi(center), slope - quadratic, -slope - quadratic)
}

fn far_tail_mean_variance(a: f64, b: f64, estimate: f64, se: f64) -> Option<(f64, f64)> {
    let (anchor, direction, distance) = if estimate < a {
        (a, 1.0, (a - estimate) / se)
    } else {
        (b, -1.0, (estimate - b) / se)
    };
    if !distance.is_finite() || distance < 32.0 {
        return None;
    }
    // In units of se/distance from the nearer endpoint the conditional
    // density is proportional to exp(-x - x²/(2*distance²)). The moments
    // beyond 48 contribute less than 4e-18 in these units.
    let end = (distance * ((b - a) / se)).min(48.0);
    let density = |x: f64| (-x - 0.5 * (x / distance).powi(2)).exp_p();
    let mut moments = [0.0; 3];
    let mut left: f64 = 0.0;
    while left < end {
        let right = (left + 4.0).min(end);
        for (moment, value) in moments.iter_mut().zip(weight_moments(density, left, right)) {
            *moment += value;
        }
        left = right;
    }
    let [mass, first, second] = moments;
    if mass <= 0.0 {
        return None;
    }
    let mean = first / mass;
    let variance = second / mass - mean * mean;
    let scale = se / distance;
    Some((
        anchor + direction * scale * mean,
        scale * (scale * variance),
    ))
}

pub(super) fn portable_mean_variance(a: f64, b: f64, estimate: f64, se: f64) -> (f64, f64) {
    if a == b {
        return (a, 0.0);
    }
    if let Some(interval) = narrow_interval(a, b, estimate, se) {
        return interval.mean_variance(a, b);
    }
    if let Some(moments) = far_tail_mean_variance(a, b, estimate, se) {
        return moments;
    }
    let alpha = (a - estimate) / se;
    let beta = (b - estimate) / se;
    let mean = truncated_mean(alpha, beta);
    let variance = truncated_second(alpha, beta) - mean * mean;
    (estimate + se * mean, se * (se * variance))
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
    larger + (1.0 - language_exp(smaller - larger)).ln_p()
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
