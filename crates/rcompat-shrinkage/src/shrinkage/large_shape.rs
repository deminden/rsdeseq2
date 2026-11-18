//! Scaled gamma/beta integrals and densities without subtracting large log gammas.
//! Integral definitions: DLMF 8.2 and 8.17; normalization: DLMF 5.11.
//! These are coordinate changes of the distributions, not normal approximations.
use crate::numeric::math::ProfileMath;
use crate::probability::{log1pmx, quadrature::integrate_positive_tail, stirling_correction};
const LOG_2_PI: f64 = 1.8378770664093453;

pub(super) fn complement(log: f64) -> f64 {
    if log < -std::f64::consts::LN_2 {
        (-log.exp_p()).ln_1p_p()
    } else {
        (-log.exp_m1_p()).ln_p()
    }
}
fn tail(log_density: impl Fn(f64) -> f64, at: f64, radius: f64, lower: bool) -> f64 {
    let scale = log_density(at);
    if scale == f64::NEG_INFINITY {
        return scale;
    }
    let integral = integrate_positive_tail(|s| {
        let t = if lower {
            at - s * radius
        } else {
            at + s * radius
        };
        (log_density(t) - scale).exp_p() * radius
    });
    scale + integral.ln_p()
}
pub(super) fn gamma_log_probability(a: f64, x: f64, lower: bool) -> f64 {
    let root = a.sqrt();
    let z = (x - a) / root;
    let correction = stirling_correction(a);
    let density = |z: f64| {
        let u = z / root;
        if u <= -1.0 {
            f64::NEG_INFINITY
        } else {
            -0.5 * LOG_2_PI - correction + a * log1pmx(u) - u.ln_1p_p()
        }
    };
    let direct_lower = x < a;
    let slope = ((x - a + 1.0) / x).abs() * root;
    let log = tail(density, z, 1.0 / (1.0 + slope), direct_lower);
    if lower == direct_lower {
        log
    } else {
        complement(log)
    }
}
pub(super) fn gamma_log_density(x: f64, a: f64, rate: f64) -> f64 {
    let y = x * rate;
    if y.is_infinite() {
        return f64::NEG_INFINITY;
    }
    let u = (y - a) / a;
    if u <= -1.0 {
        return a * rate.ln_p() + (a - 1.0) * x.ln_p()
            - rate * x
            - crate::numeric::special::ln_gamma(a);
    }
    rate.ln_p() - 0.5 * (LOG_2_PI + a.ln_p()) - stirling_correction(a) + a * log1pmx(u)
        - u.ln_1p_p()
}
pub(super) fn beta_large(a: f64, b: f64) -> bool {
    a.min(b) >= 64.0 && a.max(b) >= 1e5
}
fn expm1mx(x: f64) -> f64 {
    let mut term = 0.5 * x * x;
    let mut total = term;
    for k in 2..128 {
        term *= x / (f64::from(k) + 1.0);
        total += term;
        if term.abs() <= f64::EPSILON * total.abs() {
            break;
        }
    }
    total
}
// p<=1/2 prevents cancellation between almost equal quadratic terms.
fn loss(eta: f64, p: f64) -> f64 {
    if eta.abs() < 0.1 {
        log1pmx(p * eta.exp_m1_p()) + p * expm1mx(eta)
    } else if eta > 700.0 {
        eta + (p + (1.0 - p) * (-eta).exp_p()).ln_p() - p * eta
    } else {
        (p * eta.exp_m1_p()).ln_1p_p() - p * eta
    }
}
fn correction(a: f64, b: f64) -> f64 {
    stirling_correction(a + b) - stirling_correction(a) - stirling_correction(b)
}
pub(super) fn beta_logit_density(mut a: f64, mut b: f64, mut eta: f64) -> f64 {
    if a > b {
        std::mem::swap(&mut a, &mut b);
        eta = -eta;
    }
    let total = a + b;
    let p = a / total;
    0.5 * (a.ln_p() + b.ln_p() - total.ln_p() - LOG_2_PI) + correction(a, b) - total * loss(eta, p)
}
fn coordinate(x: f64, a: f64, b: f64) -> f64 {
    let p = a / (a + b);
    let q = b / (a + b);
    ((x - p) / p).ln_1p_p() - (((1.0 - x) - q) / q).ln_1p_p()
}
pub(super) fn beta_log_probability(mut a: f64, mut b: f64, mut x: f64, mut lower: bool) -> f64 {
    if a > b {
        std::mem::swap(&mut a, &mut b);
        x = 1.0 - x;
        lower = !lower;
    }
    let total = a + b;
    let p = a / total;
    let sigma = (1.0 / a + 1.0 / b).sqrt();
    let eta = coordinate(x, a, b);
    let z = eta / sigma;
    let correction = correction(a, b);
    let density = |z: f64| -0.5 * LOG_2_PI + correction - total * loss(sigma * z, p);
    let slope = (total * x - a).abs() * sigma;
    let direct_lower = x < p;
    let log = tail(density, z, 1.0 / (1.0 + slope), direct_lower);
    if lower == direct_lower {
        log
    } else {
        complement(log)
    }
}
pub(super) fn beta_log_density(mut x: f64, mut a: f64, mut b: f64) -> f64 {
    if a > b {
        std::mem::swap(&mut a, &mut b);
        x = 1.0 - x;
    }
    beta_logit_density(a, b, coordinate(x, a, b)) - x.ln_p() - (-x).ln_1p_p()
}
// A small shift of a large shape, from the gamma ratio in Stirling's formula.
pub(super) fn log_beta(a: f64, b: f64) -> f64 {
    let small = a.min(b);
    let large = a.max(b);
    let r = small / large;
    let ratio = small * large.ln_p()
        + (large + small - 0.5) * log1pmx(r)
        + (small - 0.5) * r
        + stirling_correction(large + small)
        - stirling_correction(large);
    crate::numeric::special::ln_gamma(small) - ratio
}
