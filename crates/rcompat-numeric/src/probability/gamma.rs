//! Stable gamma normalization and chi-square tails.
use super::quadrature::integrate_positive_tail;
use crate::numeric::math::ProfileMath;
use crate::numeric::special::{gamma_lr, gamma_ur};

/// log Gamma(a+1/2) - log Gamma(a) - (1/2)log(a), a>=32.
/// Bernoulli-polynomial gamma-ratio expansion, DLMF 5.11.13.
pub fn log_half_gamma_ratio_correction(a: f64) -> f64 {
    let inverse = 1.0 / a;
    let square = inverse * inverse;
    inverse
        * (-0.125
            + square
                * (1.0 / 192.0
                    + square
                        * (-1.0 / 640.0 + square * (17.0 / 14336.0 - square * (31.0 / 18432.0)))))
}

/// log Gamma(a) minus its leading Stirling expression, a>=1/2. Shift small
/// shapes using a positive atanh series, rather than subtracting large logs.
pub fn stirling_correction(mut a: f64) -> f64 {
    let mut correction = 0.0;
    let mut compensation = 0.0;
    while a < 32.0 {
        let t = 1.0 / (2.0 * a + 1.0);
        let square = t * t;
        let mut term = square / 3.0;
        let mut total = term;
        for k in 1..128 {
            let k = f64::from(k);
            term *= square * (2.0 * k + 1.0) / (2.0 * k + 3.0);
            total += term;
            if term <= total * f64::EPSILON {
                break;
            }
        }
        let delta = total - compensation;
        let next = correction + delta;
        compensation = (next - correction) - delta;
        correction = next;
        a += 1.0;
    }
    let inv = 1.0 / a;
    let square = inv * inv;
    correction
        + inv
            * (1.0 / 12.0
                + square
                    * (-1.0 / 360.0
                        + square
                            * (1.0 / 1260.0
                                + square
                                    * (-1.0 / 1680.0
                                        + square * (1.0 / 1188.0 - square * (691.0 / 360360.0))))))
}

/// Stable remainder log(1+x)-x, needed in gamma and beta log densities.
pub fn log1pmx(x: f64) -> f64 {
    if x.abs() > 0.125 {
        return x.ln_1p_p() - x;
    }
    let mut term = -0.5 * x * x;
    let mut sum = term;
    for k in 2..128 {
        let k = f64::from(k);
        term *= -x * k / (k + 1.0);
        sum += term;
        if term.abs() <= sum.abs() * f64::EPSILON {
            break;
        }
    }
    sum
}

/// Upper tail of a chi-square statistic for large degrees of freedom.
pub fn large_shape_chi_square_tail(statistic: f64, df: f64) -> f64 {
    let a = 0.5 * df;
    let x = 0.5 * statistic.max(0.0);
    if x == 0.0 {
        return 1.0;
    }
    let root = a.sqrt();
    let z = (x - a) / root;
    if z < -12.0 {
        return 1.0;
    }
    let correction = stirling_correction(a);
    // Integrate outward from the requested tail. Reflect the left tail so
    // its density peak never falls inside a coarse distant quadrature panel.
    let lower = z < 0.0;
    let integral = integrate_positive_tail(|s| {
        let u = (if lower { z - s } else { z + s }) / root;
        if u <= -1.0 {
            return 0.0;
        }
        0.3989422804014327 * (-correction + a * log1pmx(u) - u.ln_1p_p()).exp_p()
    });
    (if lower { 1.0 - integral } else { integral }).clamp(0.0, 1.0)
}

/// One tail of the chi-square distribution with `df` degrees of freedom at
/// `x > 0`: `P(X <= x)` when `lower_tail`, otherwise `P(X > x)`.
fn chi_square_tail(x: f64, df: f64, lower_tail: bool) -> f64 {
    if df == 2.0 {
        let half = -0.5 * x;
        return if lower_tail {
            -half.exp_m1_p()
        } else {
            half.exp_p()
        };
    }
    if df >= 64.0 {
        let upper = large_shape_chi_square_tail(x, df);
        return if lower_tail { 1.0 - upper } else { upper };
    }
    if lower_tail {
        gamma_lr(0.5 * df, 0.5 * x)
    } else {
        gamma_ur(0.5 * df, 0.5 * x)
    }
}

/// Quantile of the chi-square distribution with `df > 0` degrees of freedom:
/// the `x` with `P(X <= x) = probability` when `lower_tail`, otherwise the `x`
/// with `P(X > x) = probability`, for `0 < probability < 1`.
///
/// Bisects on `log x` until the bracket collapses to adjacent doubles. The
/// tail that is compared is whichever is smaller at the target, so its
/// relative accuracy carries over to the quantile.
pub fn chi_square_quantile(probability: f64, df: f64, lower_tail: bool) -> f64 {
    let (target, compare_lower) = if probability <= 0.5 {
        (probability, lower_tail)
    } else {
        (1.0 - probability, !lower_tail)
    };
    // Lower tail increases with x, upper tail decreases.
    let above = |x: f64| {
        let tail = chi_square_tail(x, df, compare_lower);
        if compare_lower {
            tail >= target
        } else {
            tail <= target
        }
    };
    let centre = df.ln_p();
    let mut low = centre - 1.0;
    let mut high = centre + 1.0;
    let mut step = 1.0;
    while above(low.exp_p()) {
        high = low;
        step *= 2.0;
        low -= step;
    }
    step = 1.0;
    while !above(high.exp_p()) {
        low = high;
        step *= 2.0;
        high += step;
    }
    // Halve the log bracket to a narrow relative width, then finish in x
    // itself so the last steps resolve adjacent doubles of the quantile.
    while high - low > 1e-3 {
        let mid = low + 0.5 * (high - low);
        if above(mid.exp_p()) {
            high = mid;
        } else {
            low = mid;
        }
    }
    let mut low = low.exp_p();
    let mut high = high.exp_p();
    loop {
        let mid = low + 0.5 * (high - low);
        if mid <= low || mid >= high {
            break;
        }
        if above(mid) {
            high = mid;
        } else {
            low = mid;
        }
    }
    low + 0.5 * (high - low)
}
