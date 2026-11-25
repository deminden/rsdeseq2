//! Fixed Cook's F quantile from a stable beta-logistic density.
use super::gamma::{log1pmx, stirling_correction};
use super::quadrature::integrate_positive_tail;
use crate::numeric::math::ProfileMath;

fn expm1mx(x: f64) -> f64 {
    let mut term = 0.5 * x * x;
    let mut sum = term;
    for k in 2..128 {
        term *= x / (f64::from(k) + 1.0);
        sum += term;
        if term.abs() <= sum.abs() * f64::EPSILON {
            break;
        }
    }
    sum
}

fn logistic_loss(y: f64, p: f64) -> f64 {
    if y.abs() < 0.1 {
        log1pmx(p * y.exp_m1_p()) + p * expm1mx(y)
    } else {
        (p * y.exp_m1_p()).ln_1p_p() - p * y
    }
}

/// The 0.99 quantile of the F distribution with `df1` and `df2` degrees of
/// freedom, as DESeq2 uses for its Cook's distance cutoff.
pub fn f_99_quantile(df1: f64, df2: f64) -> f64 {
    // Closed exponential/beta cases avoid an unnecessary inversion.
    if df1 == 2.0 {
        let value = (0.5 * df2) * (-2.0 * (1.0_f64 - 0.99).ln_p() / df2).exp_m1_p();
        return from_inverse_beta(value, df1, df2);
    }
    if df2 == 2.0 {
        let value = 1.0 / ((0.5 * df1) * (-2.0 * 0.99_f64.ln_p() / df1).exp_m1_p());
        return from_inverse_beta(value, df1, df2);
    }
    let mut a = 0.5 * df1;
    let mut b = 0.5 * df2;
    let reflected = a > b;
    if reflected {
        std::mem::swap(&mut a, &mut b);
    }
    let total = a + b;
    let p = a / total;
    let sigma = (1.0 / a + 1.0 / b).sqrt();
    let correction = stirling_correction(total) - stirling_correction(a) - stirling_correction(b);
    // For p<=1/2 the compensated loss avoids cancellation even for very
    // large and unequal df. Reflection changes only the log-density argument.
    let density = |z: f64| {
        let y = sigma * if reflected { -z } else { z };
        0.3989422804014327 * (correction - total * logistic_loss(y, p)).exp_p()
    };
    let upper_tail = |z: f64| integrate_positive_tail(|s| density(z + s));
    let target = 1.0 - 0.99;
    let mut low = -8.0;
    let mut high = 8.0;
    while upper_tail(low) < target {
        low *= 2.0;
    }
    while upper_tail(high) > target {
        high *= 2.0;
    }
    for _ in 0..80 {
        let mid = low + (high - low) * 0.5;
        if mid == low || mid == high {
            break;
        }
        if upper_tail(mid) > target {
            low = mid;
        } else {
            high = mid;
        }
    }
    from_inverse_beta((sigma * (low + 0.5 * (high - low))).exp_p(), df1, df2)
}

// Assemble the F quantile through its inverse-beta representation, as qf
// does. Preserve this common operation order even when the tail inversion
// above is more accurate: near-one beta quantiles affect the last output bits.
fn from_inverse_beta(value: f64, df1: f64, df2: f64) -> f64 {
    let odds = (df1 / df2) * value;
    let beta = if odds < 1.0 {
        1.0 - odds / (1.0 + odds)
    } else {
        1.0 / (1.0 + odds)
    };
    (1.0 / beta - 1.0) * (df2 / df1)
}

