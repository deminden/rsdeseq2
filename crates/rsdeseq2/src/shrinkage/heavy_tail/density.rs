//! Independent density equations for imported regression likelihoods.
use super::super::expression_quadrature::integrate;
use crate::probability::{log_half_gamma_ratio_correction, log1pmx, stirling_correction};
const LOG_2_PI: f64 = 1.8378770664093453;
pub(super) fn normal(x: f64, mean: f64, sd: f64) -> f64 {
    if x.is_nan() || mean.is_nan() || sd.is_nan() || sd < 0.0 {
        return f64::NAN;
    }
    if x == mean && !x.is_finite() {
        return f64::NAN;
    }
    if sd == 0.0 {
        return if x == mean {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        };
    }
    let z = (x - mean) / sd;
    -0.5 * z * z - 0.5 * LOG_2_PI - sd.ln()
}
fn stirling(x: f64) -> f64 {
    if x >= 0.5 {
        stirling_correction(x)
    } else {
        libm::lgamma(x) - (x - 0.5) * x.ln() + x - 0.5 * LOG_2_PI
    }
}
fn deviance(x: f64, mean: f64) -> f64 {
    if x == 0.0 {
        return mean;
    }
    if mean == 0.0 {
        return f64::INFINITY;
    }
    if (x - mean).abs() < 0.1 * (x + mean) {
        let v = (x - mean) / (x + mean);
        let square = v * v;
        let mut term = 2.0 * x * v;
        let mut sum = (x - mean) * v;
        for j in 1..1024 {
            term *= square;
            let next = sum + term / f64::from(2 * j + 1);
            if next == sum {
                return next;
            }
            sum = next;
        }
        sum
    } else {
        x * (x / mean).ln() + mean - x
    }
}
fn count(x: f64) -> Option<f64> {
    if !x.is_finite() || x < 0.0 || (x - x.round()).abs() > 1e-7 * x.abs().max(1.0) {
        None
    } else {
        Some(x.round())
    }
}
pub(super) fn binomial(x: f64, size: f64, p: f64) -> f64 {
    if x.is_nan() || size.is_nan() || p.is_nan() || !(0.0..=1.0).contains(&p) {
        return f64::NAN;
    }
    let Some(n) = count(size) else {
        return f64::NAN;
    };
    let Some(x) = count(x) else {
        return f64::NEG_INFINITY;
    };
    if x > n {
        return f64::NEG_INFINITY;
    }
    if x == 0.0 {
        return if n == 0.0 { 0.0 } else { n * (-p).ln_1p() };
    }
    if x == n {
        return n * p.ln();
    }
    stirling(n)
        - stirling(x)
        - stirling(n - x)
        - deviance(x, n * p)
        - deviance(n - x, n * (1.0 - p))
        - 0.5 * (LOG_2_PI + x.ln() + (-x / n).ln_1p())
}
pub(super) fn negative_binomial(x: f64, size: f64, mu: f64) -> f64 {
    if x.is_nan() || size.is_nan() || mu.is_nan() || size < 0.0 || mu < 0.0 {
        return f64::NAN;
    }
    let Some(x) = count(x) else {
        return f64::NEG_INFINITY;
    };
    if size == 0.0 || mu == 0.0 {
        return if x == 0.0 { 0.0 } else { f64::NEG_INFINITY };
    }
    if mu.is_infinite() {
        return f64::NEG_INFINITY;
    }
    if size.is_infinite() {
        return if x == 0.0 {
            -mu
        } else {
            -stirling(x) - deviance(x, mu) - 0.5 * (LOG_2_PI + x.ln())
        };
    }
    if x == 0.0 {
        return -size * (mu / size).ln_1p();
    }
    let n = x + size;
    let ratio = mu / size;
    let p = 1.0 / (1.0 + ratio);
    let q = if ratio.is_infinite() {
        1.0
    } else {
        ratio / (1.0 + ratio)
    };
    stirling(n)
        - stirling(x)
        - stirling(size)
        - deviance(x, n * q)
        - deviance(size, n * p)
        - 0.5 * (LOG_2_PI + x.ln())
        + 0.5 * (size / n).ln()
}
fn central_student(x: f64, df: f64) -> f64 {
    if df >= 64.0 {
        let ratio = x * x / df;
        let logarithm = if ratio.is_finite() {
            ratio.ln_1p()
        } else {
            2.0 * x.abs().ln() - df.ln()
        };
        log_half_gamma_ratio_correction(df / 2.0) - 0.5 * LOG_2_PI - 0.5 * (df + 1.0) * logarithm
    } else {
        super::super::student_t::log_density(x, df)
    }
}
fn peak_mass(df: f64, mode: f64) -> f64 {
    let f = |u: f64| {
        if u <= -mode {
            0.0
        } else {
            (df * log1pmx(u / mode) - 0.5 * u * u).exp()
        }
    };
    // The log density has curvature <= -1. Beyond 16 mode-distance units
    // its mass is bounded by the negligible standard-normal tail.
    integrate(f, -mode.min(16.0), 0.0).unwrap_or(f64::NAN)
        + integrate(f, 0.0, 16.0).unwrap_or(f64::NAN)
}
pub(super) fn student(x: f64, df: f64, ncp: f64) -> f64 {
    if x.is_nan() || df.is_nan() || ncp.is_nan() || df <= 0.0 {
        return f64::NAN;
    }
    if df.is_infinite() {
        return normal(x, ncp, 1.0);
    }
    if !x.is_finite() {
        return f64::NEG_INFINITY;
    }
    let central = central_student(x, df);
    if ncp == 0.0 {
        return central;
    }
    let root = df.sqrt();
    let k = ncp * (x / (x.hypot(root)));
    let length = k.hypot(2.0 * root);
    let mode = if k < 0.0 {
        df / (length / 2.0 - k / 2.0)
    } else {
        k / 2.0 + length / 2.0
    };
    let peak_difference = df * (k / (2.0 * root)).asinh() + 0.5 * k * mode;
    central - 0.5 * ncp * ncp + peak_difference + (peak_mass(df, mode) / peak_mass(df, root)).ln()
}
