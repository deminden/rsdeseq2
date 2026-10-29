//! Stable gamma normalization and chi-square tails.
use super::quadrature::integrate_positive_tail;

// log Gamma(a+1/2) - log Gamma(a) - (1/2)log(a), a>=32.
// Bernoulli-polynomial gamma-ratio expansion, DLMF 5.11.13.
pub(super) fn log_half_gamma_ratio_correction(a: f64) -> f64 {
    let inverse = 1.0 / a;
    let square = inverse * inverse;
    inverse
        * (-0.125
            + square
                * (1.0 / 192.0
                    + square
                        * (-1.0 / 640.0 + square * (17.0 / 14336.0 - square * (31.0 / 18432.0)))))
}

// log Gamma(a) minus its leading Stirling expression, a>=1/2. Shift small
// shapes using a positive atanh series, rather than subtracting large logs.
pub(super) fn stirling_correction(mut a: f64) -> f64 {
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

// Stable remainder log(1+x)-x, needed in gamma and beta log densities.
pub(super) fn log1pmx(x: f64) -> f64 {
    if x.abs() > 0.125 {
        return x.ln_1p() - x;
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

pub(crate) fn large_shape_chi_square_tail(statistic: f64, df: f64) -> f64 {
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
        0.3989422804014327 * (-correction + a * log1pmx(u) - u.ln_1p()).exp()
    });
    (if lower { 1.0 - integral } else { integral }).clamp(0.0, 1.0)
}

// Chi-square inverse from the same standalone tails used for inference. The
// argument is an upper-tail probability, including the small complementary
// probability used by DESeq2/R's large-df F approximation.
pub(super) fn chi_square_upper_quantile(upper: f64, df: f64) -> f64 {
    use statrs::distribution::{ChiSquared, ContinuousCDF};
    if df == 2.0 {
        return -2.0 * upper.ln();
    }
    let distribution = match ChiSquared::new(df) {
        Ok(value) => value,
        Err(_) => return f64::NAN,
    };
    let survival = |x: f64| {
        if df >= 64.0 {
            large_shape_chi_square_tail(x, df)
        } else {
            distribution.sf(x)
        }
    };
    let mut low = -740.0;
    let mut high = (df + 1.0).ln() + 8.0;
    for _ in 0..96 {
        let mid = low + (high - low) * 0.5;
        if mid == low || mid == high {
            break;
        }
        if survival(mid.exp()) > upper {
            low = mid;
        } else {
            high = mid;
        }
    }
    (low + 0.5 * (high - low)).exp()
}
