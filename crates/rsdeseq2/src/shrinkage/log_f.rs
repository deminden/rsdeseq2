//! Log-F probabilities and native truncated first moments.
use super::{beta_normalization::log_regularized_beta, expression_quadrature::integrate};
fn softplus(x: f64) -> f64 {
    if x > 0.0 {
        x + (-x).exp().ln_1p()
    } else {
        x.exp().ln_1p()
    }
}
pub(super) fn log_cdf(x: f64, d1: f64, d2: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let f = x.exp();
    if f == 0.0 {
        return f64::NEG_INFINITY;
    }
    if f == f64::INFINITY {
        return 0.0;
    }
    if d1.is_infinite() && d2.is_infinite() {
        return if f < 1.0 {
            f64::NEG_INFINITY
        } else if f == 1.0 {
            -std::f64::consts::LN_2
        } else {
            0.0
        };
    }
    if d2.is_infinite() {
        return super::gamma_probability::log_gamma_probability(d1 / 2.0, (d1 / 2.0) * f, true);
    }
    if d1.is_infinite() {
        return super::gamma_probability::log_gamma_probability(d2 / 2.0, (d2 / 2.0) / f, false);
    }
    let z = (d1 / d2).ln() + x;
    // Use the smaller logistic coordinate to retain representable complements.
    if z <= 0.0 {
        log_regularized_beta(d1 / 2.0, d2 / 2.0, z.exp() / (1.0 + z.exp()))
    } else {
        let log_tail = log_regularized_beta(d2 / 2.0, d1 / 2.0, (-z).exp() / (1.0 + (-z).exp()));
        (-log_tail.exp()).ln_1p()
    }
}
pub(super) fn log_density(x: f64, d1: f64, d2: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x == f64::INFINITY {
        return f64::NAN;
    }
    let f = x.exp();
    if d1.is_infinite() && d2.is_infinite() {
        return if f == 1.0 {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        };
    }
    if f == 0.0 {
        return if d1 < 2.0 {
            if x.is_infinite() {
                f64::NAN
            } else {
                f64::INFINITY
            }
        } else if d1 == 2.0 {
            x
        } else {
            f64::NEG_INFINITY
        };
    }
    if f == f64::INFINITY {
        return f64::NEG_INFINITY;
    }
    if d2.is_infinite() {
        let a = d1 / 2.0;
        if a >= 1e5 {
            return super::large_shape::gamma_log_density(f, a, a) + x;
        }
        return a * a.ln() - statrs::function::gamma::ln_gamma(a) + a * x - a * f;
    }
    if d1.is_infinite() {
        let b = d2 / 2.0;
        if b >= 1e5 {
            return super::large_shape::gamma_log_density(1.0 / f, b, b) - x;
        }
        return b * b.ln() - statrs::function::gamma::ln_gamma(b) - b * x - b / f;
    }
    let a = d1 / 2.0;
    let b = d2 / 2.0;
    let z = (d1 / d2).ln() + x;
    if super::large_shape::beta_large(a, b) {
        return super::large_shape::beta_logit_density(a, b, x);
    }
    if a.max(b) >= 1e5 {
        let numerator = if z > 0.0 {
            -b * z - (a + b) * (-z).exp().ln_1p()
        } else {
            a * z - (a + b) * z.exp().ln_1p()
        };
        return numerator - super::large_shape::log_beta(a, b);
    }
    a * z - (a + b) * softplus(z) - statrs::function::beta::ln_beta(a, b)
}
pub(super) fn truncated_mean(a: f64, b: f64, d1: f64, d2: f64) -> f64 {
    if a == b {
        return a;
    }
    if a.is_nan() || b.is_nan() || a > b {
        return f64::NAN;
    }
    let cdf = |x: f64| {
        if d1.is_infinite() && !d2.is_infinite() {
            super::gamma_probability::gamma_probability(d2 / 2.0, (d2 / 2.0) / x.exp(), false)
        } else if d2.is_infinite() && !d1.is_infinite() {
            super::gamma_probability::gamma_probability(d1 / 2.0, (d1 / 2.0) * x.exp(), true)
        } else if d1.is_finite() && d2.is_finite() {
            let numerator = d1 * x.exp();
            if numerator > d2 {
                let tail = log_regularized_beta(d2 / 2.0, d1 / 2.0, d2 / (d2 + numerator)).exp();
                // The ordinary beta complement uses two subtractions. Its
                // rounding is observable when a moment subtracts nearby CDFs.
                (0.5 - tail) + 0.5
            } else {
                log_regularized_beta(d1 / 2.0, d2 / 2.0, numerator / (d2 + numerator)).exp()
            }
        } else {
            log_cdf(x, d1, d2).exp()
        }
    };
    let probability = cdf(b) - cdf(a);
    if probability > 0.0 {
        let density = |x: f64| {
            let f = x.exp();
            if f.is_infinite()
                || (d1.is_infinite() && f > 0.0 && f * f == 0.0)
                || (f == 0.0 && d1 < 2.0)
            {
                f64::NAN
            } else {
                log_density(x, d1, d2).exp()
            }
        };
        let integrand = |x: f64| x * (density(x) / probability);
        let value = if a.is_infinite() || b.is_infinite() {
            super::log_f_quadrature::integrate_infinite(integrand, a, b)
        } else {
            integrate(integrand, a, b)
        };
        if let Some(value) = value {
            return value;
        }
    }
    if a.is_infinite() || b.is_infinite() {
        return (a + b) / 2.0;
    }
    // The source retries a scaled numerator/denominator before midpoint fallback.
    let scale = log_density(a, d1, d2).min(log_density(b, d1, d2));
    if scale.is_finite() {
        let numerator = integrate(|x| x * (log_density(x, d1, d2) - scale).exp(), a, b);
        let denominator = integrate(|x| (log_density(x, d1, d2) - scale).exp(), a, b);
        if let (Some(n), Some(d)) = (numerator, denominator)
            && d > 0.0
        {
            return n / d;
        }
    }
    (a + b) / 2.0
}
