//! Student t tails with stable small and large shape evaluation.
use super::gamma::{log_half_gamma_ratio_correction, stirling_correction};
use super::quadrature::integrate_positive_tail;

// Modified Lentz evaluation of the incomplete-beta continued fraction,
// DLMF 8.17(v): https://dlmf.nist.gov/8.17.v . Both complementary arguments
// are passed explicitly so a small t^2/(df+t^2) is not rounded away.
fn beta_fraction(a: f64, b: f64, x: f64, y: f64) -> Option<f64> {
    let guard = |v: f64| {
        if v.abs() < f64::MIN_POSITIVE / f64::EPSILON {
            (f64::MIN_POSITIVE / f64::EPSILON).copysign(v)
        } else {
            v
        }
    };
    let mut numerator = 1.0;
    let mut denominator = 1.0 / guard((a * y + 1.0 - b * x) / (a + 1.0));
    let mut value = denominator;
    for m in 1..=512 {
        let m = f64::from(m);
        let twice = 2.0 * m;
        let coefficients = [
            m * (b - m) * x / ((a + twice - 1.0) * (a + twice)),
            -(a + m) * (a + b + m) * x / ((a + twice) * (a + twice + 1.0)),
        ];
        let mut delta = 1.0;
        for coefficient in coefficients {
            numerator = guard(1.0 + coefficient / numerator);
            denominator = 1.0 / guard(1.0 + coefficient * denominator);
            delta = numerator * denominator;
            value *= delta;
        }
        if (delta - 1.0).abs() <= f64::EPSILON {
            return Some(value);
        }
    }
    None
}

pub(crate) fn student_t_upper_tail(q: f64, df: f64) -> Option<f64> {
    if !q.is_finite() || !df.is_finite() || df <= 0.0 {
        return None;
    }
    if q == 0.0 {
        return Some(0.5);
    }
    let a = 0.5 * df;
    if a == 0.0 {
        return Some(0.5);
    }
    let square = q * q;
    let ratio = square / df;
    let log_ratio = if ratio > 0.0 && ratio.is_finite() {
        ratio.ln()
    } else {
        2.0 * q.abs().ln() - df.ln()
    };
    let log_one_plus_ratio = if ratio.is_finite() && ratio > 0.0 {
        let square_error = q.mul_add(q, -square);
        let ratio_error = ((-ratio).mul_add(df, square) + square_error) / df;
        ratio.ln_1p() + ratio_error / (1.0 + ratio)
    } else if log_ratio > 0.0 {
        log_ratio + (-log_ratio).exp().ln_1p()
    } else {
        log_ratio.exp().ln_1p()
    };
    if df >= 64.0 {
        let tail = large_df_t_positive_tail(q.abs(), df, log_one_plus_ratio);
        return Some(if q < 0.0 { 1.0 - tail } else { tail }.clamp(0.0, 1.0));
    }
    let x = (-log_one_plus_ratio).exp();
    let y = if ratio.is_finite() && ratio > 0.0 {
        ratio / (1.0 + ratio)
    } else {
        -(-log_one_plus_ratio).exp_m1()
    };
    let log_beta = if a.fract() == 0.0 && a <= 10.0 {
        // B(m,1/2) = 2 product(k/(k+1/2)). These integer products
        // are exact in f64 through m=10; no log-gamma subtraction.
        let mut numerator = 2.0;
        let mut denominator = 1.0;
        for k in 1..a as usize {
            numerator *= (2 * k) as f64;
            denominator *= (2 * k + 1) as f64;
        }
        (numerator / denominator).ln()
    } else if a >= 0.5 {
        let correction = a.mul_add((0.5 / a).ln_1p(), -0.5) + stirling_correction(a + 0.5)
            - stirling_correction(a);
        0.5723649429247001 - 0.5 * a.ln() - correction
    } else {
        libm::lgamma(a) + libm::lgamma(0.5) - libm::lgamma(a + 0.5)
    };
    if x <= 0.9 {
        // Positive hypergeometric series for I_x(a,1/2), DLMF 8.17.7.
        // Compensated summation avoids cancellation in the tail branch.
        let mut term = 1.0;
        let mut sum = 1.0;
        let mut compensation = 0.0;
        for n in 0..4096 {
            let n = f64::from(n);
            term *= ((a + n) / (a + n + 1.0)) * ((n + 0.5) / (n + 1.0)) * x;
            let delta = term - compensation;
            let next = sum + delta;
            compensation = (next - sum) - delta;
            sum = next;
            if term <= sum * (0.25 * f64::EPSILON) {
                let tail = 0.5 * ((-a * log_one_plus_ratio - log_beta).exp() / a) * sum;
                return Some(if q < 0.0 { 1.0 - tail } else { tail }.clamp(0.0, 1.0));
            }
        }
        return None;
    }
    let front = (-a * log_one_plus_ratio + 0.5 * y.ln() - log_beta).exp();
    let tail = if x < (a + 1.0) / (a + 2.5) {
        0.5 * front * beta_fraction(a, 0.5, x, y)? / a
    } else {
        0.5 * (1.0 - front * beta_fraction(0.5, a, y, x)? / 0.5)
    };
    Some(if q < 0.0 { 1.0 - tail } else { tail }.clamp(0.0, 1.0))
}

// For df>=64, integrate the positive tail after scaling by its log-density
// slope. This avoids cancellation of gamma logs and near-one beta arguments
// as df grows. The density ratio stays near unity in the important region.
// Gauss quadrature: https://dlmf.nist.gov/3.5.v .
fn large_df_t_positive_tail(q: f64, df: f64, log_one_plus_ratio: f64) -> f64 {
    let ratio = (q / df) * q;
    let rate = if q <= 1.0 {
        1.0
    } else {
        ((1.0 + 1.0 / df) * q) / (1.0 + ratio)
    };
    let exponent = log_half_gamma_ratio_correction(0.5 * df)
        - (0.5 * df + 0.5) * log_one_plus_ratio
        - rate.ln();
    if exponent < -760.0 || !exponent.is_finite() {
        return 0.0;
    }
    let density_ratio = |s: f64| {
        let u = s / rate;
        let shift = ((u / df) * (2.0 * q + u)) / (1.0 + ratio);
        (-(0.5 * df + 0.5) * shift.ln_1p()).exp()
    };
    let integral = integrate_positive_tail(density_ratio);
    // For df>=64 the omitted relative tail beyond s=128 is below 1e-29
    // even in the limiting power-law regime; no large gamma values are formed.
    (0.3989422804014327 * integral) * exponent.exp()
}
