//! Ordinary beta normalization for conditioned posterior moments.
//! Exact endpoint checks retain representable complements close to one.
use super::gamma_normalization::gamma_positive;
pub(super) fn beta_positive(a: f64, b: f64) -> f64 {
    if a + b < 170.0 {
        (1.0 / gamma_positive(a + b)) * gamma_positive(a) * gamma_positive(b)
    } else {
        statrs::function::beta::beta(a, b)
    }
}
// Moment normalization intentionally follows the source's gamma-product
// arithmetic. Probabilities use the stable recurrence for the integer/half
// shapes of Student-t moments, avoiding cancellation between large gammas.
fn probability_beta(mut a: f64, mut b: f64) -> f64 {
    if a < b {
        std::mem::swap(&mut a, &mut b);
    }
    if a <= 256.0 && (b == 0.5 || b == 1.5) {
        let twice = 2.0 * a;
        if twice.fract() == 0.0 {
            let integer = a.fract() == 0.0;
            let mut value = if integer { 2.0 } else { std::f64::consts::PI };
            let start = if integer { 1.0 } else { 0.5 };
            let mut shape = start;
            while shape < a {
                value *= shape / (shape + 0.5);
                shape += 1.0;
            }
            return if b == 1.5 {
                value / (twice + 1.0)
            } else {
                value
            };
        }
    }
    beta_positive(a, b)
}
fn fraction(a: f64, b: f64, x: f64) -> f64 {
    let guard = |v: f64| {
        if v.abs() < 1e-300 {
            1e-300_f64.copysign(v)
        } else {
            v
        }
    };
    let mut c = 1.0;
    let mut d = 1.0 / guard(1.0 - (a + b) * x / (a + 1.0));
    let mut value = d;
    for m in 1..=1024 {
        let m = f64::from(m);
        let twice = 2.0 * m;
        for (position, coefficient) in [
            m * (b - m) * x / ((a + twice - 1.0) * (a + twice)),
            -(a + m) * (a + b + m) * x / ((a + twice) * (a + twice + 1.0)),
        ]
        .into_iter()
        .enumerate()
        {
            d = 1.0 / guard(1.0 + coefficient * d);
            c = guard(1.0 + coefficient / c);
            let delta = c * d;
            value *= delta;
            if position == 1 && (delta - 1.0).abs() <= f64::EPSILON {
                return value;
            }
        }
    }
    f64::NAN
}
fn small_shape_series(a: f64, b: f64, x: f64) -> f64 {
    let mut coefficient = 1.0;
    let mut sum = 1.0 / a;
    let mut correction = 0.0;
    for n in 1..=10000 {
        let n = f64::from(n);
        coefficient *= (n - b) * x / n;
        let term = coefficient / (a + n);
        let adjusted = term - correction;
        let next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
        // Bound the remaining geometric tail, rather than only its next term.
        if term.abs() <= 0.125 * f64::EPSILON * (1.0 - x) * sum.abs() {
            return (x.powf(a) * sum) / probability_beta(a, b);
        }
    }
    f64::NAN
}
pub(super) fn regularized_beta(a: f64, b: f64, x: f64) -> f64 {
    if a.is_nan() || b.is_nan() || x.is_nan() || a <= 0.0 || b <= 0.0 {
        return f64::NAN;
    }
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    if b == 1.0 {
        return x.powf(a);
    }
    if a == 1.0 {
        return -(b * (-x).ln_1p()).exp_m1();
    }
    if (b == 0.5 || b == 1.5) && x > 0.5 && x < 0.98 {
        return small_shape_series(a, b, x);
    }
    let reflect = x >= (a + 1.0) / (a + b + 2.0);
    let (a, b, x) = if reflect { (b, a, 1.0 - x) } else { (a, b, x) };
    let normalization = probability_beta(a, b);
    let mut front = (x.powf(a) * (1.0 - x).powf(b)) / normalization;
    if !front.is_finite() || front == 0.0 {
        front = (a * x.ln() + b * (-x).ln_1p() - normalization.ln()).exp();
    }
    let probability = (front * fraction(a, b, x)) / a;
    if reflect {
        1.0 - probability
    } else {
        probability
    }
}

pub(super) fn log_regularized_beta(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if x >= 1.0 {
        return 0.0;
    }
    let reflect = x >= (a + 1.0) / (a + b + 2.0);
    let (a, b, x) = if reflect { (b, a, 1.0 - x) } else { (a, b, x) };
    let log = a * x.ln() + b * (-x).ln_1p() - statrs::function::beta::ln_beta(a, b)
        + fraction(a, b, x).ln()
        - a.ln();
    if reflect { (-log.exp()).ln_1p() } else { log }
}
