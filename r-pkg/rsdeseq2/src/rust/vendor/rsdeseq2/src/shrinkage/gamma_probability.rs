//! Positive incomplete-gamma log probabilities from series and continued fractions.
fn log_complement(log: f64) -> f64 {
    if log < -std::f64::consts::LN_2 {
        (-log.exp()).ln_1p()
    } else {
        (-log.exp_m1()).ln()
    }
}

pub(super) fn log_gamma_probability(a: f64, x: f64, lower: bool) -> f64 {
    if a <= 0.0 || a.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    if x <= 0.0 {
        return if lower { f64::NEG_INFINITY } else { 0.0 };
    }
    if x.is_infinite() {
        return if lower { 0.0 } else { f64::NEG_INFINITY };
    }
    if a >= 1e5 {
        return super::large_shape::gamma_log_probability(a, x, lower);
    }
    let front = a * x.ln() - x - statrs::function::gamma::ln_gamma(a);
    let (value, direct_lower) = if x < a + 1.0 {
        let mut term = 1.0 / a;
        let mut sum = term;
        for n in 1..=10000 {
            term *= x / (a + f64::from(n));
            sum += term;
            if term.abs() <= f64::EPSILON * sum.abs() {
                return if lower {
                    front + sum.ln()
                } else {
                    log_complement(front + sum.ln())
                };
            }
        }
        (f64::NAN, true)
    } else {
        let guard = |x: f64| {
            if x.abs() < 1e-300 {
                1e-300_f64.copysign(x)
            } else {
                x
            }
        };
        let mut b = x + 1.0 - a;
        let mut c = 1e300;
        let mut d = 1.0 / guard(b);
        let mut fraction = d;
        for n in 1..=10000 {
            let n = f64::from(n);
            let coefficient = -n * (n - a);
            b += 2.0;
            d = 1.0 / guard(coefficient * d + b);
            c = guard(b + coefficient / c);
            let delta = d * c;
            fraction *= delta;
            if (delta - 1.0).abs() <= f64::EPSILON {
                let log = front + fraction.ln();
                return if lower { log_complement(log) } else { log };
            }
        }
        (f64::NAN, false)
    };
    if lower == direct_lower {
        value
    } else {
        log_complement(value)
    }
}

/// Ordinary small-argument gamma probabilities retain the series' operation
/// order. In particular, forming Q through exp(log(Q)) changes cancellation in
/// the original log-F moment denominator near one.
pub(super) fn gamma_probability(a: f64, x: f64, lower: bool) -> f64 {
    if x > 0.0 && x < 1.0 {
        // P(a,x) = x^a/Gamma(a+1) * (1 + sum_n a*(-x)^n/(n!*(a+n))).
        let mut coefficient = a;
        let mut sum = 0.0;
        for n in 1..=10000 {
            let n = f64::from(n);
            coefficient *= -x / n;
            let term = coefficient / (a + n);
            sum += term;
            if term.abs() <= f64::EPSILON * sum.abs() {
                let log_scale = a * x.ln() - statrs::function::gamma::ln_gamma(a + 1.0);
                return if lower {
                    (1.0 + sum) * log_scale.exp()
                } else {
                    let delta = log_scale.exp_m1();
                    -(sum + delta + sum * delta)
                };
            }
        }
        f64::NAN
    } else {
        log_gamma_probability(a, x, lower).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn gamma_log_tail_retains_the_analytic_exponential_at_extreme_arguments() {
        for x in [0.001, 0.1, 2.0, 1000.0, 1e22] {
            let q = log_gamma_probability(1.0, x, false);
            assert!((q + x).abs() <= 8.0 * f64::EPSILON * x.max(1.0));
        }
        let x = 0.001;
        let p = log_gamma_probability(1.0, x, true);
        assert!((p - (-(-x).exp_m1()).ln()).abs() < 1e-13);
    }
}
