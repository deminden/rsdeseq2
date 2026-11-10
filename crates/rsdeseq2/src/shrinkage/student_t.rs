//! Student-t likelihood leaves for finite uniform components.
use super::beta_normalization::{beta_positive as beta, regularized_beta};
#[cfg(test)]
use super::gamma_normalization::gamma_positive;
use crate::probability::student_t_upper_tail;

pub(super) fn cdf(q: f64, df: f64) -> f64 {
    if q.is_nan() {
        return q;
    }
    if q == f64::NEG_INFINITY {
        return 0.0;
    }
    if q == f64::INFINITY {
        return 1.0;
    }
    student_t_upper_tail(-q, df).unwrap_or(f64::NAN)
}
// The log-series and ordinary moment integrals retain distinct normalization
// schedules. Replacing this with the direct moment-beta product changes the
// subtraction of nearby log-CDF endpoints.
fn log_series_normalization(a: f64) -> f64 {
    let b = 0.5;
    let normalizer = if a + b < 170.0 {
        (1.0 / libm::tgamma(a + b)) * libm::tgamma(a) * libm::tgamma(b)
    } else {
        statrs::function::beta::beta(a, b)
    };
    (1.0 / (a * normalizer)).ln()
}
pub(super) fn log_cdf(q: f64, df: f64) -> f64 {
    if q == 0.0 {
        return -std::f64::consts::LN_2;
    }
    let x = 1.0 / (1.0 + (q / df) * q);
    let a = df / 2.0;
    let log_beta_probability = if x <= 0.9 {
        // Integrate the binomial series for (1-t)^(-1/2), retaining
        // its leading beta normalization separately from the remainder.
        let mut remainder = 0.0;
        let mut term = 1.0;
        for n in 1..4096 {
            let n = f64::from(n);
            term *= ((n - 0.5) / n) * x;
            let contribution = term / (a + n);
            remainder += contribution;
            if contribution.abs() <= f64::EPSILON / a {
                break;
            }
        }
        a * x.ln() + log_series_normalization(a) + (a * remainder).ln_1p()
    } else {
        (2.0 * cdf(-q.abs(), df)).ln()
    };
    if q < 0.0 {
        log_beta_probability - std::f64::consts::LN_2
    } else {
        (-0.5 * log_beta_probability.exp()).ln_1p()
    }
}
pub(super) fn log_density(q: f64, df: f64) -> f64 {
    let ratio = q * q / df;
    let log_ratio = if ratio.is_finite() {
        ratio.ln_1p()
    } else {
        2.0 * q.abs().ln() - df.ln()
    };
    libm::lgamma((df + 1.0) / 2.0)
        - libm::lgamma(df / 2.0)
        - 0.5 * (df.ln() + std::f64::consts::PI.ln())
        - ((df + 1.0) / 2.0) * log_ratio
}
fn integral(q: f64, df: f64, order: i32) -> f64 {
    let a = 0.5 * (df - f64::from(order));
    let b = 0.5 * f64::from(order + 1);
    if a <= 0.0 {
        return f64::NAN;
    }
    let x = df / (df + q * q);
    let full = beta(a, b);
    let incomplete = regularized_beta(a, b, x) * full;
    if q > 0.0 || order % 2 != 0 {
        incomplete
    } else {
        2.0 * full - incomplete
    }
}
fn one_sided_moment(q: f64, df: f64, order: i32) -> f64 {
    (df.powf(0.5 * f64::from(order)) * integral(q, df, order)) / integral(q, df, 0)
}
pub(super) fn truncated_moment(mut a: f64, mut b: f64, df: f64, order: i32) -> f64 {
    let reflect = a > 0.0 && b > 0.0;
    if reflect {
        (a, b) = (-b, -a);
    }
    let moment = if a == b {
        a.powi(order)
    } else {
        (cdf(-a, df) * one_sided_moment(a, df, order)
            - cdf(-b, df) * one_sided_moment(b, df, order))
            / (cdf(b, df) - cdf(a, df))
    };
    if reflect && order % 2 != 0 {
        -moment
    } else {
        moment
    }
}

#[cfg(test)]
mod normalization_tests {
    use super::*;
    use serde::Deserialize;
    #[derive(Deserialize)]
    struct Reference {
        gamma_arguments: Vec<String>,
        gamma_values: Vec<String>,
        cases: Vec<Case>,
    }
    #[derive(Deserialize)]
    struct Case {
        a: String,
        b: String,
        p: String,
        #[serde(default)]
        beta: String,
        regularized: String,
    }
    fn number(s: &str) -> f64 {
        if let Some(s) = s.strip_prefix("bits:") {
            f64::from_bits(u64::from_str_radix(s, 16).unwrap())
        } else {
            s.parse().unwrap()
        }
    }
    fn close(a: f64, b: f64, context: String, failures: &mut Vec<String>) {
        if a == b {
            return;
        }
        if !a.is_finite() || !b.is_finite() || a.to_bits().abs_diff(b.to_bits()) > 16 {
            failures.push(format!(
                "{context}: {a:.17e} != {b:.17e}; ULP {}",
                a.to_bits().abs_diff(b.to_bits())
            ));
        }
    }
    #[test]
    fn component_normalization_matches_independent_original_leaves() {
        let r: Reference = serde_json::from_str(include_str!(
            "../../tests/data/student_component_leaves_r462.json"
        ))
        .unwrap();
        let mut failures = vec![];
        for (x, y) in r.gamma_arguments.iter().zip(&r.gamma_values) {
            let x = number(x);
            close(
                gamma_positive(x),
                number(y),
                format!("Gamma({x})"),
                &mut failures,
            );
        }
        for c in r.cases {
            let a = number(&c.a);
            let b = number(&c.b);
            let p = number(&c.p);
            close(
                beta(a, b),
                number(&c.beta),
                format!("B({a},{b})"),
                &mut failures,
            );
            // The original pbeta differs from independent 90-digit truth by
            // up to 274 ULP in this panel. Preserve relative source parity,
            // and check tight probability accuracy against that truth below.
            let original = number(&c.regularized);
            assert!((regularized_beta(a, b, p) - original).abs() <= 1e-13 * original.abs());
        }
        assert!(
            failures.is_empty(),
            "{} normalization mismatches:\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
    #[test]
    fn beta_probabilities_match_independent_mathematical_reference() {
        #[derive(Deserialize)]
        struct MathematicalReference {
            cases: Vec<Case>,
        }
        let reference: MathematicalReference = serde_json::from_str(include_str!(
            "../../tests/data/student_beta_mathematical.json"
        ))
        .unwrap();
        let mut failures = vec![];
        for c in reference.cases {
            let a = number(&c.a);
            let b = number(&c.b);
            let x = number(&c.p);
            close(
                regularized_beta(a, b, x),
                number(&c.regularized),
                format!("I({x},{a},{b})"),
                &mut failures,
            );
        }
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    }
}
