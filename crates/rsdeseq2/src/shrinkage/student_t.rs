//! Student-t likelihood leaves for finite uniform components.
use crate::probability::student_t_upper_tail;
use statrs::function::beta::checked_beta_reg;

fn beta(a: f64, b: f64) -> f64 {
    if a + b < 170.0 {
        (1.0 / libm::tgamma(a + b)) * libm::tgamma(a) * libm::tgamma(b)
    } else {
        statrs::function::beta::beta(a, b)
    }
}

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
        a * x.ln() + (1.0 / (a * beta(a, 0.5))).ln() + (a * remainder).ln_1p()
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
    let incomplete = checked_beta_reg(a, b, x).unwrap_or(f64::NAN) * full;
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
