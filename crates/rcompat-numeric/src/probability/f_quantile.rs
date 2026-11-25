//! Fixed Cook's F quantile from a stable beta-logistic density.
use super::gamma::{chi_square_quantile, log1pmx, stirling_correction};
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

/// Degrees of freedom above which `f_99_quantile` uses a chi-square limit.
const LARGE_DF: f64 = 400_000.0;

/// The 0.99 quantile of the F distribution with `df1` and `df2` degrees of
/// freedom, as DESeq2 uses for its Cook's distance cutoff.
pub fn f_99_quantile(df1: f64, df2: f64) -> f64 {
    // Beyond 400000 degrees of freedom the quantile is taken from a
    // chi-square limit of the F distribution (docs/provenance/specs/
    // f_quantile_large_df.md): chi-square(df1) / df1 as df2 grows, and
    // df2 / chi-square(df2) as df1 grows. When both are large, the first
    // limit is used whenever df1 <= df2.
    if df1 <= df2 && df2 > LARGE_DF {
        return chi_square_quantile(0.99, df1, true) / df1;
    }
    if df1 > LARGE_DF {
        return df2 / chi_square_quantile(0.99, df2, false);
    }
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

#[cfg(test)]
mod tests {
    use super::f_99_quantile;

    // qf(0.99, df1, df2) from R 4.6.2 (scripts/references/generate_f_quantile_reference.R),
    // for df1 up to 200 and df2 up to 100,000.
    const REFERENCE: &str = include_str!("../../tests/data/f_quantile_r462.tsv");

    // qf(0.99, df1, df2) from R 4.6.2
    // (scripts/references/generate_f_quantile_large_df_reference.R) with one
    // of df1, df2 between 400001 and 1e7, where the chi-square limits apply.
    const LARGE_DF_REFERENCE: &str = include_str!("../../tests/data/f_quantile_large_df_r462.tsv");

    /// Largest relative difference from the reference, with its df, and the
    /// number of rows compared.
    fn worst_relative_difference(reference: &str) -> ((f64, f64, f64), usize) {
        let mut worst = (0.0_f64, 0.0, 0.0);
        let mut rows = 0;
        for line in reference.lines().filter(|l| !l.starts_with('#')).skip(1) {
            let fields: Vec<&str> = line.split('\t').collect();
            let df1: f64 = fields[0].parse().unwrap();
            let df2: f64 = fields[1].parse().unwrap();
            let expected = f64::from_bits(u64::from_str_radix(fields[2], 16).unwrap());
            let actual = f_99_quantile(df1, df2);
            let relative = ((actual - expected) / expected).abs();
            if relative > worst.0 {
                worst = (relative, df1, df2);
            }
            rows += 1;
        }
        eprintln!(
            "largest relative difference {:e} at df1 {} df2 {}",
            worst.0, worst.1, worst.2
        );
        (worst, rows)
    }

    #[test]
    fn matches_r_qf_over_the_supported_degrees_of_freedom() {
        let (worst, rows) = worst_relative_difference(REFERENCE);
        assert_eq!(rows, 22 * 26);
        assert!(worst.0 <= TOLERANCE, "{worst:?}");
    }

    #[test]
    fn matches_r_qf_in_the_large_degrees_of_freedom_limits() {
        let (worst, rows) = worst_relative_difference(LARGE_DF_REFERENCE);
        assert_eq!(rows, 2 * 9 * 6 + 4);
        assert!(worst.0 <= LARGE_DF_TOLERANCE, "{worst:?}");
    }

    const TOLERANCE: f64 = 2e-14;
    const LARGE_DF_TOLERANCE: f64 = 1.2e-14;
}
