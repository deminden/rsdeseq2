use crate::support as common;
use crate::support::formula_reference;
use common::{r_number as number, r_numbers as numbers};
use rsdeseq2::shrinkage::{FixedMixture, TruncatedNormalMixture, UniformMixture};
use rsdeseq2::{NumericProfile, numeric_profile};
use serde::Deserialize;
use std::collections::BTreeMap;

#[derive(Deserialize)]
struct References {
    source: Source,
    cases: Vec<Case>,
}
#[derive(Deserialize)]
struct Source {
    dnorm_probe: Option<Vec<String>>,
}
#[derive(Deserialize)]
struct Case {
    name: String,
    df: Option<String>,
    weights: Vec<String>,
    means: Option<Vec<String>>,
    sds: Option<Vec<String>>,
    lower: Vec<String>,
    upper: Vec<String>,
    estimates: Vec<String>,
    standard_errors: Vec<String>,
    thresholds: Vec<String>,
    log_likelihoods: Vec<String>,
    posterior_weights: Vec<String>,
    posterior_mean: Vec<String>,
    posterior_sd: Vec<String>,
    cdf: Vec<String>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum FixtureKind {
    Uniform,
    TruncatedNormal,
    StudentT,
}

#[derive(Deserialize)]
struct MathematicalReference {
    cases: BTreeMap<String, Vec<MathematicalRow>>,
}

#[derive(Deserialize)]
struct MathematicalRow {
    estimate: String,
    standard_error: String,
    log_likelihoods: Vec<String>,
    posterior_weights: Vec<String>,
    posterior_mean: String,
    posterior_sd: String,
    cdf: Vec<String>,
}
fn close(a: f64, b: f64, absolute: f64, relative: f64, context: &str, failures: &mut Vec<String>) {
    if (b.is_nan() && a.is_nan()) || a == b {
        return;
    }
    if a.is_finite() && b.is_finite() && (a - b).abs() <= absolute + relative * b.abs() {
        return;
    }
    failures.push(format!("{context}: {a:.17e} != {b:.17e}"));
}
#[test]
fn bounded_mixture_likelihoods_posteriors_and_cdfs_match_independent_ashr() {
    check_bounded_fixtures(&[
        (
            FixtureKind::Uniform,
            formula_reference!("uniform_posterior_r462"),
        ),
        (
            FixtureKind::TruncatedNormal,
            include_str!("data/truncated_posterior_r462.json"),
        ),
    ]);
}

#[test]
fn student_uniform_mixtures_match_independent_posterior_references() {
    check_bounded_fixtures(&[(
        FixtureKind::StudentT,
        include_str!("data/student_uniform_posterior_r462.json"),
    )]);
}

fn check_bounded_fixtures(fixtures: &[(FixtureKind, &str)]) {
    let mut failures = vec![];
    #[derive(Deserialize)]
    struct CdfReference {
        case: String,
        row: usize,
        cutoff: usize,
        estimate: String,
        se: String,
        threshold: String,
        cdf: String,
    }
    let portable = numeric_profile() == NumericProfile::Portable;
    let calibration = (!portable
        && fixtures
            .iter()
            .any(|(kind, _)| *kind == FixtureKind::Uniform))
    .then(|| std::env::var_os("RSDESEQ2_UNIFORM_CDF_REFERENCE"))
    .flatten()
    .map(|path| {
        csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(path)
            .unwrap()
            .deserialize::<CdfReference>()
            .map(|r| {
                let row = r.unwrap();
                ((row.case.clone(), row.row, row.cutoff), row)
            })
            .collect::<BTreeMap<_, _>>()
    });
    let mut calibrated = 0usize;
    for &(kind, fixture) in fixtures {
        let mathematical = match (portable, kind) {
            (true, FixtureKind::Uniform) => {
                Some(include_str!("data/gaussian_uniform_mathematical.json"))
            }
            (true, FixtureKind::StudentT) => {
                Some(include_str!("data/student_uniform_mathematical.json"))
            }
            _ => None,
        }
        .map(|source| serde_json::from_str::<MathematicalReference>(source).unwrap());
        let mut mathematical_cases = 0usize;
        let r: References = serde_json::from_str(fixture).unwrap();
        if !portable && let Some(probe) = r.source.dnorm_probe {
            assert_eq!(probe.len(), 1);
            let point = UniformMixture::new(vec![1.0], vec![0.0], vec![0.0]).unwrap();
            let density = point.log_likelihoods(37.0 / 0.15, 1.0).unwrap()[0];
            assert_eq!(
                density.to_bits(),
                number(&probe[0]).to_bits(),
                "original dnorm probe"
            );
        }
        for c in r.cases {
            let mathematical_rows = mathematical
                .as_ref()
                .map(|reference| &reference.cases[&c.name]);
            if let Some(rows) = mathematical_rows {
                assert_eq!(rows.len(), c.estimates.len());
                mathematical_cases += 1;
            }
            let prior = if let Some(means) = &c.means {
                FixedMixture::TruncatedNormal(
                    TruncatedNormalMixture::new(
                        numbers(&c.weights),
                        numbers(means),
                        numbers(c.sds.as_ref().unwrap()),
                        numbers(&c.lower),
                        numbers(&c.upper),
                    )
                    .unwrap(),
                )
            } else {
                let mut p =
                    UniformMixture::new(numbers(&c.weights), numbers(&c.lower), numbers(&c.upper))
                        .unwrap();
                if let Some(df) = &c.df {
                    p = p.with_student_t(number(df)).unwrap();
                }
                FixedMixture::Uniform(p)
            };
            let thresholds = numbers(&c.thresholds);
            let posteriors = prior
                .posteriors(&numbers(&c.estimates), &numbers(&c.standard_errors))
                .unwrap();
            for (i, (&estimate, se)) in numbers(&c.estimates)
                .iter()
                .zip(numbers(&c.standard_errors))
                .enumerate()
            {
                let context = format!("{} row{i} x{estimate} se{se}", c.name);
                let mathematical_row = mathematical_rows.map(|rows| &rows[i]);
                let (expected_mean, expected_sd) = if let Some(row) = mathematical_row {
                    assert_eq!(number(&row.estimate).to_bits(), estimate.to_bits());
                    assert_eq!(number(&row.standard_error).to_bits(), se.to_bits());
                    assert_eq!(row.log_likelihoods.len(), c.weights.len());
                    assert_eq!(row.posterior_weights.len(), c.weights.len());
                    assert_eq!(row.cdf.len(), thresholds.len());
                    (number(&row.posterior_mean), number(&row.posterior_sd))
                } else {
                    (number(&c.posterior_mean[i]), number(&c.posterior_sd[i]))
                };
                let likelihood = prior.log_likelihoods(estimate, se).unwrap();
                for (j, &a) in likelihood.iter().enumerate() {
                    let expected = mathematical_row.map_or_else(
                        || &c.log_likelihoods[i * c.weights.len() + j],
                        |row| &row.log_likelihoods[j],
                    );
                    close(
                        a,
                        number(expected),
                        1e-8,
                        1e-12,
                        &format!("{context} likelihood{j}"),
                        &mut failures,
                    );
                }
                let p = &posteriors[i];
                close(
                    p.mean(),
                    expected_mean,
                    2e-6,
                    0.,
                    &format!("{context} mean"),
                    &mut failures,
                );
                close(
                    p.standard_deviation(),
                    expected_sd,
                    2e-6,
                    0.,
                    &format!("{context} SD"),
                    &mut failures,
                );
                for (j, &a) in p.weights().iter().enumerate() {
                    let expected = mathematical_row.map_or_else(
                        || &c.posterior_weights[i * c.weights.len() + j],
                        |row| &row.posterior_weights[j],
                    );
                    close(
                        a,
                        number(expected),
                        2e-8,
                        0.,
                        &format!("{context} weight{j}"),
                        &mut failures,
                    );
                }
                for (j, &threshold) in thresholds.iter().enumerate() {
                    let expected = if let Some(row) = mathematical_row {
                        Some(number(&row.cdf[j]))
                    } else if kind == FixtureKind::Uniform {
                        calibration.as_ref().map(|rows| {
                            let row = &rows[&(c.name.clone(), i, j)];
                            assert_eq!(number(&row.estimate).to_bits(), estimate.to_bits());
                            assert_eq!(number(&row.se).to_bits(), se.to_bits());
                            assert_eq!(number(&row.threshold).to_bits(), threshold.to_bits());
                            calibrated += 1;
                            number(&row.cdf)
                        })
                    } else {
                        None
                    }
                    .unwrap_or_else(|| number(&c.cdf[i * thresholds.len() + j]));
                    close(
                        p.cdf(threshold),
                        expected,
                        2e-8,
                        0.,
                        &format!("{context} cdf({threshold})"),
                        &mut failures,
                    );
                }
            }
        }
        if let Some(reference) = mathematical {
            assert_eq!(mathematical_cases, reference.cases.len());
        }
    }
    if let Some(rows) = calibration {
        assert_eq!(
            calibrated,
            rows.len(),
            "every independently generated CDF must be checked"
        );
    }
    assert!(
        failures.is_empty(),
        "{} parity mismatches:\n{}",
        failures.len(),
        failures
            .iter()
            .take(30)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn uniform_point_masses_boundaries_and_invalid_inputs_remain_explicit() {
    for (w, a, b) in [
        (vec![], vec![], vec![]),
        (vec![1.], vec![], vec![1.]),
        (vec![1.], vec![0.], vec![]),
        (vec![0.], vec![0.], vec![1.]),
        (vec![-1.], vec![0.], vec![1.]),
        (vec![f64::NAN], vec![0.], vec![1.]),
        (vec![1.], vec![2.], vec![1.]),
        (vec![1.], vec![0.], vec![f64::INFINITY]),
        (vec![1.], vec![-1e308], vec![1e308]),
    ] {
        assert!(UniformMixture::new(w, a, b).is_err());
    }
    for df in [f64::NAN, 0.0, -1.0, f64::NEG_INFINITY] {
        assert!(
            UniformMixture::new(vec![1.0], vec![-1.0], vec![1.0])
                .unwrap()
                .with_student_t(df)
                .is_err()
        );
    }
    let prior = UniformMixture::new(vec![2.], vec![1.25], vec![1.25]).unwrap();
    for (x, s) in [
        (f64::NAN, 1.),
        (f64::INFINITY, 1.),
        (1., 0.),
        (1., -1.),
        (1., f64::INFINITY),
    ] {
        assert!(prior.posterior(x, s).is_err());
    }
    let p = prior.posterior(2., 0.5).unwrap();
    assert_eq!(p.mean(), 1.25);
    assert_eq!(p.standard_deviation(), 0.);
    assert_eq!(p.cdf(1.25 - f64::EPSILON), 0.);
    assert_eq!(p.cdf(1.25), 1.);
    assert!(p.cdf(f64::NAN).is_nan());
    let prior = UniformMixture::new(vec![1.], vec![-1.], vec![1.]).unwrap();
    let p = prior.posterior(40., 0.01).unwrap();
    if numeric_profile() == NumericProfile::Portable {
        assert!((p.mean() - 0.999_997_435_897_773).abs() < 1e-12);
        assert!((p.standard_deviation() - 2.564_102_058_362_641e-6).abs() < 1e-12);
        assert_eq!(p.cdf(0.), 0.);

        let estimate = 1e6;
        let se = 0.001_f64;
        let cutoff = 1.0 - 1e-12;
        // The squared-distance correction to the limiting exponential is
        // below f64 precision at one endpoint scale in this case.
        let log_tail = -((estimate - 1.0) / se.powi(2)) * (1.0 - cutoff);
        for (x, threshold, expected) in [
            (estimate, cutoff, log_tail.exp()),
            (-estimate, -cutoff, -log_tail.exp_m1()),
        ] {
            let cdf = prior.cdf_values(x, se, &[threshold]).unwrap()[0];
            assert!((cdf - expected).abs() < 5e-15, "{cdf} != {expected}");
        }
    } else {
        assert!(p.mean().is_nan());
        assert!(p.standard_deviation().is_nan());
        assert!(p.cdf(0.).is_nan());
    }
    let p = prior.posterior(f64::MAX, 0.01).unwrap();
    assert!(p.mean().is_nan());
    assert!(p.standard_deviation().is_nan());
    assert!(p.cdf(0.).is_nan());
}
