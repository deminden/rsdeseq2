use rsdeseq2::shrinkage::{
    LocationLikelihood, LocationShrinkageOptions, MixtureComponents, MixtureGrid,
    MixtureGridOptions, MixturePriorPenalty, shrink_location_batch,
};
use serde::Deserialize;
#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/formula.rs"]
mod formula;
#[path = "../tools/validation/evaluate_r_interface/likelihood.rs"]
mod likelihood;
#[derive(Deserialize)]
struct Reference {
    cases: Vec<Case>,
}
#[derive(Deserialize)]
struct Prior {
    component_class: String,
    weights: Vec<String>,
    means: Option<Vec<String>>,
    sds: Option<Vec<String>>,
    lower: Option<Vec<String>>,
    upper: Option<Vec<String>>,
}
#[derive(Deserialize)]
struct Case {
    name: String,
    estimates: Vec<String>,
    standard_errors: Vec<String>,
    thresholds: Vec<String>,
    family: String,
    point_mass: bool,
    nullweight: Option<String>,
    scales: Option<Vec<String>>,
    mode: String,
    range: Vec<String>,
    alpha: String,
    prune_threshold: String,
    observation_weights: Option<Vec<String>>,
    fixed_prior: bool,
    supplied_prior: Option<Prior>,
    df: Option<String>,
    custom_likelihood: Option<serde_json::Value>,
    expected_error: Option<String>,
    #[serde(default)]
    posterior_mean: Vec<String>,
    #[serde(default)]
    posterior_sd: Vec<String>,
    #[serde(default)]
    negative_probability: Vec<String>,
    #[serde(default)]
    positive_probability: Vec<String>,
    #[serde(default)]
    local_false_discovery_rate: Vec<String>,
    #[serde(default)]
    local_false_sign_rate: Vec<String>,
    #[serde(default)]
    svalues: Vec<String>,
    #[serde(default)]
    qvalues: Vec<String>,
    #[serde(default)]
    cdf: Vec<String>,
    log_likelihood: Option<String>,
    warnings: Vec<String>,
}
fn number(s: &str) -> f64 {
    if let Some(s) = s.strip_prefix("bits:") {
        return f64::from_bits(u64::from_str_radix(s, 16).unwrap());
    }
    match s {
        "NA" => f64::from_bits(0x7ff0_0000_0000_07a2),
        "NaN" => f64::NAN,
        "Inf" => f64::INFINITY,
        "-Inf" => f64::NEG_INFINITY,
        _ => s.parse().unwrap(),
    }
}
fn numbers(v: &[String]) -> Vec<f64> {
    v.iter().map(|s| number(s)).collect()
}
fn is_na(x: f64) -> bool {
    x.is_nan() && x.to_bits() as u32 == 1954
}
fn compare(
    actual: &[f64],
    expected: &[String],
    bound: f64,
    context: &str,
    failures: &mut Vec<String>,
) {
    if actual.len() != expected.len() {
        failures.push(format!(
            "{context}: lengths {} != {}",
            actual.len(),
            expected.len()
        ));
        return;
    }
    for (i, (&a, b)) in actual.iter().zip(expected).enumerate() {
        let b = number(b);
        if a.is_nan() && b.is_nan() && is_na(a) == is_na(b) {
            continue;
        }
        if a == b {
            continue;
        }
        if a.is_finite() && b.is_finite() && (a - b).abs() <= bound {
            continue;
        }
        failures.push(format!(
            "{context}[{i}] {a:.17e} != {b:.17e}; NA {} != {}",
            is_na(a),
            is_na(b)
        ));
    }
}
#[test]
fn complete_location_batches_match_original_fit_and_summary_contract() {
    let mut reference: Reference =
        serde_json::from_str(include_str!("data/location_shrinkage_r462.json")).unwrap();
    let custom: Reference =
        serde_json::from_str(include_str!("data/custom_likelihood_r462.json")).unwrap();
    reference.cases.extend(custom.cases);
    let mut failures = vec![];
    assert_eq!(reference.cases.len(), 23);
    for c in reference.cases {
        let penalty = c
            .nullweight
            .as_ref()
            .map(|v| MixturePriorPenalty::NullBiased(number(v)))
            .unwrap_or(MixturePriorPenalty::Uniform);
        let supplied_grid = c.supplied_prior.as_ref().map(|p| {
            let weights = numbers(&p.weights);
            let mut dirichlet = vec![1.0; weights.len()];
            if let Some(w) = &c.nullweight {
                dirichlet[0] = number(w);
            }
            let components = match p.component_class.as_str() {
                "normalmix" => MixtureComponents::Normal {
                    means: numbers(p.means.as_ref().unwrap()),
                    sds: numbers(p.sds.as_ref().unwrap()),
                },
                "unimix" => MixtureComponents::Uniform {
                    lower: numbers(p.lower.as_ref().unwrap()),
                    upper: numbers(p.upper.as_ref().unwrap()),
                },
                "tnormalmix" => MixtureComponents::TruncatedNormal {
                    means: numbers(p.means.as_ref().unwrap()),
                    sds: numbers(p.sds.as_ref().unwrap()),
                    lower: numbers(p.lower.as_ref().unwrap()),
                    upper: numbers(p.upper.as_ref().unwrap()),
                },
                _ => panic!("unexpected prior"),
            };
            MixtureGrid {
                components,
                initial_weights: weights,
                dirichlet,
                warnings: vec![],
            }
        });
        let options = LocationShrinkageOptions {
            grid: MixtureGridOptions {
                family: c.family.parse().unwrap(),
                scales: c.scales.as_ref().map(|v| numbers(v)),
                point_mass: c.point_mass,
                mode: number(&c.mode),
                range: (number(&c.range[0]), number(&c.range[1])),
                penalty,
                ..Default::default()
            },
            supplied_grid,
            fixed_prior: c.fixed_prior,
            likelihood: if let Some(value) = &c.custom_likelihood {
                LocationLikelihood::Restricted(std::sync::Arc::new(
                    likelihood::restricted(value).unwrap(),
                ))
            } else {
                c.df.as_ref()
                    .map(|df| LocationLikelihood::StudentT(number(df)))
                    .unwrap_or(LocationLikelihood::Normal)
            },
            alpha: number(&c.alpha),
            prune_threshold: number(&c.prune_threshold),
            observation_weights: c.observation_weights.as_ref().map(|v| numbers(v)),
            thresholds: numbers(&c.thresholds),
            ..Default::default()
        };
        let result = shrink_location_batch(
            &numbers(&c.estimates),
            &numbers(&c.standard_errors),
            &options,
        );
        if c.expected_error.is_some() {
            if result.is_ok() {
                failures.push(format!("{}: original rejects this batch", c.name));
            }
            continue;
        }
        let r = match result {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{} unexpected error: {e}", c.name));
                continue;
            }
        };
        // The original R matrix is frozen in cutoff-major storage.
        let reference_cdf = (0..c.estimates.len())
            .flat_map(|i| (0..c.thresholds.len()).map(move |j| (i, j)))
            .map(|(i, j)| c.cdf[j * c.estimates.len() + i].clone())
            .collect::<Vec<_>>();
        for (field, actual, expected) in [
            ("mean", r.posterior_mean, &c.posterior_mean),
            ("sd", r.posterior_sd, &c.posterior_sd),
            ("np", r.negative_probability, &c.negative_probability),
            ("pp", r.positive_probability, &c.positive_probability),
            (
                "lfdr",
                r.local_false_discovery_rate,
                &c.local_false_discovery_rate,
            ),
            ("lfsr", r.local_false_sign_rate, &c.local_false_sign_rate),
            ("svalue", r.svalues, &c.svalues),
            ("qvalue", r.qvalues, &c.qvalues),
            ("cdf", r.cdf, &reference_cdf),
        ] {
            compare(
                &actual,
                expected,
                2e-6,
                &format!("{} {field}", c.name),
                &mut failures,
            );
        }
        compare(
            &[r.log_likelihood],
            &[c.log_likelihood.unwrap()],
            2e-6,
            &format!("{} loglik", c.name),
            &mut failures,
        );
        if r.warnings != c.warnings {
            failures.push(format!(
                "{} warnings {:?} != {:?}",
                c.name, r.warnings, c.warnings
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} batch mismatches:\n{}",
        failures.len(),
        failures
            .iter()
            .take(60)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}
