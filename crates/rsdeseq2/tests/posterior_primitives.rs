use rsdeseq2::shrinkage::{NormalMixture, svalues};
use serde::Deserialize;

#[derive(Deserialize)]
struct References {
    cases: Vec<Case>,
    svalues: Vec<SvalueCase>,
}
#[derive(Deserialize)]
struct Case {
    name: String,
    weights: Vec<String>,
    means: Vec<String>,
    sds: Vec<String>,
    estimates: Vec<String>,
    standard_errors: Vec<String>,
    thresholds: Vec<String>,
    posterior_weights: Vec<String>,
    posterior_mean: Vec<String>,
    posterior_sd: Vec<String>,
    cdf: Vec<String>,
    upper_tail: Vec<String>,
}
#[derive(Deserialize)]
struct SvalueCase {
    name: String,
    rates: Vec<String>,
    values: Vec<String>,
}
fn number(s: &str) -> f64 {
    match s {
        "Inf" => f64::INFINITY,
        "-Inf" => f64::NEG_INFINITY,
        "NA" | "NaN" => f64::NAN,
        _ => s.parse().unwrap(),
    }
}
fn numbers(x: &[String]) -> Vec<f64> {
    x.iter().map(|s| number(s)).collect()
}
fn references() -> References {
    serde_json::from_str(include_str!("data/normal_mixture_posterior_r462.json")).unwrap()
}
fn close(a: f64, b: f64, absolute: f64, relative: f64, context: &str) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= absolute + relative * b.abs(),
        "{context}: {a:.17e} != {b:.17e}"
    );
}
#[test]
fn normal_mixture_moments_weights_and_tails_match_independent_ashr() {
    let references = references();
    assert_eq!(references.cases.len(), 8);
    for case in references.cases {
        let prior = NormalMixture::new(
            numbers(&case.weights),
            numbers(&case.means),
            numbers(&case.sds),
        )
        .unwrap();
        let thresholds = numbers(&case.thresholds);
        for (i, (&estimate, se)) in numbers(&case.estimates)
            .iter()
            .zip(numbers(&case.standard_errors))
            .enumerate()
        {
            let posterior = prior.posterior(estimate, se).unwrap();
            close(
                posterior.mean(),
                number(&case.posterior_mean[i]),
                2e-12,
                2e-12,
                &format!("{}[{i}] mean", case.name),
            );
            close(
                posterior.standard_deviation(),
                number(&case.posterior_sd[i]),
                2e-12,
                2e-12,
                &format!("{}[{i}] SD", case.name),
            );
            for (j, &w) in posterior.weights().iter().enumerate() {
                close(
                    w,
                    number(&case.posterior_weights[i * case.weights.len() + j]),
                    1e-300,
                    2e-11,
                    &format!("{}[{i},{j}] weight", case.name),
                );
            }
            for (j, &c) in thresholds.iter().enumerate() {
                close(
                    posterior.cdf(c),
                    number(&case.cdf[i * thresholds.len() + j]),
                    1e-300,
                    2e-11,
                    &format!("{}[{i},{c}] CDF", case.name),
                );
                close(
                    posterior.survival(c),
                    number(&case.upper_tail[i * thresholds.len() + j]),
                    1e-300,
                    2e-11,
                    &format!("{}[{i},{c}] upper tail", case.name),
                );
            }
        }
    }
}
#[test]
fn ordered_svalues_match_independent_apeglm_including_ties_and_missing() {
    let references = references();
    assert_eq!(references.svalues.len(), 8);
    for case in references.svalues {
        let rates = numbers(&case.rates)
            .into_iter()
            .map(|p| p.is_finite().then_some(p))
            .collect::<Vec<_>>();
        let actual = svalues(&rates).unwrap();
        assert_eq!(actual.len(), case.values.len());
        for (a, b) in actual.into_iter().zip(numbers(&case.values)) {
            if b.is_nan() {
                assert_eq!(a, None, "{}", case.name);
            } else {
                close(a.unwrap(), b, 2e-15, 0.0, &case.name);
            }
        }
    }
}
#[test]
fn malformed_posterior_inputs_fail_and_point_masses_are_right_continuous() {
    for (w, m, s) in [
        (vec![], vec![], vec![]),
        (vec![1.], vec![], vec![1.]),
        (vec![1.], vec![0.], vec![]),
        (vec![0.], vec![0.], vec![1.]),
        (vec![-1.], vec![0.], vec![1.]),
        (vec![f64::NAN], vec![0.], vec![1.]),
        (vec![f64::INFINITY], vec![0.], vec![1.]),
        (vec![1.], vec![f64::NAN], vec![1.]),
        (vec![1.], vec![0.], vec![-1.]),
        (vec![1.], vec![0.], vec![f64::INFINITY]),
        (vec![1.], vec![0.], vec![1e200]),
    ] {
        assert!(NormalMixture::new(w, m, s).is_err());
    }
    let prior = NormalMixture::new(vec![2.], vec![1.25], vec![0.]).unwrap();
    for (b, s) in [
        (f64::NAN, 1.),
        (f64::INFINITY, 1.),
        (1., 0.),
        (1., -1.),
        (1., f64::INFINITY),
        (1., f64::NAN),
        (1., 1e-200),
        (1., 1e200),
    ] {
        assert!(prior.posterior(b, s).is_err());
    }
    let p = prior.posterior(2., 0.5).unwrap();
    assert_eq!(p.mean(), 1.25);
    assert_eq!(p.standard_deviation(), 0.);
    assert_eq!(p.cdf(1.25 - f64::EPSILON), 0.);
    assert_eq!(p.cdf(1.25), 1.);
    assert!(p.cdf(f64::NAN).is_nan());
    for p in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(svalues(&[Some(p)]).is_err());
    }
}
