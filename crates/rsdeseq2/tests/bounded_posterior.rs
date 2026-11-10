use rsdeseq2::shrinkage::{FixedMixture, TruncatedNormalMixture, UniformMixture};
use serde::Deserialize;
#[derive(Deserialize)]
struct References {
    cases: Vec<Case>,
}
#[derive(Deserialize)]
struct Case {
    name: String,
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
fn number(v: &str) -> f64 {
    if let Some(bits) = v.strip_prefix("bits:") {
        return f64::from_bits(u64::from_str_radix(bits, 16).unwrap());
    }
    match v {
        "NaN" | "NA" => f64::NAN,
        "Inf" => f64::INFINITY,
        "-Inf" => f64::NEG_INFINITY,
        _ => v.parse().unwrap(),
    }
}
fn numbers(v: &[String]) -> Vec<f64> {
    v.iter().map(|v| number(v)).collect()
}
fn close(a: f64, b: f64, absolute: f64, relative: f64, context: &str) {
    if b.is_nan() {
        assert!(a.is_nan(), "{context}");
        return;
    }
    if a == b {
        return;
    }
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= absolute + relative * b.abs(),
        "{context}: {a:.17e} != {b:.17e}"
    );
}
#[test]
fn bounded_mixture_likelihoods_posteriors_and_cdfs_match_independent_ashr() {
    for fixture in [
        include_str!("data/uniform_posterior_r462.json"),
        include_str!("data/truncated_posterior_r462.json"),
    ] {
        let r: References = serde_json::from_str(fixture).unwrap();
        for c in r.cases {
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
                FixedMixture::Uniform(
                    UniformMixture::new(numbers(&c.weights), numbers(&c.lower), numbers(&c.upper))
                        .unwrap(),
                )
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
                let expected_mean = number(&c.posterior_mean[i]);
                let expected_sd = number(&c.posterior_sd[i]);
                let likelihood = prior.log_likelihoods(estimate, se).unwrap();
                for (j, &a) in likelihood.iter().enumerate() {
                    close(
                        a,
                        number(&c.log_likelihoods[i * c.weights.len() + j]),
                        1e-8,
                        1e-12,
                        &format!("{context} likelihood{j}"),
                    );
                }
                let p = &posteriors[i];
                close(
                    p.mean(),
                    expected_mean,
                    2e-6,
                    0.,
                    &format!("{context} mean"),
                );
                close(
                    p.standard_deviation(),
                    expected_sd,
                    2e-6,
                    0.,
                    &format!("{context} SD"),
                );
                for (j, &a) in p.weights().iter().enumerate() {
                    close(
                        a,
                        number(&c.posterior_weights[i * c.weights.len() + j]),
                        2e-8,
                        0.,
                        &format!("{context} weight{j}"),
                    );
                }
                for (j, &threshold) in thresholds.iter().enumerate() {
                    close(
                        p.cdf(threshold),
                        number(&c.cdf[i * thresholds.len() + j]),
                        2e-8,
                        0.,
                        &format!("{context} cdf({threshold})"),
                    );
                }
            }
        }
    }
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
    assert!(p.mean().is_nan());
    assert!(p.standard_deviation().is_nan());
    assert!(p.cdf(0.).is_nan());
}
