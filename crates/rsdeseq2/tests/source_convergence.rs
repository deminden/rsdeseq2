use rsdeseq2::shrinkage::{LocationShrinkageOptions, MixtureComponents, shrink_location_batch};
use serde::Deserialize;
#[derive(Deserialize)]
struct Reference {
    estimates: Vec<String>,
    standard_errors: Vec<String>,
    fitted_weights: Vec<String>,
    fitted_sds: Vec<String>,
    posterior_mean: Vec<String>,
    posterior_sd: Vec<String>,
}
fn numbers(v: &[String]) -> Vec<f64> {
    v.iter()
        .map(|s| {
            if let Some(s) = s.strip_prefix("bits:") {
                f64::from_bits(u64::from_str_radix(s, 16).unwrap())
            } else {
                s.parse().unwrap()
            }
        })
        .collect()
}
#[test]
fn automatic_prior_retains_original_positive_coordinate_convergence() {
    let r: Reference =
        serde_json::from_str(include_str!("data/mixture_convergence_r462.json")).unwrap();
    let fit = shrink_location_batch(
        &numbers(&r.estimates),
        &numbers(&r.standard_errors),
        &LocationShrinkageOptions::default(),
    )
    .unwrap();
    assert!(fit.fit.as_ref().unwrap().converged);
    assert_eq!(fit.fitted_weights, numbers(&r.fitted_weights));
    let MixtureComponents::Normal { sds, .. } = &fit.fitted_grid.components else {
        panic!("normal grid required")
    };
    for (actual, expected, bound) in [
        (sds.as_slice(), &r.fitted_sds, 2e-14),
        (fit.posterior_mean.as_slice(), &r.posterior_mean, 2e-6),
        (fit.posterior_sd.as_slice(), &r.posterior_sd, 2e-6),
    ] {
        assert_eq!(actual.len(), expected.len());
        for (&a, b) in actual.iter().zip(numbers(expected)) {
            assert!((a - b).abs() <= bound, "{a:.17e} != {b:.17e}");
        }
    }
}
