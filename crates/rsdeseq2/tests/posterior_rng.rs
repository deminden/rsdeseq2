//! Independent original-R seed, uniform, normal and state-transition references.
use rsdeseq2::shrinkage::PosteriorRng;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/random.rs"]
mod random;
fn bits(v: f64) -> String {
    format!("{:016x}", v.to_bits())
}
fn ordered(v: f64) -> u64 {
    let x = v.to_bits();
    if x >> 63 == 0 { x | 1 << 63 } else { !x }
}
fn number(v: &serde_json::Value) -> f64 {
    f64::from_bits(u64::from_str_radix(v.as_str().unwrap(), 16).unwrap())
}
#[test]
fn native_rng_matches_original_r_uniforms_normals_and_serialized_states() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/posterior_rng_r462.json")).unwrap();
    let mut failures = vec![];
    for c in reference["random_cases"].as_array().unwrap() {
        let input: random::Case = serde_json::from_value(c.clone()).unwrap();
        let original = random::generator(&input).unwrap();
        assert_eq!(
            random::encode_state(&original.r_seed()),
            input.state,
            "{} initial state",
            input.name
        );
        let mut u = original.clone();
        let mut n = original;
        for (i, (ue, ne)) in c["uniforms"]
            .as_array()
            .unwrap()
            .iter()
            .zip(c["normals"].as_array().unwrap())
            .enumerate()
        {
            let actual = u.uniform();
            let expected = number(ue);
            if actual.to_bits() != expected.to_bits() {
                failures.push(format!(
                    "{} uniform {i}: {actual:.17e} != {expected:.17e}",
                    input.name
                ));
            }
            let actual = n.start_coordinate();
            let expected = number(ne);
            let ulp = ordered(actual).abs_diff(ordered(expected));
            if !actual.is_finite() || ulp > 16 {
                failures.push(format!(
                    "{} normal {i}: {actual:.17e} != {expected:.17e}; {ulp} ULP",
                    input.name
                ));
            }
        }
        let encode = |v: &serde_json::Value| {
            v.as_array()
                .unwrap()
                .iter()
                .map(|s| s.as_str().unwrap().to_owned())
                .collect::<Vec<_>>()
        };
        assert_eq!(
            random::encode_state(&u.r_seed()),
            encode(&c["uniform_state"]),
            "{} uniform state",
            input.name
        );
        assert_eq!(
            random::encode_state(&n.r_seed()),
            encode(&c["normal_state"]),
            "{} normal state",
            input.name
        );
    }
    assert!(
        failures.is_empty(),
        "{} mismatches:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
#[test]
fn unsupported_or_malformed_rng_states_fail_without_draws() {
    assert!(PosteriorRng::from_seed(i32::MIN).is_err());
    for state in [
        vec![],
        vec![10403, 624],
        vec![10407, 1, 2],
        vec![10203, 624],
        vec![10407, 0, 0, 0, 1, 2, 3],
        vec![10407, -1, 1, 2, 1, 2, 3],
    ] {
        assert!(PosteriorRng::from_r_seed(&state).is_err(), "{state:?}");
    }
    let mut state = PosteriorRng::from_seed(42).unwrap().r_seed();
    state[2..].fill(0);
    assert!(PosteriorRng::from_r_seed(&state).is_err());
    let rng = PosteriorRng::from_entropy();
    assert!(!rng.has_drawn());
    assert!(PosteriorRng::from_r_seed(&rng.r_seed()).is_ok());
}

#[test]
fn failed_fit_retains_drawn_state_but_early_validation_does_not_consume_it() {
    use rsdeseq2::{
        matrix::{MatrixOrder, MatrixView},
        shrinkage::{
            NegativeBinomialPosteriorInput, NegativeBinomialPosteriorOptions, PosteriorOffsets,
            fit_negative_binomial_posterior_with_rng,
        },
    };
    let design = [1., 0., 1., 1.];
    let extreme = [f64::MAX, f64::MAX];
    let options = NegativeBinomialPosteriorOptions {
        coefficient: Some(1),
        ..Default::default()
    };
    let input = NegativeBinomialPosteriorInput {
        counts: MatrixView::new(1, 2, &extreme, MatrixOrder::RowMajor).unwrap(),
        design: MatrixView::new(2, 2, &design, MatrixOrder::RowMajor).unwrap(),
        dispersions: &[0.1],
        offsets: PosteriorOffsets::None,
        weights: None,
        mle: None,
    };
    let original = PosteriorRng::from_seed(42).unwrap();
    let mut expected = original.clone();
    expected.start_vector(2);
    let mut rng = original.clone();
    assert!(fit_negative_binomial_posterior_with_rng(input, &options, &mut rng).is_err());
    assert!(rng.has_drawn());
    assert_eq!(rng.r_seed(), expected.r_seed());
    let invalid = [f64::NAN, 1.];
    let mut rng = original.clone();
    assert!(
        fit_negative_binomial_posterior_with_rng(
            NegativeBinomialPosteriorInput {
                counts: MatrixView::new(1, 2, &invalid, MatrixOrder::RowMajor).unwrap(),
                ..input
            },
            &options,
            &mut rng
        )
        .is_err()
    );
    assert!(!rng.has_drawn());
    assert_eq!(rng.r_seed(), original.r_seed());
}
