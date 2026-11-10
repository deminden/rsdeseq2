//! Independent original-ashr optimizer states and complete posterior batches.
#[path = "common/posterior.rs"]
mod support;
use support::*;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/em_fit.rs"]
mod em_fit;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/location_shrinkage.rs"]
mod location_shrinkage;
fn reference() -> serde_json::Value {
    serde_json::from_str(include_str!("data/ashr_optimizer_r462.json")).unwrap()
}
#[test]
fn accelerated_em_and_variational_states_match_original_ashr() {
    let mut failures = vec![];
    let mut stale_source_objectives = vec![];
    for c in reference()["em_fit_cases"].as_array().unwrap() {
        let input: em_fit::Case = serde_json::from_value(c.clone()).unwrap();
        let f = match em_fit::run(&input) {
            Ok(f) => f,
            Err(e) => {
                failures.push(format!("{}: {e}", input.name));
                continue;
            }
        };
        compare(
            &f.weights,
            c.get("corrected_weights")
                .filter(|v| !v.is_null())
                .unwrap_or(&c["expected_weights"]),
            2e-6,
            0.0,
            &input.name,
            &mut failures,
        );
        assert_eq!(
            f.converged,
            c["converged"].as_bool().unwrap(),
            "{} convergence",
            input.name
        );
        let expected_iterations = c["iterations"].as_u64();
        assert_eq!(
            f.iterations.filter(|v| *v != 0),
            expected_iterations.map(|v| v as usize),
            "{} iterations",
            input.name
        );
        if let Some(objective) = f.optimizer_objective {
            if input.optimizer == "cxxMixSquarem"
                && (parse_numeric_input(c["optimizer_objective"].as_str().unwrap()).unwrap()
                    - parse_numeric_input(c["fresh_original_objective"].as_str().unwrap()).unwrap())
                .abs()
                    > 2e-6
            {
                stale_source_objectives.push(input.name.clone());
            }
            compare(
                &[objective],
                &serde_json::json!([c["native_objective_reference"]]),
                2e-6,
                0.0,
                &format!("{} objective", input.name),
                &mut failures,
            );
        }
        if let Some(post) = &f.posterior_concentrations {
            compare(
                post,
                &c["posterior_concentrations"],
                2e-6,
                0.0,
                &format!("{} concentrations", input.name),
                &mut failures,
            );
        }
    }
    assert_eq!(
        stale_source_objectives,
        [
            "budget_two_cxxMixSquarem",
            "budget_three_cxxMixSquarem",
            "budget_four_cxxMixSquarem",
            "budget_ten_cxxMixSquarem",
            "early_second_cxxMixSquarem",
            "categorical_compiled_stop"
        ]
    );
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
#[test]
fn alternative_optimizers_keep_complete_posterior_and_warning_coverage() {
    let mut failures = vec![];
    for c in reference()["cases"].as_array().unwrap() {
        let input: location_shrinkage::Case = serde_json::from_value(c.clone()).unwrap();
        let f = match location_shrinkage::run(&input) {
            Ok(f) => f,
            Err(e) => {
                failures.push(format!("{}: {e}", input.name));
                continue;
            }
        };
        for (field, actual) in [
            ("posterior_mean", &f.posterior_mean),
            ("posterior_sd", &f.posterior_sd),
            ("svalues", &f.svalues),
            ("local_false_sign_rate", &f.local_false_sign_rate),
        ] {
            compare(
                actual,
                &c[field],
                2e-6,
                0.0,
                &format!("{} {field}", input.name),
                &mut failures,
            );
        }
        let expected = c["fitted_prior"]["weights"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| {
                serde_json::Value::String(format!(
                    "{:.17e}",
                    parse_numeric_input(v.as_str().unwrap()).unwrap()
                ))
            })
            .collect::<Vec<_>>();
        compare(
            &f.fitted_weights,
            &serde_json::Value::Array(expected),
            2e-6,
            0.0,
            &format!("{} weights", input.name),
            &mut failures,
        );
        compare(
            &[f.log_likelihood],
            &serde_json::json!([c["log_likelihood"]]),
            2e-6,
            0.0,
            &format!("{} likelihood", input.name),
            &mut failures,
        );
        // The original compiled optimizer reads an empty diagnostic vector.
        // Safe native code does not reproduce an out-of-bounds access/warning.
        let expected = c["warnings"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap().to_owned())
            .filter(|v| v != "subscript out of bounds (index 0 >= vector size 0)")
            .collect::<Vec<_>>();
        assert_eq!(f.warnings, expected, "{} warnings", input.name);
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn em_validation_does_not_silently_change_optimizer_or_prior_semantics() {
    let r = reference();
    let base = &r["em_fit_cases"][0];
    let mut variants = vec![];
    let mut c = base.clone();
    c["observation_weights"] = serde_json::json!(["1"]);
    variants.push((c, "weighted EM"));
    let mut c = base.clone();
    c["control"] = serde_json::json!({"K":2});
    variants.push((c, "higher-order"));
    let mut c = base.clone();
    c["control"] = serde_json::json!({"tol":-1});
    variants.push((c, "controls"));
    let mut c = base.clone();
    c["log_likelihoods"][0] = serde_json::json!("NaN");
    variants.push((c, "log likelihoods"));
    let mut c = base.clone();
    c["dirichlet"][0] = serde_json::json!("0.5");
    variants.push((c, "Dirichlet"));
    for (c, message) in variants {
        let c: em_fit::Case = serde_json::from_value(c).unwrap();
        assert!(
            em_fit::run(&c).unwrap_err().to_string().contains(message),
            "{}: {message}",
            c.name
        );
    }
    // No fixed-point calls are allowed when the requested budget is zero.
    let mut c = base.clone();
    c["control"] = serde_json::json!({"maxiter":0});
    let f = em_fit::run(&serde_json::from_value(c).unwrap()).unwrap();
    assert!(!f.converged);
    assert_eq!(f.function_evaluations, 0);
    let initial = base["initial_weights"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| parse_numeric_input(v.as_str().unwrap()).unwrap())
        .collect::<Vec<_>>();
    let sum = initial.iter().sum::<f64>();
    assert!(
        f.weights
            .iter()
            .zip(initial)
            .all(|(a, b)| (a - b / sum).abs() < 1e-15)
    );
}

#[test]
fn em_methods_have_independent_categorical_and_dirichlet_solutions() {
    use rsdeseq2::{
        matrix::{MatrixOrder, MatrixView},
        shrinkage::{
            MixtureFitOptions, MixtureOptimizer as M, SquaremOptions, fit_mixture_proportions,
        },
    };
    let n = f64::NEG_INFINITY;
    let ll = [0.0, n, 0.0, n, 0.0, n, n, 0.0];
    let o = SquaremOptions::default();
    for (method, weights, prior, expected) in [
        (M::Em(o), None, [2.0, 3.0], [4.0 / 7.0, 3.0 / 7.0]),
        (M::CxxEm(o), None, [2.0, 3.0], [4.0 / 7.0, 3.0 / 7.0]),
        (
            M::WeightedEm(o),
            Some([2.0, 0.0, 1.0, 4.0]),
            [2.0, 3.0],
            [0.4, 0.6],
        ),
        (
            M::VariationalEm(o),
            None,
            [0.2, 0.4],
            [3.2 / 4.6, 1.4 / 4.6],
        ),
    ] {
        let f = fit_mixture_proportions(
            MatrixView::new(4, 2, &ll, MatrixOrder::RowMajor).unwrap(),
            weights.as_ref().map(|v| &v[..]),
            &prior,
            None,
            MixtureFitOptions {
                optimizer: method,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(f.converged);
        for (a, b) in f.weights.iter().zip(expected) {
            assert!(
                (a - b).abs() < 1e-14,
                "{method:?}: {:?} expected {expected:?}",
                f.weights
            );
        }
        if method.is_variational() {
            let post = f.posterior_concentrations.unwrap();
            assert!((post[0] - 3.2).abs() < 1e-14);
            assert!((post[1] - 1.4).abs() < 1e-14);
        }
    }
}
