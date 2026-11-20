//! Independent original-ashr optimizer states and complete posterior batches.
use crate::support::json_fixture;
use crate::support::posterior as support;
use rsdeseq2::DeseqError;
use rsdeseq2::errors::{CalculationLocation, ErrorKind};
use std::assert_matches;
use support::*;
#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/location_shrinkage.rs"]
mod location_shrinkage;
#[test]
fn accelerated_em_and_variational_states_match_original_ashr() {
    let mut failures = vec![];
    let mut stale_source_objectives = vec![];
    for c in json_fixture!("ashr_optimizer_r462.json")["em_fit_cases"]
        .as_array()
        .unwrap()
    {
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
    for c in json_fixture!("ashr_optimizer_r462.json")["cases"]
        .as_array()
        .unwrap()
    {
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
    let r = json_fixture!("ashr_optimizer_r462.json");
    let base = &r["em_fit_cases"][0];
    let reject = |path: &str, value: serde_json::Value| {
        let mut c = base.clone();
        *c.pointer_mut(path).unwrap() = value;
        *em_fit::run(&serde_json::from_value(c).unwrap())
            .unwrap_err()
            .downcast::<DeseqError>()
            .unwrap()
    };
    assert_matches!(
        reject("/observation_weights", serde_json::json!(["1"])),
        DeseqError::WeightedEmRequired { weights: 1 }
    );
    assert_matches!(
        reject("/control", serde_json::json!({"K":5})),
        DeseqError::EmPolynomialOrderExceeded {
            order: 5,
            components: 4
        }
    );
    assert_matches!(
        reject("/control", serde_json::json!({"K":-1})),
        DeseqError::InvalidOptimizerCount { control, value: -1.0 } if control == "K"
    );
    let invalid_tolerance = reject("/control", serde_json::json!({"tol":-1}));
    assert_eq!(
        invalid_tolerance.to_string(),
        "invalid options: mixture prior fit: invalid squared-EM controls",
    );
    assert_matches!(
        invalid_tolerance,
        DeseqError::InvalidCalculation {
            domain: ErrorKind::Options,
            operation: "squared-EM tolerance",
            location: None,
            value: -1.0,
            ..
        }
    );
    assert_matches!(
        reject("/log_likelihoods/6", serde_json::json!("NaN")),
        DeseqError::InvalidCalculation {
            domain: ErrorKind::Options, operation: "EM log likelihood",
            location: Some(CalculationLocation::MixtureLikelihood { observation: 1, component: 2 }),
            value, ..
        } if value.is_nan()
    );
    assert_matches!(
        reject("/dirichlet/1", serde_json::json!("0.5")),
        DeseqError::InvalidCalculation {
            domain: ErrorKind::Options,
            operation: "EM Dirichlet parameter",
            location: Some(CalculationLocation::Component(1)),
            value: 0.5,
            ..
        }
    );
    assert_matches!(
        reject("/dirichlet", serde_json::json!(["1"])),
        DeseqError::MixtureInputLengthMismatch {
            input: "EM Dirichlet parameters",
            expected: 4,
            actual: 1,
            ..
        }
    );
    let mut c = base.clone();
    c["optimizer"] = serde_json::json!("w_mixEM");
    let mut weights = vec!["1"; 48];
    weights[5] = "-0.5";
    c["observation_weights"] = serde_json::json!(weights);
    let error = *em_fit::run(&serde_json::from_value(c).unwrap())
        .unwrap_err()
        .downcast::<DeseqError>()
        .unwrap();
    assert_matches!(
        error,
        DeseqError::InvalidCalculation {
            domain: ErrorKind::Options,
            operation: "EM observation weight",
            location: Some(CalculationLocation::Observation(5)),
            value: -0.5,
            ..
        }
    );
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
    let row_major = [0.0, n, n, 0.0, n, n, 0.0, n, n, n, n, 0.0];
    let column_major = [0.0, 0.0, 0.0, n, n, n, n, n, n, n, n, 0.0];
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
        let prior = [prior[0], 1.0, prior[1]];
        let expected = [expected[0], 0.0, expected[1]];
        for (order, values) in [
            (MatrixOrder::RowMajor, &row_major),
            (MatrixOrder::ColumnMajor, &column_major),
        ] {
            let f = fit_mixture_proportions(
                MatrixView::new(4, 3, values, order).unwrap(),
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
                    "{method:?} {order:?}: {:?} expected {expected:?}",
                    f.weights
                );
            }
            if method.is_variational() {
                let post = f.posterior_concentrations.unwrap();
                assert!((post[0] - 3.2).abs() < 1e-14);
                assert_eq!(post[1], 0.0);
                assert!((post[2] - 1.4).abs() < 1e-14);
            }
        }
    }
}
