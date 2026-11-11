//! Original SQP control profiles and independent direction/validation contracts.
#[path = "common/posterior.rs"]
mod support;
use support::*;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/em_fit.rs"]
mod em_fit;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/location_shrinkage.rs"]
mod location_shrinkage;
use rsdeseq2::{
    matrix::{MatrixOrder, MatrixView},
    shrinkage::*,
};
#[test]
fn sqp_control_profiles_match_original_numerical_outputs() {
    let r: serde_json::Value =
        serde_json::from_str(include_str!("data/sqp_controls_r462.json")).unwrap();
    let mut failures = vec![];
    for c in r["cases"].as_array().unwrap() {
        let input: location_shrinkage::Case = serde_json::from_value(c.clone()).unwrap();
        let f = location_shrinkage::run(&input);
        if c["expected_error"].is_string() {
            assert!(f.is_err(), "{}", input.name);
            continue;
        }
        let f = match f {
            Ok(f) => f,
            Err(e) => {
                failures.push(format!("{}: {e}", input.name));
                continue;
            }
        };
        let target = c;
        for (field, value) in [
            ("posterior_mean", &f.posterior_mean),
            ("posterior_sd", &f.posterior_sd),
            ("svalues", &f.svalues),
        ] {
            compare(
                value,
                &target[field],
                2e-6,
                0.0,
                &format!("{} {field}", input.name),
                &mut failures,
            );
        }
        if let Some(exact) = c.get("deterministic_reference") {
            for (field, value) in [
                ("posterior_mean", &f.posterior_mean),
                ("posterior_sd", &f.posterior_sd),
                ("svalues", &f.svalues),
            ] {
                compare(
                    value,
                    &exact[field],
                    2e-6,
                    0.0,
                    &format!("{} independent exact-SVD {field}", input.name),
                    &mut failures,
                );
            }
        }
        let ll = parse_numeric_input(target["log_likelihood"].as_str().unwrap()).unwrap();
        if (f.log_likelihood - ll).abs() > 2e-6 {
            failures.push(format!(
                "{} loglik {} vs {}",
                input.name, f.log_likelihood, ll
            ));
        }
    }
    for c in r["em_fit_cases"].as_array().unwrap() {
        let input: em_fit::Case = serde_json::from_value(c.clone()).unwrap();
        let f = em_fit::run(&input);
        if c["expected_error"].is_string() {
            assert!(f.is_err(), "{}", input.name);
            continue;
        }
        let f = match f {
            Ok(f) => f,
            Err(e) => {
                failures.push(format!("{} primitive: {e}", input.name));
                continue;
            }
        };
        if input.name == "wide_prior_truncated" {
            println!(
                "RAW-NATIVE {:?} SUM {} RES {}",
                f.weights, f.optimizer_weight_sum, f.dual_residual
            );
        }
        let target = c;
        compare(
            &f.weights,
            &target["expected_weights"],
            2e-6,
            0.0,
            &format!("{} weights", input.name),
            &mut failures,
        );
        for (name, value) in [
            ("expected_log_likelihood", f.log_likelihood),
            (
                "expected_penalized_log_likelihood",
                f.penalized_log_likelihood,
            ),
        ] {
            let expected = parse_numeric_input(target[name].as_str().unwrap()).unwrap();
            if (value - expected).abs() > 2e-6 {
                failures.push(format!("{} {name}: {value} vs {expected}", input.name));
            }
        }
        if f.converged != c["converged"].as_bool().unwrap() {
            failures.push(format!(
                "{} convergence {} vs {}",
                input.name, f.converged, c["converged"]
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
#[test]
fn direction_cutoff_has_a_known_em_state_and_invalid_controls_reject() {
    let l = [1.0_f64, 0.3, 0.2, 1.0, 0.7, 0.4, 0.5, 1.0].map(f64::ln);
    let matrix = MatrixView::new(4, 2, &l, MatrixOrder::RowMajor).unwrap();
    let options = MixtureFitOptions {
        max_iterations: 1,
        em_iterations: 0,
        gradient_tolerance: 0.0,
        sqp: MixtureSqpOptions {
            search_direction_threshold: Some(2.0),
            ..Default::default()
        },
        ..Default::default()
    };
    let f = fit_mixture_proportions(matrix, None, &[1.0, 1.0], None, options).unwrap();
    let expected = (1.0 / 1.3 + 0.2 / 1.2 + 0.7 / 1.1 + 0.5 / 1.5) / 4.0;
    assert!((f.weights[0] - expected).abs() < 1e-14);
    let uncut = fit_mixture_proportions(
        matrix,
        None,
        &[1.0, 1.0],
        None,
        MixtureFitOptions {
            sqp: Default::default(),
            ..options
        },
    )
    .unwrap();
    assert!((uncut.weights[0] - f.weights[0]).abs() > 1e-5);
    for sqp in [
        MixtureSqpOptions {
            active_tolerance: -1.0,
            ..Default::default()
        },
        MixtureSqpOptions {
            step_reduction: 1.0,
            ..Default::default()
        },
        MixtureSqpOptions {
            minimum_step: 0.0,
            ..Default::default()
        },
        MixtureSqpOptions {
            svd_tolerance: Some(f64::NAN),
            ..Default::default()
        },
    ] {
        assert!(
            fit_mixture_proportions(
                matrix,
                None,
                &[1.0, 1.0],
                None,
                MixtureFitOptions {
                    sqp,
                    ..Default::default()
                }
            )
            .is_err()
        );
    }
}

#[test]
fn nearly_unit_backtracking_factors_terminate_with_visible_nonconvergence() {
    let l = [1.0_f64, 0.3, 0.2, 1.0, 0.7, 0.4, 0.5, 1.0].map(f64::ln);
    let fit = fit_mixture_proportions(
        MatrixView::new(4, 2, &l, MatrixOrder::RowMajor).unwrap(),
        None,
        &[1.0, 1.0],
        None,
        MixtureFitOptions {
            max_iterations: 1,
            em_iterations: 0,
            gradient_tolerance: 0.0,
            sqp: MixtureSqpOptions {
                sufficient_decrease: 0.99,
                step_reduction: 1.0 - f64::EPSILON,
                ..Default::default()
            },
            ..Default::default()
        },
    )
    .unwrap();
    assert!(!fit.converged);
    assert!(fit.message.contains("line search iteration limit"));
    assert!(fit.function_evaluations < 20_000);
    assert!((fit.weights[0] - (1.0 / 1.3 + 0.2 / 1.2 + 0.7 / 1.1 + 0.5 / 1.5) / 4.0).abs() < 1e-14);
}
