//! Independent cycling references and certified interior-point optima.
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
    serde_json::from_str(include_str!("data/ashr_advanced_optimizer_r462.json")).unwrap()
}
#[test]
fn cycling_and_interior_point_queries_match_independent_references() {
    let r = reference();
    let mut failures = vec![];
    for c in r["cases"].as_array().unwrap() {
        let input: location_shrinkage::Case = serde_json::from_value(c.clone()).unwrap();
        let result = location_shrinkage::run(&input);
        if c["expected_error"].is_string() {
            assert!(
                result.unwrap_err().to_string().contains("numeric method"),
                "{}",
                input.name
            );
            continue;
        }
        let f = match result {
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
        let weights = c["fitted_prior"]["weights"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| {
                serde_json::Value::String(format!(
                    "{:.17e}",
                    parse_numeric_input(v.as_str().unwrap()).unwrap()
                ))
            })
            .collect();
        compare(
            &f.fitted_weights,
            &serde_json::Value::Array(weights),
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
        let expected_warnings = c["warnings"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap().to_owned())
            .filter(|v| v != "subscript out of bounds (index 0 >= vector size 0)")
            .collect::<Vec<_>>();
        assert_eq!(f.warnings, expected_warnings, "{} warnings", input.name);
        assert_eq!(
            f.fit.as_ref().unwrap().converged,
            c["converged"].as_bool().unwrap(),
            "{} convergence",
            input.name
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
#[test]
fn primitive_cycling_states_keep_source_budgets_and_variational_counts() {
    let r = reference();
    let mut failures = vec![];
    for c in r["em_fit_cases"].as_array().unwrap() {
        let input: em_fit::Case = serde_json::from_value(c.clone()).unwrap();
        let f = em_fit::run(&input);
        if c["expected_error"].is_string() {
            assert!(f.unwrap_err().to_string().contains("numeric method"));
            continue;
        }
        let f = f.unwrap();
        compare(
            &f.weights,
            &c["expected_weights"],
            2e-6,
            0.0,
            &input.name,
            &mut failures,
        );
        if input.optimizer != "mixIP" {
            assert_eq!(
                f.converged,
                c["converged"].as_bool().unwrap(),
                "{} convergence",
                input.name
            );
            assert_eq!(
                f.iterations,
                c["iterations"].as_u64().map(|v| v as usize),
                "{} iterations",
                input.name
            );
            if let Some(p) = f.posterior_concentrations {
                compare(
                    &p,
                    &c["posterior_concentrations"],
                    2e-6,
                    0.0,
                    &format!("{} concentrations", input.name),
                    &mut failures,
                );
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
#[test]
fn interior_point_has_independent_high_precision_kkt_certificates() {
    let truth: serde_json::Value =
        serde_json::from_str(include_str!("data/interior_mixture_truth.json")).unwrap();
    let reference = reference();
    let mut failures = vec![];
    for t in truth["cases"].as_array().unwrap() {
        let c = reference["em_fit_cases"]
            .as_array()
            .unwrap()
            .iter()
            .find(|v| v["name"] == t["name"])
            .unwrap();
        let f = em_fit::run(&serde_json::from_value(c.clone()).unwrap()).unwrap();
        assert!(f.converged);
        assert!(f.message.starts_with("interior-point"));
        assert!(f.dual_residual < 1e-8);
        compare(
            &f.weights,
            &t["weights"],
            2e-7,
            0.0,
            &format!("{} weights", c["name"]),
            &mut failures,
        );
        compare(
            &[f.log_likelihood],
            &serde_json::json!([t["log_likelihood"]]),
            2e-7,
            0.0,
            &format!("{} likelihood", c["name"]),
            &mut failures,
        );
        compare(
            &[f.penalized_log_likelihood],
            &serde_json::json!([t["penalized_log_likelihood"]]),
            2e-7,
            0.0,
            &format!("{} penalized", c["name"]),
            &mut failures,
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn native_barrier_obeys_simplex_solutions_limits_and_validation() {
    use rsdeseq2::{
        matrix::{MatrixOrder, MatrixView},
        shrinkage::{
            InteriorPointOptions, MixtureFitOptions, MixtureOptimizer, fit_mixture_proportions,
        },
    };
    let n = f64::NEG_INFINITY;
    let ll = [0.0, n, 0.0, n, 0.0, n, n, 0.0];
    let fit = |o, weights: Option<&[f64]>, prior: &[f64]| {
        fit_mixture_proportions(
            MatrixView::new(4, 2, &ll, MatrixOrder::RowMajor).unwrap(),
            weights,
            prior,
            None,
            MixtureFitOptions {
                optimizer: MixtureOptimizer::InteriorPoint(o),
                ..Default::default()
            },
        )
    };
    for (weights, prior, target) in [
        (None, [2.0, 3.0], [4.0 / 7.0, 3.0 / 7.0]),
        (Some(&[2.0, 0.0, 1.0, 4.0][..]), [2.0, 3.0], [0.4, 0.6]),
        (Some(&[0.0; 4][..]), [2.0, 4.0], [0.25, 0.75]),
    ] {
        let f = fit(InteriorPointOptions::default(), weights, &prior).unwrap();
        assert!(f.converged);
        assert!(
            f.weights
                .iter()
                .zip(target)
                .all(|(a, b)| (a - b).abs() < 1e-8)
        );
        assert!(f.dual_residual < 1e-8);
    }
    let f = fit(
        InteriorPointOptions {
            max_iterations: 0,
            ..Default::default()
        },
        None,
        &[1.0, 1.0],
    )
    .unwrap();
    assert!(!f.converged);
    assert_eq!(f.iterations, Some(0));
    assert!(
        fit(
            InteriorPointOptions {
                relative_tolerance: 0.0,
                ..Default::default()
            },
            None,
            &[1.0, 1.0]
        )
        .unwrap_err()
        .to_string()
        .contains("tolerance")
    );
    // Priors below one select variational fitting, never a convex barrier fit.
    assert!(fit(InteriorPointOptions::default(), None, &[0.2, 0.4]).is_err());
}
