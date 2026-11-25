//! Original ashr fits whose mixsqp subproblems have numerically singular
//! reduced Hessians; a wrong search direction ends at a worse vertex.
use crate::support::posterior as support;
use rsdeseq2::shrinkage::*;
use support::*;

#[test]
fn singular_subproblems_reach_the_original_vertex() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/ashr_singular_direction_r462.json")).unwrap();
    let cases = reference["cases"].as_array().unwrap();
    assert_eq!(cases.len(), 4);
    let strings = |v: &serde_json::Value| {
        v.as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_owned())
            .collect::<Vec<_>>()
    };
    let mut failures = vec![];
    for c in cases {
        let name = c["name"].as_str().unwrap();
        let estimates = parse_numeric_inputs(&strings(&c["estimates"])).unwrap();
        let errors = parse_numeric_inputs(&strings(&c["standard_errors"])).unwrap();
        let query = AdaptiveShrinkageQuery {
            estimates: &estimates,
            standard_errors: &errors,
            lfc_threshold: 0.0,
            svalues: false,
        };
        let result = adaptive_shrinkage(
            &estimates,
            &errors,
            &LocationShrinkageOptions::default(),
            query,
        )
        .unwrap();
        assert!(result.fit.warnings.is_empty(), "{name}");
        assert!(c["warnings"].as_array().unwrap().is_empty(), "{name}");
        let fit = result.fit.fit.as_ref().unwrap();
        assert!(fit.converged, "{name}");
        let expected = parse_numeric_input(c["log_likelihood"].as_str().unwrap()).unwrap();
        if (fit.log_likelihood - expected).abs() > 2e-6 {
            failures.push(format!(
                "{name} loglik {} vs {expected}",
                fit.log_likelihood
            ));
        }
        for (field, value) in [
            ("fitted_weights", &result.fit.fitted_weights),
            ("posterior_mean", &result.posterior_mean),
            ("posterior_sd", &result.posterior_sd),
        ] {
            compare(
                value,
                &c[field],
                2e-6,
                0.0,
                &format!("{name} {field}"),
                &mut failures,
            );
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
