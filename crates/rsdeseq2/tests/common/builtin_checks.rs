//! Shared independent numerical gates for continuous and count likelihood fixtures.
use crate::{location_shrinkage, support::*};
pub fn kernels(reference: serde_json::Value) {
    let mut failures = vec![];
    for c in reference["kernel_cases"].as_array().unwrap() {
        let input: builtin_likelihood::Case = serde_json::from_value(c.clone()).unwrap();
        let result = builtin_likelihood::run(&input);
        if c["expected_error"].is_string() {
            assert!(result.is_err(), "{}", input.name);
            continue;
        }
        let r = result.unwrap();
        for (i, name) in ["log_cdf", "log_density", "first", "second"]
            .iter()
            .enumerate()
        {
            compare(
                &r[i],
                &c[name],
                if i < 2 { 2e-8 } else { 2e-6 },
                2e-12,
                &format!("{} {name}", input.name),
                &mut failures,
            );
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
pub fn complete_fits(reference: serde_json::Value) {
    let mut failures = vec![];
    for c in reference["cases"].as_array().unwrap() {
        let input: location_shrinkage::Case = serde_json::from_value(c.clone()).unwrap();
        let r = location_shrinkage::run(&input);
        if c["expected_error"].is_string() {
            assert!(r.is_err(), "{}", input.name);
            continue;
        }
        let r = r.unwrap_or_else(|e| panic!("{}: {e}", input.name));
        for (name, v) in [
            ("posterior_mean", &r.posterior_mean),
            ("posterior_sd", &r.posterior_sd),
            ("svalues", &r.svalues),
            ("cdf", &r.cdf),
        ] {
            compare(
                v,
                &c[name],
                2e-6,
                0.0,
                &format!("{} {name}", input.name),
                &mut failures,
            );
        }
        compare(
            &[r.log_likelihood],
            &serde_json::json!([c["log_likelihood"]]),
            2e-6,
            0.0,
            &format!("{} likelihood", input.name),
            &mut failures,
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
