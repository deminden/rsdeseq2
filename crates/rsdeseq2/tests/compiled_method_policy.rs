//! Original compiled-prefit versus refinement control and zero-row contracts.
use crate::support::regression_checks as checks;
#[test]
fn compiled_method_controls_preserve_original_posteriors_errors_and_rng_state() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/compiled_method_policy_r462.json")).unwrap();
    for c in reference["cases"].as_array().unwrap() {
        let calls = c["original_refinement_calls"].as_u64().unwrap();
        let name = c["name"].as_str().unwrap();
        if name.starts_with("valid_prefit") || name.contains("map") {
            assert_eq!(calls, 0, "{name}");
        }
        if name.starts_with("refinement_") {
            assert!(calls > 0, "{name}");
        }
    }
    checks::compare_posteriors(reference, true);
}
