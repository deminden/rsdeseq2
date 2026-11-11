//! Independent ignored/deferred controls and scalar-bound policies across methods.
#[path = "common/regression_checks.rs"]
mod checks;
#[path = "common/posterior.rs"]
mod support;

#[test]
fn regression_controls_match_original_accepted_outputs_and_rejection_categories() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/regression_method_policy_r462.json")).unwrap();
    for case in reference["cases"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        if name.starts_with("betabinC") {
            let calls = case["original_beta_optimizer_calls"].as_u64().unwrap();
            if name.contains("refinement") {
                assert!(calls > 0, "{name}");
            } else {
                assert_eq!(calls, 0, "{name}");
            }
        }
    }
    checks::compare_posteriors(reference, true);
}
