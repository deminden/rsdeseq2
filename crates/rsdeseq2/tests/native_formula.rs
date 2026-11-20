use crate::support::{formula_reference, matches_within_ulps};
use rsdeseq2::{NumericProfile, numeric_profile};
#[path = "../tools/validation/evaluate_r_interface/factor_expansion.rs"]
mod factor_expansion;
#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/formula.rs"]
mod formula;
#[path = "../tools/validation/evaluate_r_interface/native_formula.rs"]
mod input;

#[test]
fn actual_native_formula_parser_matches_independent_r_matrices_and_offsets() {
    let reference: serde_json::Value =
        serde_json::from_str(formula_reference!("native_formula_r462")).unwrap();
    let input: input::Input = serde_json::from_value(reference.clone()).unwrap();
    let actual = input::evaluate(&input);
    let cases = reference["native_formula_cases"].as_array().unwrap();
    assert_eq!(cases.len(), 54);
    let portable = numeric_profile() == NumericProfile::Portable;
    for (case, actual) in cases.iter().zip(actual["cases"].as_array().unwrap()) {
        let valid = !case["expected"].is_null();
        assert_eq!(actual["name"], case["name"]);
        assert_eq!(actual["ok"], valid, "{}", case["name"]);
        if valid {
            assert_eq!(
                actual["names"], case["native_names"],
                "{} labels",
                case["name"]
            );
            let name = case["name"].as_str().unwrap();
            let ulps = match (portable, name) {
                (true, "p_log_base_3") => 2,
                (true, "power_7" | "power_8" | "p_log" | "z_log1p") => 1,
                _ => 0,
            };
            if ulps == 0 {
                assert_eq!(actual["values"], case["expected"], "{name} values");
            } else {
                let actual = actual["values"].as_array().unwrap();
                let expected = case["expected"].as_array().unwrap();
                assert_eq!(actual.len(), expected.len(), "{name} values length");
                let value = |hex: &serde_json::Value| {
                    f64::from_bits(u64::from_str_radix(hex.as_str().unwrap(), 16).unwrap())
                };
                for (i, (a, b)) in actual.iter().zip(expected).enumerate() {
                    assert!(
                        matches_within_ulps(value(a), value(b), ulps),
                        "{name} element {i}: {a} vs {b}"
                    );
                }
            }
            assert_eq!(
                actual["offsets"], case["offsets"],
                "{} offsets",
                case["name"]
            );
        }
    }
}
