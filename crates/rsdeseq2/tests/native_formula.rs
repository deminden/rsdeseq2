#[macro_use]
#[path = "support/formula_reference.rs"]
mod reference;
#[path = "../tools/validation/evaluate_r_interface/factor_expansion.rs"]
mod factor_expansion;
#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/formula.rs"]
mod formula;
#[path = "../tools/validation/evaluate_r_interface/native_formula.rs"]
mod input;

#[test]
fn actual_native_formula_parser_matches_independent_r_matrix_and_offset_bits() {
    let reference: serde_json::Value =
        serde_json::from_str(formula_reference!("native_formula_r462")).unwrap();
    let input: input::Input = serde_json::from_value(reference.clone()).unwrap();
    let actual = input::evaluate(&input);
    let cases = reference["native_formula_cases"].as_array().unwrap();
    assert_eq!(cases.len(), 54);
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
            assert_eq!(
                actual["values"], case["expected"],
                "{} values",
                case["name"]
            );
            assert_eq!(
                actual["offsets"], case["offsets"],
                "{} offsets",
                case["name"]
            );
        }
    }
}
