#[macro_use]
#[path = "support/formula_reference.rs"]
mod reference;
#[path = "../examples/evaluate_r_interface/factor_expansion.rs"]
mod factor_expansion;
#[allow(dead_code)]
#[path = "../examples/evaluate_r_interface/formula.rs"]
mod formula;
#[path = "../examples/evaluate_r_interface/native_formula.rs"]
mod input;

#[test]
fn factor_classes_and_matrix_column_maps_match_original_r() {
    let reference: serde_json::Value =
        serde_json::from_str(formula_reference!("native_factor_r462")).unwrap();
    let input: input::Input = serde_json::from_value(reference.clone()).unwrap();
    let actual = input::evaluate(&input);
    let cases = reference["native_formula_cases"].as_array().unwrap();
    for (case, actual) in cases.iter().zip(actual["cases"].as_array().unwrap()) {
        assert_eq!(
            actual["ok"],
            !case["expected"].is_null(),
            "{} decision",
            case["name"]
        );
        if !case["expected"].is_null() {
            for (field, expected) in [
                ("names", "native_names"),
                ("values", "expected"),
                ("expanded_names", "expanded_names"),
                ("expanded_values", "expanded_expected"),
                ("groups", "groups"),
                ("frames", "frames"),
            ] {
                assert_eq!(actual[field], case[expected], "{} {field}", case["name"]);
            }
        }
    }
    for (case, actual) in reference["factor_expansions"]
        .as_array()
        .unwrap()
        .iter()
        .zip(actual["expansions"].as_array().unwrap())
    {
        assert_eq!(
            actual["ok"],
            !case["expected"].is_null(),
            "{} decision",
            case["name"]
        );
        if !case["expected"].is_null() {
            assert_eq!(
                actual["value"], case["expected"],
                "{} primitive",
                case["name"]
            );
        }
    }
}
