use super::*;

#[test]
fn expanded_formula_design_supports_numeric_function_transforms() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [1.0_f64, 4.0, 9.0, 16.0];
    let factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "A",
        levels: None,
    }];
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];

    let design = expanded_formula_design(
        "~ condition + log(dose) + log2(dose) + log10(dose) + log1p(dose) + sqrt(dose) + scale(dose)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "dose_log".to_string(),
            "dose_log2".to_string(),
            "dose_log10".to_string(),
            "dose_log1p".to_string(),
            "dose_sqrt".to_string(),
            "dose_scale".to_string(),
        ]
    );
    let mean = dose.iter().sum::<f64>() / dose.len() as f64;
    let sd = (dose
        .iter()
        .map(|value| {
            let centered = value - mean;
            centered * centered
        })
        .sum::<f64>()
        / (dose.len() as f64 - 1.0))
        .sqrt();
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = design.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[2], value.ln());
        assert_eq!(row[3], value.log2());
        assert_eq!(row[4], value.log10());
        assert_eq!(row[5], value.ln_1p());
        assert_eq!(row[6], value.sqrt());
        assert!((row[7] - ((value - mean) / sd)).abs() < 1e-12);
    }
    assert_eq!(
        design.numeric_covariates,
        vec![
            "dose_log".to_string(),
            "dose_log2".to_string(),
            "dose_log10".to_string(),
            "dose_log1p".to_string(),
            "dose_sqrt".to_string(),
            "dose_scale".to_string(),
        ]
    );

    let interaction =
        expanded_formula_design("~ condition * log2(dose)", &factors, &numeric).unwrap();
    assert!(
        interaction
            .standard_design
            .coefficient_index("condition_B_vs_A:dose_log2")
            .is_ok()
    );

    let spaced = expanded_formula_design(
        "~ log ( dose ) + log2 ( dose ) + sqrt ( x = dose ) + scale ( dose, center = FALSE, scale = FALSE )",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        spaced.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose_log".to_string(),
            "dose_log2".to_string(),
            "dose_sqrt".to_string(),
            "dose_identity".to_string(),
        ]
    );

    let named_function_args = expanded_formula_design(
        "~ log(x=dose) + log(dose, base=2) + log(x=dose, base=10) + log(dose, 4) + log2(x=dose) + sqrt(x=dose)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        named_function_args
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "dose_log".to_string(),
            "dose_log_base_2".to_string(),
            "dose_log_base_10".to_string(),
            "dose_log_base_4".to_string(),
            "dose_log2".to_string(),
            "dose_sqrt".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = named_function_args
            .standard_design
            .matrix()
            .row(sample)
            .unwrap();
        assert_eq!(row[1], value.ln());
        assert!((row[2] - value.log2()).abs() < 1e-12);
        assert!((row[3] - value.log10()).abs() < 1e-12);
        assert!((row[4] - value.log(4.0)).abs() < 1e-12);
        assert_eq!(row[5], value.log2());
        assert_eq!(row[6], value.sqrt());
    }

    let numeric_coercions = expanded_formula_design(
        "~ as.numeric(dose) + as.double(dose) + as.integer(x=dose) + condition:as.numeric(dose)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        numeric_coercions
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "dose_as.numeric".to_string(),
            "dose_as.double".to_string(),
            "dose_as.integer".to_string(),
            "condition_B_vs_A:dose_as.numeric".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = numeric_coercions
            .standard_design
            .matrix()
            .row(sample)
            .unwrap();
        assert_eq!(row[1], value);
        assert_eq!(row[2], value);
        assert_eq!(row[3], value);
        assert_eq!(row[4], if condition[sample] == "B" { value } else { 0.0 });
    }

    let spaced_coercions = expanded_formula_design(
        "~ as.numeric ( dose ) + as.double ( dose ) + as.integer ( x = dose )",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        spaced_coercions
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "dose_as.numeric".to_string(),
            "dose_as.double".to_string(),
            "dose_as.integer".to_string(),
        ]
    );

    let removed =
        expanded_formula_design("~ condition + sqrt(dose) - sqrt(dose)", &factors, &numeric)
            .unwrap();
    let expected = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected);

    let removed_coercion = expanded_formula_design(
        "~ condition + as.numeric(dose) - as.numeric(dose)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(removed_coercion, expected);

    let condition_levels = vec!["B".to_string(), "A".to_string()];
    let aliased_factor = [ExpandedFactorSpec {
        factor: "cell type",
        sample_levels: &condition,
        reference: "A",
        levels: Some(&condition_levels),
    }];
    let factor_coercions = expanded_formula_design(
        "~ as.numeric(cell.type) + as.double(cell.type) + as.integer(`cell type`) + as.numeric(cell.type):dose",
        &aliased_factor,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        factor_coercions
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "cell type_as.numeric".to_string(),
            "cell type_as.double".to_string(),
            "cell type_as.integer".to_string(),
            "cell type_as.numeric:dose".to_string(),
        ]
    );
    for (sample, dose_value) in dose.iter().copied().enumerate() {
        let code = if condition[sample] == "B" { 1.0 } else { 2.0 };
        let row = factor_coercions
            .standard_design
            .matrix()
            .row(sample)
            .unwrap();
        assert_eq!(row[1], code);
        assert_eq!(row[2], code);
        assert_eq!(row[3], code);
        assert_eq!(row[4], code * dose_value);
    }

    let removed_factor_coercion = expanded_formula_design(
        "~ `cell type` + as.numeric(`cell type`) - as.numeric(`cell type`)",
        &aliased_factor,
        &numeric,
    )
    .unwrap();
    let expected_factor =
        expanded_formula_design("~ `cell type`", &aliased_factor, &numeric).unwrap();
    assert_eq!(removed_factor_coercion, expected_factor);

    let duplicate_name_numeric = [
        ExpandedNumericSpec {
            name: "dose",
            values: &dose,
        },
        ExpandedNumericSpec {
            name: "dose_as.numeric",
            values: &dose,
        },
    ];
    let duplicate_name_error =
        expanded_formula_design("~ as.numeric(dose)", &factors, &duplicate_name_numeric)
            .unwrap_err()
            .to_string();
    assert!(
        duplicate_name_error.contains("formula variable 'dose_as.numeric' appears more than once"),
        "{duplicate_name_error}"
    );

    let scaled_interaction =
        expanded_formula_design("~ condition * scale(dose)", &factors, &numeric).unwrap();
    assert!(
        scaled_interaction
            .standard_design
            .coefficient_index("condition_B_vs_A:dose_scale")
            .is_ok()
    );

    let scaled_options = expanded_formula_design(
        "~ scale(dose, center=FALSE) + scale(dose, scale=FALSE) + scale(dose, center=FALSE, scale=FALSE) + scale(dose, FALSE, FALSE) + scale(dose, F, scale=T)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        scaled_options.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose_scale_uncentered".to_string(),
            "dose_center".to_string(),
            "dose_identity".to_string(),
        ]
    );
    let rms =
        (dose.iter().map(|value| value * value).sum::<f64>() / (dose.len() as f64 - 1.0)).sqrt();
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = scaled_options.standard_design.matrix().row(sample).unwrap();
        assert!((row[1] - value / rms).abs() < 1e-12);
        assert!((row[2] - (value - mean)).abs() < 1e-12);
        assert_eq!(row[3], value);
    }

    let positional_scaled_options = expanded_formula_design(
        "~ scale(dose, FALSE, FALSE) + scale(dose, F, T)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        positional_scaled_options
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "dose_identity".to_string(),
            "dose_scale_uncentered".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = positional_scaled_options
            .standard_design
            .matrix()
            .row(sample)
            .unwrap();
        assert_eq!(row[1], value);
        assert!((row[2] - value / rms).abs() < 1e-12);
    }

    let named_x_scaled_options = expanded_formula_design(
        "~ scale(x=dose) + scale(x=dose, center=FALSE, scale=FALSE) + scale(center=FALSE, x=dose, scale=TRUE)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        named_x_scaled_options
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "dose_scale".to_string(),
            "dose_identity".to_string(),
            "dose_scale_uncentered".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = named_x_scaled_options
            .standard_design
            .matrix()
            .row(sample)
            .unwrap();
        assert!((row[1] - ((value - mean) / sd)).abs() < 1e-12);
        assert_eq!(row[2], value);
        assert!((row[3] - value / rms).abs() < 1e-12);
    }

    let scaled_constants = expanded_formula_design(
        "~ scale(dose, center=5, scale=FALSE) + scale(dose, center=FALSE, scale=2) + scale(dose, center=5, scale=2)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        scaled_constants
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "dose_centered".to_string(),
            "dose_scaled".to_string(),
            "dose_centered_scaled".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = scaled_constants
            .standard_design
            .matrix()
            .row(sample)
            .unwrap();
        assert_eq!(row[1], value - 5.0);
        assert_eq!(row[2], value / 2.0);
        assert_eq!(row[3], (value - 5.0) / 2.0);
    }

    let parenthesized_arguments = expanded_formula_design(
        "~ scale(dose, center=(5), scale=FALSE) + scale(dose, center=FALSE, scale=(2))",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        parenthesized_arguments
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "dose_centered".to_string(),
            "dose_scaled".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = parenthesized_arguments
            .standard_design
            .matrix()
            .row(sample)
            .unwrap();
        assert_eq!(row[1], value - 5.0);
        assert_eq!(row[2], value / 2.0);
    }
}

#[test]
fn expanded_formula_design_transform_parser_handles_nested_parentheses_and_quotes() {
    let numeric_with_parens = [1.0_f64, 2.0, 4.0, 8.0];
    let factors = [];
    let numeric = [ExpandedNumericSpec {
        name: "dose(value)",
        values: &numeric_with_parens,
    }];

    let design = expanded_formula_design(
        "~ log2(`dose(value)`) + I(`dose(value)` + 1)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose(value)_log2".to_string(),
            "dose(value)_plus_1".to_string(),
        ]
    );
    assert_eq!(
        resolve_coefficient_index(&design.standard_design, "dose.value._log2").unwrap(),
        1
    );
    assert_eq!(
        resolve_coefficient_index(&design.standard_design, "dose.value._plus_1").unwrap(),
        2
    );
    for (sample, value) in numeric_with_parens.iter().copied().enumerate() {
        let row = design.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[1], value.log2());
        assert_eq!(row[2], value + 1.0);
    }

    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let factor_with_space = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "A",
        levels: None,
    }];
    let factor_design = expanded_formula_design(
        "~ factor (condition) + relevel (condition, ref = 'B')",
        &factor_with_space,
        &[],
    )
    .unwrap();
    assert_eq!(
        factor_design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "factor(condition)_B_vs_A".to_string(),
            "relevel(condition, ref = \"B\")_A_vs_B".to_string(),
        ]
    );

    let literal_log = [1.0_f64, 3.0, 5.0, 7.0];
    let literal_poly = [2.0_f64, 4.0, 6.0, 8.0];
    let literal_identity = [0.5_f64, 1.5, 2.5, 3.5];
    let equals_name = [1.0_f64, 2.0, 3.0, 4.0];
    let literal_numeric = [
        ExpandedNumericSpec {
            name: "log2(dose)",
            values: &literal_log,
        },
        ExpandedNumericSpec {
            name: "poly(dose, 2)",
            values: &literal_poly,
        },
        ExpandedNumericSpec {
            name: "I(dose + 1)",
            values: &literal_identity,
        },
        ExpandedNumericSpec {
            name: "dose=a",
            values: &equals_name,
        },
    ];
    let literal = expanded_formula_design(
        "~ `log2(dose)` + `poly(dose, 2)` + `I(dose + 1)` + log2(`dose=a`) + scale(`dose=a`, FALSE, FALSE)",
        &factors,
        &literal_numeric,
    )
    .unwrap();
    assert_eq!(
        literal.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "log2(dose)".to_string(),
            "poly(dose, 2)".to_string(),
            "I(dose + 1)".to_string(),
            "dose=a_log2".to_string(),
            "dose=a_identity".to_string(),
        ]
    );
    for sample in 0..literal_log.len() {
        let row = literal.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[1], literal_log[sample]);
        assert_eq!(row[2], literal_poly[sample]);
        assert_eq!(row[3], literal_identity[sample]);
        assert_eq!(row[4], equals_name[sample].log2());
        assert_eq!(row[5], equals_name[sample]);
    }
}
