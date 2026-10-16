use super::*;

#[test]
fn expanded_formula_design_supports_raw_polynomial_transforms() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0_f64, 1.0, 2.0, 3.0];
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

    let design =
        expanded_formula_design("~ condition + poly(dose, 3, raw=TRUE)", &factors, &numeric)
            .unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "dose_poly_1".to_string(),
            "dose_poly_2".to_string(),
            "dose_poly_3".to_string(),
        ]
    );
    assert_eq!(
        design.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 1.0, 1.0, //
            1.0, 1.0, 2.0, 4.0, 8.0, //
            1.0, 1.0, 3.0, 9.0, 27.0,
        ]
    );
    assert_eq!(
        design.numeric_covariates,
        vec![
            "dose_poly_1".to_string(),
            "dose_poly_2".to_string(),
            "dose_poly_3".to_string(),
        ]
    );

    let interaction = expanded_formula_design(
        "~ condition:poly(dose, degree=2, raw = TRUE)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        interaction.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "conditionA:dose_poly_1".to_string(),
            "conditionB:dose_poly_1".to_string(),
            "conditionA:dose_poly_2".to_string(),
            "conditionB:dose_poly_2".to_string(),
        ]
    );

    let treatment_interaction = expanded_formula_design(
        "~ condition + poly(dose, degree=2, raw = TRUE) + condition:poly(dose, degree=2, raw = TRUE)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert!(
        treatment_interaction
            .standard_design
            .coefficient_index("condition_B_vs_A:dose_poly_1")
            .is_ok()
    );
    assert!(
        treatment_interaction
            .standard_design
            .coefficient_index("condition_B_vs_A:dose_poly_2")
            .is_ok()
    );

    let named_order =
        expanded_formula_design("~ poly(dose, raw = TRUE, degree = 2)", &factors, &numeric)
            .unwrap();
    assert_eq!(
        named_order.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose_poly_1".to_string(),
            "dose_poly_2".to_string(),
        ]
    );

    let named_x =
        expanded_formula_design("~ poly(x=dose, degree=2, raw=TRUE)", &factors, &numeric).unwrap();
    assert_eq!(named_x, named_order);

    let removed = expanded_formula_design(
        "~ condition + poly(dose, 2, raw=TRUE) - poly(dose, 2, raw=TRUE)",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected);
}

#[test]
fn expanded_formula_design_supports_orthogonal_polynomial_transforms() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0_f64, 1.0, 2.0, 3.0];
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

    let design =
        expanded_formula_design("~ condition + poly(dose, 3)", &factors, &numeric).unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "poly(dose, 3)1".to_string(),
            "poly(dose, 3)2".to_string(),
            "poly(dose, 3)3".to_string(),
        ]
    );
    let expected = [
        -0.6708203932499369,
        0.5,
        -0.22360679774997888,
        -0.22360679774997896,
        -0.5,
        0.6708203932499369,
        0.22360679774997896,
        -0.5,
        -0.6708203932499369,
        0.6708203932499369,
        0.5,
        0.22360679774997896,
    ];
    for (row_idx, row) in design
        .standard_design
        .matrix()
        .as_slice()
        .as_chunks::<5>()
        .0
        .iter()
        .enumerate()
    {
        for (observed, expected) in row[2..].iter().zip(&expected[row_idx * 3..row_idx * 3 + 3]) {
            assert!((observed - expected).abs() < 1e-12);
        }
    }

    let raw_false =
        expanded_formula_design("~ poly(dose, degree=2, raw=FALSE)", &factors, &numeric).unwrap();
    let default = expanded_formula_design("~ poly(dose, degree=2)", &factors, &numeric).unwrap();
    assert_eq!(raw_false, default);

    let named_x = expanded_formula_design("~ poly(x=dose, degree=2)", &factors, &numeric).unwrap();
    assert_eq!(named_x, default);
    let reordered_named =
        expanded_formula_design("~ poly(degree=2, x=dose)", &factors, &numeric).unwrap();
    assert_eq!(reordered_named, default);
    let simple =
        expanded_formula_design("~ poly(dose, 2, simple=TRUE)", &factors, &numeric).unwrap();
    assert_eq!(simple, default);

    let interaction = expanded_formula_design(
        "~ condition + poly(dose, 2) + condition:poly(dose, 2)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert!(
        interaction
            .standard_design
            .coefficient_index("condition_B_vs_A:poly(dose, 2)1")
            .is_ok()
    );
    assert!(
        interaction
            .standard_design
            .coefficient_index("condition_B_vs_A:poly(dose, 2)2")
            .is_ok()
    );

    let removed = expanded_formula_design(
        "~ condition + poly(dose, 2) - poly(dose, 2)",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected_removed = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected_removed);

    let repeated_dose = [0.0_f64, 0.0, 1.0, 2.0];
    let repeated_numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &repeated_dose,
    }];
    let repeated = expanded_formula_design("~ poly(dose, 2)", &[], &repeated_numeric).unwrap();
    let repeated_expected = [
        -0.45226701686664544,
        0.21320071635561041,
        -0.45226701686664544,
        0.21320071635561041,
        0.15075567228888181,
        -0.8528028654224417,
        0.753778361444409,
        0.4264014327112208,
    ];
    for (observed, expected) in repeated
        .standard_design
        .matrix()
        .as_slice()
        .as_chunks::<3>()
        .0
        .iter()
        .flat_map(|row| row[1..].iter())
        .zip(repeated_expected)
    {
        assert!((observed - expected).abs() < 1e-12);
    }
}
