use super::*;

#[test]
fn expanded_formula_design_supports_numeric_power_transforms() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0, 1.0, 2.0, 3.0];
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
        expanded_formula_design("~ condition + dose + I(dose^2)", &factors, &numeric).unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "dose".to_string(),
            "dose_pow_2".to_string(),
        ]
    );
    assert_eq!(
        design.standard_design.matrix().as_slice(),
        &[
            1.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 1.0, //
            1.0, 1.0, 2.0, 4.0, //
            1.0, 1.0, 3.0, 9.0,
        ]
    );
    assert_eq!(
        design.numeric_covariates,
        vec!["dose".to_string(), "dose_pow_2".to_string()]
    );

    let interaction =
        expanded_formula_design("~ condition * I(dose^2)", &factors, &numeric).unwrap();
    assert!(
        interaction
            .standard_design
            .coefficient_index("condition_B_vs_A:dose_pow_2")
            .is_ok()
    );

    let removed =
        expanded_formula_design("~ condition + I(dose^2) - I(dose^2)", &factors, &numeric).unwrap();
    let expected = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected);
}

#[test]
fn expanded_formula_design_supports_numeric_identity_transform() {
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

    let design = expanded_formula_design("~ condition + I(dose)", &factors, &numeric).unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "dose_identity".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        assert_eq!(
            design.standard_design.matrix().row(sample).unwrap()[2],
            value
        );
    }
    assert_eq!(design.numeric_covariates, vec!["dose_identity".to_string()]);

    let interaction = expanded_formula_design("~ condition * I(dose)", &factors, &numeric).unwrap();
    assert!(
        interaction
            .standard_design
            .coefficient_index("condition_B_vs_A:dose_identity")
            .is_ok()
    );

    let removed =
        expanded_formula_design("~ condition + I(dose) - I(dose)", &factors, &numeric).unwrap();
    let expected = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected);
}

#[test]
fn expanded_formula_design_supports_signed_numeric_identity_transform() {
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

    let design = expanded_formula_design("~ I(+dose) + I(-dose)", &factors, &numeric).unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose_identity".to_string(),
            "dose_neg".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = design.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[1], value);
        assert_eq!(row[2], -value);
    }

    let interaction =
        expanded_formula_design("~ condition * I(-dose)", &factors, &numeric).unwrap();
    assert!(
        interaction
            .standard_design
            .coefficient_index("condition_B_vs_A:dose_neg")
            .is_ok()
    );

    let removed =
        expanded_formula_design("~ condition + I(-dose) - I(-dose)", &factors, &numeric).unwrap();
    let expected = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected);
}

#[test]
fn expanded_formula_design_supports_numeric_scalar_arithmetic_transform() {
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

    let design = expanded_formula_design(
        "~ I(dose + 1) + I(dose - 1) + I(dose * 2) + I(dose / 2)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose_plus_1".to_string(),
            "dose_minus_1".to_string(),
            "dose_times_2".to_string(),
            "dose_div_2".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = design.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[1], value + 1.0);
        assert_eq!(row[2], value - 1.0);
        assert_eq!(row[3], value * 2.0);
        assert_eq!(row[4], value / 2.0);
    }

    let interaction =
        expanded_formula_design("~ condition * I(dose + 1)", &factors, &numeric).unwrap();
    assert!(
        interaction
            .standard_design
            .coefficient_index("condition_B_vs_A:dose_plus_1")
            .is_ok()
    );

    let removed = expanded_formula_design(
        "~ condition + I(dose + 1) - I(dose + 1)",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected);
}

#[test]
fn expanded_formula_design_supports_scalar_left_numeric_arithmetic_transform() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [1.0_f64, 2.0, 4.0, 8.0];
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
        "~ I(1 + dose) + I(1 - dose) + I(2 * dose) + I(2 / dose)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose_plus_1".to_string(),
            "dose_rminus_1".to_string(),
            "dose_times_2".to_string(),
            "dose_rdiv_2".to_string(),
        ]
    );
    for (sample, value) in dose.iter().copied().enumerate() {
        let row = design.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[1], 1.0 + value);
        assert_eq!(row[2], 1.0 - value);
        assert_eq!(row[3], 2.0 * value);
        assert_eq!(row[4], 2.0 / value);
    }

    let interaction =
        expanded_formula_design("~ condition * I(1 + dose)", &factors, &numeric).unwrap();
    assert!(
        interaction
            .standard_design
            .coefficient_index("condition_B_vs_A:dose_plus_1")
            .is_ok()
    );
}

#[test]
fn expanded_formula_design_supports_numeric_binary_arithmetic_transform() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [1.0_f64, 2.0, 4.0, 8.0];
    let time = [1.0_f64, 2.0, 3.0, 4.0];
    let dose_over_time = [1.0_f64, 1.0, 4.0 / 3.0, 2.0];
    let factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "A",
        levels: None,
    }];
    let numeric = [
        ExpandedNumericSpec {
            name: "dose",
            values: &dose,
        },
        ExpandedNumericSpec {
            name: "time",
            values: &time,
        },
        ExpandedNumericSpec {
            name: "dose/time",
            values: &dose_over_time,
        },
    ];

    let design = expanded_formula_design(
        "~ I(dose + time) + I(dose - time) + I(dose * time) + I(dose / time)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose_plus_time".to_string(),
            "dose_minus_time".to_string(),
            "dose_times_time".to_string(),
            "dose_div_time".to_string(),
        ]
    );
    for sample in 0..dose.len() {
        let row = design.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[1], dose[sample] + time[sample]);
        assert_eq!(row[2], dose[sample] - time[sample]);
        assert_eq!(row[3], dose[sample] * time[sample]);
        assert_eq!(row[4], dose[sample] / time[sample]);
    }

    let interaction =
        expanded_formula_design("~ condition * I(dose + time)", &factors, &numeric).unwrap();
    assert!(
        interaction
            .standard_design
            .coefficient_index("condition_B_vs_A:dose_plus_time")
            .is_ok()
    );

    let named_x = expanded_formula_design(
        "~ I(x=dose) + I(x=dose + 1) + I(x=`dose/time`)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        named_x.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "dose_identity".to_string(),
            "dose_plus_1".to_string(),
            "dose/time_identity".to_string(),
        ]
    );
    for (sample, dose_value) in dose.iter().copied().enumerate() {
        let row = named_x.standard_design.matrix().row(sample).unwrap();
        assert_eq!(row[1], dose_value);
        assert_eq!(row[2], dose_value + 1.0);
        assert_eq!(row[3], numeric[2].values[sample]);
    }

    let removed = expanded_formula_design(
        "~ condition + I(dose + time) - I(dose + time)",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(removed, expected);
}
