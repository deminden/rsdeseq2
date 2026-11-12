use super::*;

#[test]
fn expanded_formula_design_supports_dot_main_effects() {
    let condition = vec![
        "B".to_string(),
        "A".to_string(),
        "B".to_string(),
        "A".to_string(),
    ];
    let batch = vec![
        "X".to_string(),
        "X".to_string(),
        "Y".to_string(),
        "Y".to_string(),
    ];
    let dose = [1.0, 2.0, 4.0, 8.0];
    let factors = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
            levels: None,
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
            levels: None,
        },
    ];
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];

    let dot = expanded_formula_design("~ .", &factors, &numeric).unwrap();
    let explicit =
        expanded_formula_design("~ condition + batch + dose", &factors, &numeric).unwrap();
    assert_eq!(dot, explicit);

    let removed = expanded_formula_design("~ . - dose - condition", &factors, &numeric).unwrap();
    let expected_removed = expanded_formula_design("~ batch", &factors, &numeric).unwrap();
    assert_eq!(removed, expected_removed);

    let star_dot = expanded_formula_design("~ condition * .", &factors, &numeric).unwrap();
    let star_explicit = expanded_formula_design(
        "~ condition + batch + dose + condition:batch + condition:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(star_dot, star_explicit);

    let interaction_dot = expanded_formula_design("~ condition:.", &factors, &numeric).unwrap();
    let interaction_explicit =
        expanded_formula_design("~ condition:batch + condition:dose", &factors, &numeric).unwrap();
    assert_eq!(interaction_dot, interaction_explicit);

    let nested_dot = expanded_formula_design("~ condition / .", &factors, &numeric).unwrap();
    let nested_explicit = expanded_formula_design(
        "~ condition + condition:batch + condition:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(nested_dot, nested_explicit);

    assert!(expanded_formula_design("~ .", &[], &[]).is_err());
}

#[test]
fn expanded_formula_design_resolves_r_cleaned_variable_aliases() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [1.0, 2.0, 4.0, 8.0];
    let factors = [ExpandedFactorSpec {
        factor: "cell type",
        sample_levels: &condition,
        reference: "A",
        levels: None,
    }];
    let numeric = [ExpandedNumericSpec {
        name: "dose value",
        values: &dose,
    }];

    let cleaned = expanded_formula_design(
        "~ cell.type + dose.value + cell.type:dose.value",
        &factors,
        &numeric,
    )
    .unwrap();
    let quoted = expanded_formula_design(
        "~ `cell type` + `dose value` + `cell type`:`dose value`",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(cleaned, quoted);

    let cleaned_removed =
        expanded_formula_design("~ cell.type + dose.value - dose.value", &factors, &numeric)
            .unwrap();
    let factor_only = expanded_formula_design("~ `cell type`", &factors, &numeric).unwrap();
    assert_eq!(cleaned_removed, factor_only);

    let cleaned_transform = expanded_formula_design(
        "~ factor(cell.type) + relevel(cell.type, ref='B') + droplevels(cell.type)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        cleaned_transform
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "factor(cell type)_B_vs_A".to_string(),
            "relevel(cell type, ref = \"B\")_A_vs_B".to_string(),
            "droplevels(cell type)_B_vs_A".to_string(),
        ]
    );

    let cleaned_numeric_transforms = expanded_formula_design(
        "~ log2(dose.value) + I(dose.value + 1) + scale(dose.value, center=FALSE, scale=FALSE) + poly(dose.value, 2, raw=TRUE)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(
        cleaned_numeric_transforms
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "dose value_log2".to_string(),
            "dose value_plus_1".to_string(),
            "dose value_identity".to_string(),
            "dose value_poly_1".to_string(),
            "dose value_poly_2".to_string(),
        ]
    );
}

#[test]
fn expanded_formula_design_rejects_ambiguous_r_cleaned_variable_aliases() {
    let first = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let second = vec![
        "X".to_string(),
        "Y".to_string(),
        "X".to_string(),
        "Y".to_string(),
    ];
    let factors = [
        ExpandedFactorSpec {
            factor: "cell type",
            sample_levels: &first,
            reference: "A",
            levels: None,
        },
        ExpandedFactorSpec {
            factor: "cell-type",
            sample_levels: &second,
            reference: "X",
            levels: None,
        },
    ];

    let error = expanded_formula_design("~ cell.type", &factors, &[])
        .unwrap_err()
        .to_string();
    assert!(error.contains("resolves ambiguously"));
    let transform_error = expanded_formula_design("~ factor(cell.type)", &factors, &[])
        .unwrap_err()
        .to_string();
    assert!(transform_error.contains("resolves ambiguously"));

    let transform_alias_levels = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let transform_alias_collision = [
        ExpandedFactorSpec {
            factor: "cell type",
            sample_levels: &transform_alias_levels,
            reference: "A",
            levels: None,
        },
        ExpandedFactorSpec {
            factor: "factor.cell.type.",
            sample_levels: &transform_alias_levels,
            reference: "A",
            levels: None,
        },
    ];
    let derived_factor_alias_error =
        expanded_formula_design("~ factor(`cell type`)", &transform_alias_collision, &[])
            .unwrap_err()
            .to_string();
    assert!(derived_factor_alias_error.contains("resolves ambiguously"));

    let dose = [0.0, 1.0];
    let exact_cross_type_error = expanded_formula_design(
        "~ condition",
        &[ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &first,
            reference: "A",
            levels: None,
        }],
        &[ExpandedNumericSpec {
            name: "condition",
            values: &dose,
        }],
    )
    .unwrap_err()
    .to_string();
    assert!(exact_cross_type_error.contains("formula variable 'condition' appears more than once"));

    let alias_cross_type_error = expanded_formula_design(
        "~ cell.type",
        &[ExpandedFactorSpec {
            factor: "cell type",
            sample_levels: &first,
            reference: "A",
            levels: None,
        }],
        &[ExpandedNumericSpec {
            name: "cell-type",
            values: &dose,
        }],
    )
    .unwrap_err()
    .to_string();
    assert!(alias_cross_type_error.contains("resolves ambiguously"));

    let numeric = [
        ExpandedNumericSpec {
            name: "dose value",
            values: &dose,
        },
        ExpandedNumericSpec {
            name: "dose-value",
            values: &dose,
        },
    ];
    let numeric_error = expanded_formula_design("~ log2(dose.value)", &[], &numeric)
        .unwrap_err()
        .to_string();
    assert!(numeric_error.contains("resolves ambiguously"));
    let numeric_identity_error = expanded_formula_design("~ I(dose.value + 1)", &[], &numeric)
        .unwrap_err()
        .to_string();
    assert!(numeric_identity_error.contains("resolves ambiguously"));

    let numeric_transform_alias_collision = [
        ExpandedNumericSpec {
            name: "dose value",
            values: &dose,
        },
        ExpandedNumericSpec {
            name: "dose.value_as.numeric",
            values: &dose,
        },
    ];
    let derived_numeric_alias_error = expanded_formula_design(
        "~ as.numeric(`dose value`)",
        &[],
        &numeric_transform_alias_collision,
    )
    .unwrap_err()
    .to_string();
    assert!(derived_numeric_alias_error.contains("resolves ambiguously"));
}

#[test]
fn expanded_formula_design_simplifies_duplicate_formula_terms() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let batch = vec![
        "X".to_string(),
        "Y".to_string(),
        "X".to_string(),
        "Y".to_string(),
    ];
    let dose = [0.0, 1.0, 2.0, 3.0];
    let factors = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
            levels: None,
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
            levels: None,
        },
    ];
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];

    let simplified = expanded_formula_design(
        "~ condition + condition + dose + dose + condition:batch + batch:condition + condition:dose + condition:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    let explicit = expanded_formula_design(
        "~ condition + dose + condition:batch + condition:dose",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(simplified, explicit);

    let star_simplified =
        expanded_formula_design("~ condition * batch + condition:batch", &factors, &numeric)
            .unwrap();
    let star_explicit = expanded_formula_design("~ condition * batch", &factors, &numeric).unwrap();
    assert_eq!(star_simplified, star_explicit);

    let higher_order = expanded_formula_design(
        "~ condition:batch:dose + dose:condition:batch",
        &factors,
        &numeric,
    )
    .unwrap();
    let higher_order_explicit =
        expanded_formula_design("~ condition:batch:dose", &factors, &numeric).unwrap();
    assert_eq!(higher_order, higher_order_explicit);
}

#[test]
fn expanded_formula_design_from_model_frame_keeps_unused_declared_reference() {
    let condition = vec!["B".to_string(), "C".to_string(), "B".to_string()];
    let levels = vec!["A".to_string(), "B".to_string(), "C".to_string()];
    let factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "A",
        levels: Some(&levels),
    }];
    let direct = expanded_formula_design("~ condition", &factors, &[]).unwrap();
    assert_eq!(
        direct.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "condition_C_vs_A".to_string(),
        ]
    );
    assert_eq!(
        direct.standard_design.matrix().as_slice(),
        &[1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0]
    );

    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: condition,
            levels: Some(vec!["A".to_string(), "B".to_string(), "C".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let from_model_frame =
        expanded_formula_design_from_model_frame("~ condition", &model_frame).unwrap();
    assert_eq!(from_model_frame, direct);
    assert_eq!(model_frame.n_samples().unwrap(), 3);
    assert_eq!(
        model_frame.resolved_factor_reference("condition").unwrap(),
        Some("A")
    );
    assert_eq!(
        model_frame.resolved_factor_references().unwrap(),
        vec![ResolvedFormulaFactorReference {
            factor: "condition",
            reference: "A",
            levels: Some(&levels),
        }]
    );
}

#[test]
fn expanded_formula_design_from_model_frame_validates_declared_factor_levels() {
    let duplicate_level = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string()],
            levels: Some(vec!["A".to_string(), "B".to_string(), "A".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    assert!(
        expanded_formula_design_from_model_frame("~ condition", &duplicate_level)
            .unwrap_err()
            .to_string()
            .contains("duplicate declared level")
    );

    let missing_sample_level = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string(), "C".to_string()],
            levels: Some(vec!["A".to_string(), "B".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    assert!(
        expanded_formula_design_from_model_frame("~ condition", &missing_sample_level)
            .unwrap_err()
            .to_string()
            .contains("declared levels")
    );

    let missing_reference = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string()],
            levels: Some(vec!["A".to_string(), "B".to_string()]),
            reference: Some("C".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };
    assert!(
        expanded_formula_design_from_model_frame("~ condition", &missing_reference)
            .unwrap_err()
            .to_string()
            .contains("reference level")
    );
}

#[test]
fn expanded_formula_design_from_model_frame_handles_offsets_and_validation() {
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string(), "A".to_string()],
            levels: None,
            reference: Some("A".to_string()),
        }],
        numeric_covariates: vec![
            FormulaNumericColumn {
                name: "dose".to_string(),
                values: vec![1.0, 2.0, 4.0],
            },
            FormulaNumericColumn {
                name: "offset_log".to_string(),
                values: vec![0.1, 0.2, 0.3],
            },
        ],
    };

    let parsed = expanded_formula_design_with_offsets_from_model_frame(
        "~ condition + log2(dose) + offset(offset_log)",
        &model_frame,
    )
    .unwrap();
    assert_eq!(parsed.offsets, vec![0.1, 0.2, 0.3]);
    assert_eq!(
        parsed.design.standard_design.coefficient_names().unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "dose_log2".to_string(),
        ]
    );

    let named_transforms = expanded_formula_design_with_offsets_from_model_frame(
        "~ condition + I(x=dose + 1) + log(x=dose, base=2) + scale(x=dose, center=FALSE, scale=FALSE) + poly(x=dose, degree=2) + offset(I(x=offset_log + dose))",
        &model_frame,
    )
    .unwrap();
    assert_eq!(named_transforms.offsets, vec![1.1, 2.2, 4.3],);
    assert_eq!(
        named_transforms
            .design
            .standard_design
            .coefficient_names()
            .unwrap(),
        &[
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "dose_plus_1".to_string(),
            "dose_log_base_2".to_string(),
            "dose_identity".to_string(),
            "poly(dose, 2)1".to_string(),
            "poly(dose, 2)2".to_string(),
        ]
    );
    for (sample, dose) in [1.0_f64, 2.0, 4.0].iter().copied().enumerate() {
        let row = named_transforms
            .design
            .standard_design
            .matrix()
            .row(sample)
            .unwrap();
        assert_eq!(row[2], dose + 1.0);
        assert!((row[3] - dose.log2()).abs() < 1e-12);
        assert_eq!(row[4], dose);
    }

    let duplicate = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string()],
            levels: None,
            reference: None,
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "condition".to_string(),
            values: vec![0.0, 1.0],
        }],
    };
    let duplicate_error = expanded_formula_design_from_model_frame("~ condition", &duplicate)
        .unwrap_err()
        .to_string();
    assert!(duplicate_error.contains("cannot be both a factor and numeric covariate"));

    let ambiguous_cross_type = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string()],
            levels: None,
            reference: None,
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "cell-type".to_string(),
            values: vec![0.0, 1.0],
        }],
    };
    let ambiguous_error =
        expanded_formula_design_from_model_frame("~ cell.type", &ambiguous_cross_type)
            .unwrap_err()
            .to_string();
    assert!(
        ambiguous_error
            .contains("formula column 'cell-type' resolves ambiguously after R-style cleanup")
    );

    let bad_reference = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: vec!["A".to_string(), "B".to_string()],
            levels: None,
            reference: Some("C".to_string()),
        }],
        numeric_covariates: vec![],
    };
    assert!(expanded_formula_design_from_model_frame("~ condition", &bad_reference).is_err());
}

#[test]
fn formula_interactions_are_idempotent_before_expansion() {
    // Independently checked with base R 4.6.2 model.matrix: all these formulas
    // have the same intercept/x columns, assignments and order as ~ x.
    let x = [1., 2., 3., 4.];
    let numeric = [ExpandedNumericSpec {
        name: "x",
        values: &x,
    }];
    let expected = expanded_formula_design("~ x", &[], &numeric).unwrap();
    for formula in [
        "~ x*x",
        "~ x:x",
        "~ (x+x)*x",
        "~ x/x",
        "~ (x+x)^2",
        "~ x:`x`",
    ] {
        assert_eq!(
            expanded_formula_design(formula, &[], &numeric).unwrap(),
            expected,
            "{formula}"
        );
    }
    let repeated = format!("~ {}", vec!["x"; 100].join("*"));
    assert_eq!(
        expanded_formula_design(&repeated, &[], &numeric).unwrap(),
        expected
    );
    let names = (0..20).map(|i| format!("v{i}")).collect::<Vec<_>>();
    let covariates = names
        .iter()
        .map(|name| ExpandedNumericSpec { name, values: &x })
        .collect::<Vec<_>>();
    let explosive = format!("~ {}", names.join("*"));
    assert!(
        expanded_formula_design(&explosive, &[], &covariates)
            .unwrap_err()
            .to_string()
            .contains("formula expansion exceeds")
    );
}
