use super::*;

#[test]
fn expanded_formula_design_with_offsets_extracts_formula_offsets() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0, 1.0, 2.0, 3.0];
    let exposure = [0.5, 1.0, 1.5, 2.0];
    let batch_offset = [0.1, 0.2, 0.3, 0.4];
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
            name: "exposure",
            values: &exposure,
        },
        ExpandedNumericSpec {
            name: "batch_offset",
            values: &batch_offset,
        },
    ];

    let with_offsets = expanded_formula_design_with_offsets(
        "~ condition + dose + offset(exposure) + offset(batch_offset)",
        &factors,
        &numeric,
    )
    .unwrap();
    let without_offsets =
        expanded_formula_design("~ condition + dose", &factors, &numeric).unwrap();
    assert_eq!(with_offsets.design, without_offsets);
    assert_eq!(
        with_offsets
            .design
            .standard_design
            .origin()
            .variables()
            .unwrap(),
        &["condition", "dose", "exposure", "batch_offset"]
    );
    assert_eq!(
        with_offsets.design.standard_design.clone().origin(),
        with_offsets.design.expanded_design.origin()
    );
    assert_eq!(
        without_offsets
            .standard_design
            .origin()
            .variables()
            .unwrap(),
        &["condition", "dose"]
    );

    assert_eq!(with_offsets.offsets, vec![0.6, 1.2, 1.8, 2.4]);

    let intercept_only =
        expanded_formula_design_with_offsets("~ offset(exposure)", &[], &numeric).unwrap();
    assert_eq!(
        intercept_only
            .design
            .standard_design
            .coefficient_names()
            .unwrap(),
        &["Intercept".to_string()]
    );
    assert_eq!(intercept_only.offsets, exposure);
}

#[test]
fn expanded_formula_design_with_offsets_supports_numeric_transform_offsets() {
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [1.0_f64, 2.0, 4.0, 8.0];
    let exposure = [0.5_f64, 1.0, 1.5, 2.0];
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
            name: "exposure value",
            values: &exposure,
        },
    ];

    let transformed = expanded_formula_design_with_offsets(
        "~ condition + offset(log2(dose)) + offset(I(x=dose + `exposure value`))",
        &factors,
        &numeric,
    )
    .unwrap();
    let expected_design = expanded_formula_design("~ condition", &factors, &numeric).unwrap();
    assert_eq!(transformed.design, expected_design);
    assert_eq!(
        transformed.offsets,
        dose.iter()
            .zip(exposure.iter())
            .map(|(dose, exposure)| dose.log2() + dose + exposure)
            .collect::<Vec<_>>()
    );

    let spaced = expanded_formula_design_with_offsets(
        "~ condition + offset ( log2 ( dose ) ) + offset ( I ( x = dose + `exposure value` ) )",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(spaced.design, expected_design);
    assert_eq!(spaced.offsets, transformed.offsets);

    let quoted =
        expanded_formula_design_with_offsets("~ offset(`exposure value`)", &[], &numeric).unwrap();
    assert_eq!(quoted.offsets, exposure);

    let scaled = expanded_formula_design_with_offsets(
        "~ offset(scale(dose, center=FALSE, scale=FALSE))",
        &[],
        &numeric,
    )
    .unwrap();
    assert_eq!(scaled.offsets, dose);

    let named_transform_offsets = expanded_formula_design_with_offsets(
        "~ offset(log(x=dose, base=2)) + offset(scale(x=`exposure value`, center=FALSE, scale=FALSE))",
        &[],
        &numeric,
    )
    .unwrap();
    assert_eq!(
        named_transform_offsets.offsets,
        dose.iter()
            .zip(exposure.iter())
            .map(|(dose, exposure)| dose.log2() + exposure)
            .collect::<Vec<_>>()
    );

    let parenthesized = expanded_formula_design_with_offsets(
        "~ (condition + offset(log2(dose))) + (offset(`exposure value`) + 1)",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(parenthesized.design, expected_design);
    assert_eq!(
        parenthesized.offsets,
        dose.iter()
            .zip(exposure.iter())
            .map(|(dose, exposure)| dose.log2() + exposure)
            .collect::<Vec<_>>()
    );

    let nested_parenthesized = expanded_formula_design_with_offsets(
        "~ condition + ((offset(log2(dose)) + 1)) + (offset(`exposure value`) + (dose - dose))",
        &factors,
        &numeric,
    )
    .unwrap();
    assert_eq!(nested_parenthesized.design, expected_design);
    assert_eq!(nested_parenthesized.offsets, parenthesized.offsets);

    let coercion_offsets = expanded_formula_design_with_offsets(
        "~ offset(as.numeric(dose)) + offset(as.double(dose)) + offset(as.integer(x=`exposure value`))",
        &[],
        &numeric,
    )
    .unwrap();
    assert_eq!(
        coercion_offsets.offsets,
        // Original R truncates the fractional exposure before adding offsets.
        [2.0, 5.0, 9.0, 18.0]
    );

    let factor_offset_levels = vec!["B".to_string(), "A".to_string()];
    let factor_offsets = expanded_formula_design_with_offsets(
        "~ condition + offset(as.numeric(condition)) + offset(as.double(condition))",
        &[ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
            levels: Some(&factor_offset_levels),
        }],
        &numeric,
    )
    .unwrap();
    assert_eq!(factor_offsets.design, expected_design);
    assert_eq!(
        factor_offsets.offsets,
        condition
            .iter()
            .map(|level| if level == "B" { 2.0 } else { 4.0 })
            .collect::<Vec<_>>()
    );

    assert!(
        expanded_formula_design_with_offsets("~ offset(poly(dose, 2, raw=TRUE))", &[], &numeric)
            .is_err()
    );
}

#[test]
fn expanded_formula_design_with_offsets_carries_formula_local_model_frame() {
    let condition = vec![
        "B".to_string(),
        "A".to_string(),
        "B".to_string(),
        "A".to_string(),
    ];
    let condition_levels = vec!["Z".to_string(), "B".to_string(), "A".to_string()];
    let dose = [1.0_f64, 2.0, 4.0, 8.0];
    let factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "Z",
        levels: Some(&condition_levels),
    }];
    let numeric = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
    }];

    let formula = "~ relevel(droplevels(condition), ref='A') + log2(dose) + offset(dose)";
    let parsed = expanded_formula_design_with_offsets(formula, &factors, &numeric).unwrap();

    assert_eq!(
        parsed
            .model_frame
            .factors
            .iter()
            .map(|factor| factor.name.as_str())
            .collect::<Vec<_>>(),
        vec!["condition", "relevel(droplevels(condition), ref = \"A\")",]
    );
    assert_eq!(
        parsed.model_frame.factors[1].levels.as_deref(),
        Some(&["B".to_string(), "A".to_string()][..])
    );
    assert_eq!(
        parsed.model_frame.factors[1].reference.as_deref(),
        Some("A")
    );
    assert_eq!(
        parsed
            .model_frame
            .numeric_covariates
            .iter()
            .map(|covariate| covariate.name.as_str())
            .collect::<Vec<_>>(),
        vec!["dose", "dose_log2"]
    );
    assert_eq!(
        parsed.model_frame.numeric_covariates[1].values,
        dose.iter().map(|value| value.log2()).collect::<Vec<_>>()
    );
}

#[test]
fn formula_offset_detection_is_syntax_level() {
    assert!(!formula_has_offset_terms("~ condition + offset_like").unwrap());
    assert!(!formula_has_offset_terms("~ condition + myoffset(exposure)").unwrap());
    assert!(formula_has_offset_terms("~ condition + offset(exposure)").unwrap());
    assert!(formula_has_offset_terms("~ condition + offset(log2(exposure))").unwrap());
    assert!(formula_has_offset_terms("~ condition + offset ( log2 ( exposure ) )").unwrap());
    assert!(formula_has_offset_terms("~ condition + offset(I(exposure + 1))").unwrap());
    assert!(formula_has_offset_terms("~ condition - offset(exposure)").unwrap());
    assert!(formula_has_offset_terms("~ condition + offset( exposure )").unwrap());
    assert!(formula_has_offset_terms("~ (condition + offset(exposure))").unwrap());
    assert!(formula_has_offset_terms("~ (condition + offset(log2(exposure)))").unwrap());
    assert!(formula_has_offset_terms("~ ((condition + (offset(log2(exposure)))))").unwrap());
    assert!(formula_has_offset_terms("~ condition + offset(exposure").is_err());
    assert!(formula_has_offset_terms("condition + offset(exposure)").is_err());
}

#[test]
fn expanded_formula_design_validates_unsupported_terms() {
    let condition = vec!["A".to_string(), "B".to_string()];
    let batch = vec!["X".to_string(), "Y".to_string()];
    let dose = [0.0, 1.0];
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
    let numeric_with_zero_time = [
        ExpandedNumericSpec {
            name: "dose",
            values: &dose,
        },
        ExpandedNumericSpec {
            name: "time",
            values: &[0.0, 1.0],
        },
    ];

    assert!(expanded_formula_design("condition + dose", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ condition + missing", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design("~ condition + condition:condition", &factors, &numeric).is_err()
    );
    assert!(
        expanded_formula_design(
            "~ as.ordered(condition, levels=c('A','B'))",
            &factors,
            &numeric
        )
        .is_err()
    );
    assert!(expanded_formula_design("~ condition - missing", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ condition -", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ (condition + batch", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design(
            "~ condition + (missing + batch - dose)^2",
            &factors,
            &numeric
        )
        .is_err()
    );
    assert!(expanded_formula_design("~ I(dose^x)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ I(dose^32)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ I(dose + missing)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ I(y=dose)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ I(dose / 0)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ I(1 / dose)", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design("~ I(dose / time)", &factors, &numeric_with_zero_time).is_err()
    );
    assert!(expanded_formula_design("~ log(missing)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ log(dose + 1)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ log(dose, base=1)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ log(dose, base=-2)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ log(dose, base=maybe)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ log(x=dose, x=dose)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ log2(dose, base=2)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ sqrt(dose, extra=1)", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design(
            "~ log1p(dose)",
            &factors,
            &[ExpandedNumericSpec {
                name: "dose",
                values: &[-1.0, 0.0],
            }]
        )
        .is_err()
    );
    assert!(expanded_formula_design("~ log2(dose", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(missing)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, center=maybe)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, scale=maybe)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, center=NaN)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, scale=0)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ scale(dose, raw=TRUE)", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design(
            "~ scale(dose, center=FALSE, center=TRUE)",
            &factors,
            &numeric
        )
        .is_err()
    );
    assert!(
        expanded_formula_design("~ scale(dose, scale=FALSE, scale=TRUE)", &factors, &numeric)
            .is_err()
    );
    assert!(expanded_formula_design("~ scale(x=dose, x=dose)", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design(
            "~ scale(dose)",
            &factors,
            &[ExpandedNumericSpec {
                name: "dose",
                values: &[1.0, 1.0],
            }]
        )
        .is_err()
    );
    assert!(
        expanded_formula_design(
            "~ poly(dose, 2)",
            &factors,
            &[ExpandedNumericSpec {
                name: "dose",
                values: &[0.0, 0.0],
            }]
        )
        .is_err()
    );
    assert!(expanded_formula_design("~ poly(dose, x, raw=TRUE)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ poly(dose, 32, raw=TRUE)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ poly(missing, 2, raw=TRUE)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ poly(dose, 2, coefs=1)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ poly(dose, 2, simple=maybe)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ poly(dose, 2, raw=maybe)", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design("~ poly(x=dose, x=dose, degree=2)", &factors, &numeric).is_err()
    );
    assert!(
        expanded_formula_design("~ poly(dose, degree=2, degree=3)", &factors, &numeric).is_err()
    );
    assert!(
        expanded_formula_design("~ poly(dose, 2, raw=TRUE, raw=FALSE)", &factors, &numeric)
            .is_err()
    );
    assert!(expanded_formula_design("~ factor(dose)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ factor(missing)", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design("~ factor(condition, levels=c('A'))", &factors, &numeric).is_err()
    );
    assert!(
        expanded_formula_design("~ factor(condition, levels=c('A','A'))", &factors, &numeric)
            .is_err()
    );
    assert!(
        expanded_formula_design("~ factor(condition, levels='A')", &factors, &numeric).is_err()
    );
    assert!(
        expanded_formula_design(
            "~ factor(condition, labels=c('control','treated'))",
            &factors,
            &numeric
        )
        .is_err()
    );
    assert!(
        expanded_formula_design(
            "~ factor(condition, levels=c('A','B'), labels=c('control'))",
            &factors,
            &numeric
        )
        .is_err()
    );
    assert!(
        expanded_formula_design(
            "~ factor(condition, levels=c('A','B'), labels=c('group','group'))",
            &factors,
            &numeric
        )
        .is_err()
    );
    assert!(
        expanded_formula_design(
            "~ as.factor(condition, levels=c('A','B'))",
            &factors,
            &numeric
        )
        .is_err()
    );
    assert!(
        expanded_formula_design(
            "~ as.factor(condition, labels=c('control','treated'))",
            &factors,
            &numeric
        )
        .is_err()
    );
    assert!(
        expanded_formula_design(
            "~ as.ordered(condition, labels=c('control','treated'))",
            &factors,
            &numeric
        )
        .is_err()
    );
    assert!(expanded_formula_design("~ droplevels(dose)", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design("~ droplevels(condition, except='A')", &factors, &numeric).is_err()
    );
    assert!(expanded_formula_design("~ factor(x=condition, x=batch)", &factors, &numeric).is_err());
    assert!(
        expanded_formula_design("~ as.factor(condition, extra=1)", &factors, &numeric).is_err()
    );
    assert!(expanded_formula_design("~ relevel(dose, ref='A')", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ relevel(condition)", &factors, &numeric).is_err());
    assert!(expanded_formula_design("~ relevel(condition, ref='Z')", &factors, &numeric).is_err());
    let ambiguous_levels = vec![
        "B cell".to_string(),
        "B-cell".to_string(),
        "T cell".to_string(),
    ];
    let ambiguous_factor = [ExpandedFactorSpec {
        factor: "cell type",
        sample_levels: &ambiguous_levels,
        reference: "T cell",
        levels: Some(&ambiguous_levels),
    }];
    let ambiguous_relevel = expanded_formula_design(
        "~ relevel(`cell type`, ref='B.cell')",
        &ambiguous_factor,
        &[],
    )
    .unwrap_err()
    .to_string();
    assert!(ambiguous_relevel.contains("resolves ambiguously"));
    assert!(
        expanded_formula_design_with_offsets("~ condition - offset(dose)", &factors, &numeric)
            .is_err()
    );
    assert!(
        expanded_formula_design_with_offsets(
            "~ condition - (offset(dose) + 1)",
            &factors,
            &numeric
        )
        .is_err()
    );
    assert!(
        expanded_formula_design_with_offsets("~ condition + offset(missing)", &factors, &numeric)
            .is_err()
    );
    assert!(
        expanded_formula_design_with_offsets("~ condition + offset(dose + 1)", &factors, &numeric)
            .is_err()
    );
    assert!(
        expanded_formula_design("~ condition * batch * dose * missing", &factors, &numeric)
            .is_err()
    );
    assert!(expanded_formula_design("~ condition:batch:dose:missing", &factors, &numeric).is_err());
}
