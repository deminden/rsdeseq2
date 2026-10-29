use super::*;

#[test]
fn fixed_dispersion_wald_list_contrast_replacement_accepts_cleaned_positive_and_negative_aliases() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        2,
        6,
        vec![
            8, 9, 14, 16, 30, 33, //
            30, 36, 50, 58, 80, 91,
        ],
        Some(vec!["g1".into(), "g2".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        6,
        3,
        vec![
            1.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, //
            1.0, 1.0, 0.0, //
            1.0, 0.0, 1.0, //
            1.0, 0.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "cell type_B cell_vs_T cell".into(),
            "cell type_NK cell_vs_T cell".into(),
        ]),
    )
    .unwrap();
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 6])
        .disable_cooks_cutoff()
        .disable_independent_filtering();
    let options = CooksReplacementOptions::new(f64::MAX);

    let cleaned = builder
        .clone()
        .fit_fixed_dispersion_wald_results_contrast_with_cooks_replacement(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::list(
                vec!["cell.type_NK.cell_vs_T.cell".into()],
                vec!["cell.type_B.cell_vs_T.cell".into()],
            ),
            None::<&[String]>,
            &options,
        )
        .unwrap();
    let raw = builder
        .fit_fixed_dispersion_wald_results_contrast_with_cooks_replacement(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::list(
                vec!["cell type_NK cell_vs_T cell".into()],
                vec!["cell type_B cell_vs_T cell".into()],
            ),
            None::<&[String]>,
            &options,
        )
        .unwrap();

    assert_eq!(cleaned.refit_plan, raw.refit_plan);
    assert_eq!(cleaned.results.rows, raw.results.rows);
    assert_eq!(
        cleaned.results.metadata.comparison.as_deref(),
        Some(
            "coefficient list contrast: cell.type_NK.cell_vs_T.cell vs cell.type_B.cell_vs_T.cell"
        )
    );
    assert_eq!(
        raw.results.metadata.comparison.as_deref(),
        Some(
            "coefficient list contrast: cell type_NK cell_vs_T cell vs cell type_B cell_vs_T cell"
        )
    );
}

#[test]
fn wald_replacement_numeric_contrast_validates_stored_model_frame_metadata() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        2,
        vec![
            10, 20, //
            20, 10,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        2,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let ambiguous = FormulaModelFrame {
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
    let options = CooksReplacementOptions::new(f64::MAX);

    let error = DeseqBuilder::new()
        .size_factors(vec![1.0; 2])
        .model_frame(ambiguous)
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald_results_contrast_with_cooks_replacement(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::numeric(vec![0.0, 1.0]),
            None::<&[String]>,
            &options,
        )
        .unwrap_err()
        .to_string();

    assert!(error.contains("resolves ambiguously after R-style cleanup"));
}

#[test]
fn fixed_dispersion_wald_replacement_uses_stored_formula_low_count_cooks_gate() {
    let counts = CountMatrix::from_row_major_u32(1, 6, vec![1, 20, 21, 20, 20, 20]).unwrap();
    let design = DesignMatrix::from_row_major(
        6,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let levels = ["A", "A", "A", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels,
            levels: Some(vec!["A".to_string(), "B".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = DeseqBuilder::new()
        .fit_type(FitType::Mean)
        .size_factors(vec![1.0; 6])
        .cooks_cutoff_threshold(0.0)
        .disable_independent_filtering();

    let generic = builder
        .clone()
        .fit_wald_glm_mu_with_cooks_replacement(
            &counts,
            &design,
            1,
            &CooksReplacementOptions::new(f64::MAX),
        )
        .unwrap();
    let formula = builder
        .model_frame(model_frame)
        .fit_wald_glm_mu_with_cooks_replacement(
            &counts,
            &design,
            1,
            &CooksReplacementOptions::new(f64::MAX),
        )
        .unwrap();

    assert!(formula.results.rows[0].max_cooks.unwrap() > 0.0);
    assert_eq!(generic.results.rows[0].cooks_outlier, Some(true));
    assert_eq!(generic.results.rows[0].pvalue, None);
    assert_eq!(formula.results.rows[0].cooks_outlier, Some(false));
    assert!(formula.results.rows[0].pvalue.is_some());
}

#[test]
fn fixed_dispersion_wald_named_replacement_uses_stored_formula_low_count_cooks_gate() {
    let counts = CountMatrix::from_row_major_u32(1, 6, vec![1, 20, 21, 20, 20, 20]).unwrap();
    let design = DesignMatrix::from_row_major(
        6,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let levels = ["A", "A", "A", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels,
            levels: Some(vec!["A".to_string(), "B".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 6])
        .cooks_cutoff_threshold(0.0)
        .disable_independent_filtering();
    let contrast = ContrastSpec::coefficient_name("condition_B_vs_A");
    let options = CooksReplacementOptions::new(f64::MAX);

    let generic = builder
        .clone()
        .fit_fixed_dispersion_wald_contrast_spec_with_cooks_replacement(
            &counts,
            &design,
            &[0.1],
            &contrast,
            &options,
        )
        .unwrap();
    let formula = builder
        .model_frame(model_frame)
        .fit_fixed_dispersion_wald_contrast_spec_with_cooks_replacement(
            &counts,
            &design,
            &[0.1],
            &contrast,
            &options,
        )
        .unwrap();

    assert!(formula.results.rows[0].max_cooks.unwrap() > 0.0);
    assert_eq!(generic.results.rows[0].cooks_outlier, Some(true));
    assert_eq!(generic.results.rows[0].pvalue, None);
    assert_eq!(formula.results.rows[0].cooks_outlier, Some(false));
    assert!(formula.results.rows[0].pvalue.is_some());
    assert_eq!(
        formula.results.metadata.comparison.as_deref(),
        Some("coefficient condition_B_vs_A")
    );
}

#[test]
fn fixed_dispersion_wald_reverse_numeric_replacement_uses_stored_formula_low_count_cooks_gate() {
    let counts = CountMatrix::from_row_major_u32(1, 6, vec![1, 20, 21, 20, 20, 20]).unwrap();
    let design = DesignMatrix::from_row_major(
        6,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let levels = ["A", "A", "A", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels,
            levels: Some(vec!["A".to_string(), "B".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 6])
        .cooks_cutoff_threshold(0.0)
        .disable_independent_filtering();
    let options = CooksReplacementOptions::new(f64::MAX);

    let generic = builder
        .clone()
        .fit_fixed_dispersion_wald_contrast_with_cooks_replacement(
            &counts,
            &design,
            &[0.1],
            &[0.0, -1.0],
            &options,
        )
        .unwrap();
    let formula = builder
        .model_frame(model_frame)
        .fit_fixed_dispersion_wald_contrast_with_cooks_replacement(
            &counts,
            &design,
            &[0.1],
            &[0.0, -1.0],
            &options,
        )
        .unwrap();

    assert!(formula.results.rows[0].max_cooks.unwrap() > 0.0);
    assert_eq!(generic.results.rows[0].cooks_outlier, Some(true));
    assert_eq!(generic.results.rows[0].pvalue, None);
    assert_eq!(formula.results.rows[0].cooks_outlier, Some(false));
    assert!(formula.results.rows[0].pvalue.is_some());
}
