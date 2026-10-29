use super::*;

#[test]
fn fixed_dispersion_wald_contrast_applies_explicit_cooks_cutoff() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 6]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let (_fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .cooks_cutoff_threshold(0.1)
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &[0.2], &[1.0])
        .unwrap();

    assert_relative_eq!(
        results.rows[0].max_cooks.unwrap(),
        4.0 / 8.16 * 0.75,
        epsilon = 1e-12
    );
    assert_eq!(results.rows[0].cooks_outlier, Some(true));
    assert_eq!(results.rows[0].pvalue, None);
    assert_eq!(results.rows[0].padj, None);
}

#[test]
fn fixed_dispersion_wald_factor_level_contrast_applies_low_count_cooks_gate() {
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

    let (_fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0; 6])
        .cooks_cutoff_threshold(0.0)
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald_factor_level_contrast(
            &counts,
            &design,
            &[0.1],
            FactorLevelContrast::new("condition", "B", "A", &levels),
        )
        .unwrap();

    assert!(results.rows[0].max_cooks.unwrap() > 0.0);
    assert_eq!(results.rows[0].cooks_outlier, Some(false));
    assert!(results.rows[0].pvalue.is_some());
}

#[test]
fn fixed_dispersion_wald_coefficient_uses_stored_formula_low_count_cooks_gate() {
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

    let (_generic_fit, generic_results) = builder
        .clone()
        .fit_fixed_dispersion_wald(&counts, &design, &[0.1], 1)
        .unwrap();
    let (_formula_fit, formula_results) = builder
        .model_frame(model_frame)
        .fit_fixed_dispersion_wald(&counts, &design, &[0.1], 1)
        .unwrap();

    assert!(formula_results.rows[0].max_cooks.unwrap() > 0.0);
    assert_eq!(generic_results.rows[0].cooks_outlier, Some(true));
    assert_eq!(generic_results.rows[0].pvalue, None);
    assert_eq!(formula_results.rows[0].cooks_outlier, Some(false));
    assert!(formula_results.rows[0].pvalue.is_some());
}

#[test]
fn fixed_dispersion_wald_coefficient_uses_cleaned_stored_reference_for_cooks_gate() {
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
        Some(vec![
            "Intercept".into(),
            "cell.type_B.cell_vs_A.cell".into(),
        ]),
    )
    .unwrap();
    let levels = ["A cell", "A cell", "A cell", "B cell", "B cell", "B cell"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: levels,
            levels: Some(vec!["A cell".to_string(), "B cell".to_string()]),
            reference: Some("A.cell".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = DeseqBuilder::new()
        .fit_type(FitType::Mean)
        .size_factors(vec![1.0; 6])
        .cooks_cutoff_threshold(0.0)
        .disable_independent_filtering();

    let (_generic_fit, generic_results) = builder
        .clone()
        .fit_fixed_dispersion_wald(&counts, &design, &[0.1], 1)
        .unwrap();
    let (_formula_fit, formula_results) = builder
        .model_frame(model_frame)
        .fit_fixed_dispersion_wald(&counts, &design, &[0.1], 1)
        .unwrap();

    assert_eq!(generic_results.rows[0].cooks_outlier, Some(true));
    assert_eq!(generic_results.rows[0].pvalue, None);
    assert_eq!(formula_results.rows[0].cooks_outlier, Some(false));
    assert!(formula_results.rows[0].pvalue.is_some());
}

#[test]
fn fixed_dispersion_wald_pipeline_applies_explicit_cooks_cutoff() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 6]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let (_fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .cooks_cutoff_threshold(0.1)
        .fit_fixed_dispersion_wald(&counts, &design, &[0.2], 0)
        .unwrap();

    assert_relative_eq!(
        results.rows[0].max_cooks.unwrap(),
        4.0 / 8.16 * 0.75,
        epsilon = 1e-12
    );
    assert_eq!(results.rows[0].cooks_outlier, Some(true));
    assert_eq!(results.rows[0].pvalue, None);
    assert_eq!(results.rows[0].padj, None);
}

#[test]
fn fixed_dispersion_wald_pipeline_can_disable_cooks_cutoff() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 6]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let (_fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .cooks_cutoff_threshold(0.1)
        .disable_cooks_cutoff()
        .fit_fixed_dispersion_wald(&counts, &design, &[0.2], 0)
        .unwrap();

    assert_eq!(results.rows[0].cooks_outlier, None);
    assert!(results.rows[0].pvalue.is_some());
    assert!(results.rows[0].padj.is_some());
}

#[test]
fn fixed_dispersion_wald_pipeline_applies_independent_filtering() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![1, 1, 1, 100, 100, 100]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let (_fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .disable_cooks_cutoff()
        .independent_filtering_theta(vec![1.0, 1.0])
        .fit_fixed_dispersion_wald(&counts, &design, &[0.2, 0.2], 0)
        .unwrap();

    let filtering = results.independent_filtering.as_ref().unwrap();
    assert!(filtering.enabled);
    assert_eq!(filtering.filter_threshold, Some(100.0));
    assert_eq!(results.rows[0].filtered, Some(true));
    assert_eq!(results.rows[0].padj, None);
    assert_eq!(results.rows[1].filtered, Some(false));
    assert!(results.rows[1].padj.is_some());
}

#[test]
fn fixed_dispersion_wald_pipeline_can_disable_independent_filtering() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![1, 1, 1, 100, 100, 100]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let (_fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald(&counts, &design, &[0.2, 0.2], 0)
        .unwrap();

    let filtering = results.independent_filtering.as_ref().unwrap();
    assert!(!filtering.enabled);
    assert_eq!(results.rows[0].filtered, None);
    assert_eq!(results.rows[1].filtered, None);
    assert!(results.rows[0].padj.is_some());
    assert!(results.rows[1].padj.is_some());
}
