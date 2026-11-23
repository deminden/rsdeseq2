use super::*;
use crate::support::run::{run_replacement, run_results, run_wald_replacement};
use rsdeseq2::prelude::{Analysis, MeanModel, ResultColumn};

#[test]
fn native_glm_mu_factor_level_contrast_applies_low_count_cooks_gate() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        3,
        6,
        vec![
            1, 20, 21, 20, 20, 20, //
            10, 11, 12, 15, 16, 17, //
            30, 32, 31, 50, 52, 51,
        ],
        Some(vec!["low_count".into(), "stable".into(), "up".into()]),
        None,
    )
    .unwrap();
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
    // This numerical fixture represents the stated original factor formula.
    let design = design.with_origin(DesignOrigin::formula(vec!["condition".into()]).unwrap());
    let levels = ["A", "A", "A", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

    let builder = DeseqBuilder::new()
        .fit_type(FitType::Mean)
        .size_factors(vec![1.0; 6])
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .cooks_cutoff_threshold(0.0)
        .disable_independent_filtering();
    let (_numeric_fit, numeric_results) = run_results(
        builder
            .clone()
            .run(&Analysis::new(&counts, &design).result(ResultColumn::Numeric(&[0.0, 1.0]))),
    )
    .unwrap();
    let (_fit, results) = run_results(builder.run(&Analysis::new(&counts, &design).result(
        ResultColumn::FactorLevel(FactorLevelContrast::new("condition", "B", "A", &levels)),
    )))
    .unwrap();

    assert_eq!(results.rows[0].gene.as_deref(), Some("low_count"));
    assert!(results.rows[0].max_cooks.unwrap() > 0.0);
    assert_eq!(numeric_results.rows[0].cooks_outlier, Some(true));
    assert_eq!(results.rows[0].cooks_outlier, Some(false));
}

#[test]
fn native_formula_wald_default_uses_stored_formula_low_count_cooks_gate() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        3,
        6,
        vec![
            1, 20, 21, 20, 20, 20, //
            10, 11, 12, 15, 16, 17, //
            30, 32, 31, 50, 52, 51,
        ],
        Some(vec!["low_count".into(), "stable".into(), "up".into()]),
        None,
    )
    .unwrap();
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
    // This numerical fixture represents the stated original factor formula.
    let design = design.with_origin(DesignOrigin::formula(vec!["condition".into()]).unwrap());
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
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .cooks_cutoff_threshold(0.0)
        .disable_independent_filtering();

    let (_generic_fit, generic_results) = run_results(
        builder.clone().run(
            &Analysis::new(&counts, &design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Coefficient(1)),
        ),
    )
    .unwrap();
    let (_formula_fit, formula_results) = run_results(
        builder
            .model_frame(model_frame.clone())
            .run(&Analysis::formula(&counts, "~ condition")),
    )
    .unwrap();

    assert_eq!(formula_results.rows[0].gene.as_deref(), Some("low_count"));
    assert!(formula_results.rows[0].max_cooks.unwrap() > 0.0);
    assert_eq!(generic_results.rows[0].cooks_outlier, Some(true));
    assert_eq!(generic_results.rows[0].pvalue, None);
    assert_eq!(formula_results.rows[0].cooks_outlier, Some(false));
    assert!(formula_results.rows[0].pvalue.is_some());
}

#[test]
fn native_formula_wald_default_uses_cleaned_stored_reference_for_cooks_gate() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        3,
        6,
        vec![
            1, 20, 21, 20, 20, 20, //
            10, 11, 12, 15, 16, 17, //
            30, 32, 31, 50, 52, 51,
        ],
        Some(vec!["low_count".into(), "stable".into(), "up".into()]),
        None,
    )
    .unwrap();
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
            "cell type_B cell_vs_A cell".into(),
        ]),
    )
    .unwrap();
    // This numerical fixture represents the stated original factor formula.
    let design = design.with_origin(DesignOrigin::formula(vec!["cell type".into()]).unwrap());
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
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .cooks_cutoff_threshold(0.0)
        .disable_independent_filtering();

    let (_generic_fit, generic_results) = run_results(
        builder.clone().run(
            &Analysis::new(&counts, &design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Coefficient(1)),
        ),
    )
    .unwrap();
    let generic_replacement = run_wald_replacement(builder.clone().run(
        &Analysis::new(&counts, &design).cooks_replacement(&CooksReplacementOptions::new(f64::MAX)),
    ))
    .unwrap();
    let (_formula_fit, formula_results) = run_results(
        builder
            .clone()
            .model_frame(model_frame.clone())
            .run(&Analysis::formula(&counts, "~ `cell type`")),
    )
    .unwrap();
    let formula_replacement = run_wald_replacement(
        builder.model_frame(model_frame).run(
            &Analysis::formula(&counts, "~ `cell type`")
                .cooks_replacement(&CooksReplacementOptions::new(f64::MAX)),
        ),
    )
    .unwrap();

    assert_eq!(formula_results.rows[0].gene.as_deref(), Some("low_count"));
    assert!(formula_results.rows[0].max_cooks.unwrap() > 0.0);
    assert_eq!(generic_results.rows[0].cooks_outlier, Some(true));
    assert_eq!(generic_results.rows[0].pvalue, None);
    assert_eq!(formula_results.rows[0].cooks_outlier, Some(false));
    assert!(formula_results.rows[0].pvalue.is_some());
    assert_eq!(
        generic_replacement.results.rows[0].cooks_outlier,
        Some(true)
    );
    assert_eq!(generic_replacement.results.rows[0].pvalue, None);
    assert_eq!(
        formula_replacement.results.rows[0].cooks_outlier,
        Some(false)
    );
    assert!(formula_replacement.results.rows[0].pvalue.is_some());
}

#[test]
fn native_replacement_named_contrasts_use_stored_formula_low_count_cooks_gate() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        3,
        6,
        vec![
            1, 20, 21, 20, 20, 20, //
            10, 11, 12, 15, 16, 17, //
            30, 32, 31, 50, 52, 51,
        ],
        Some(vec!["low_count".into(), "stable".into(), "up".into()]),
        None,
    )
    .unwrap();
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
    let reduced_design = DesignMatrix::intercept_only(counts.n_samples()).unwrap();
    // This numerical fixture represents the stated original factor formula.
    let design = design.with_origin(DesignOrigin::formula(vec!["condition".into()]).unwrap());
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
    let replacement_options = CooksReplacementOptions::new(f64::MAX);
    let builder = DeseqBuilder::new()
        .fit_type(FitType::Mean)
        .size_factors(vec![1.0; 6])
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .cooks_cutoff_threshold(0.0)
        .disable_independent_filtering();
    let contrast = ContrastSpec::coefficient_name("condition_B_vs_A");

    let generic_wald = run_wald_replacement(
        builder.clone().run(
            &Analysis::new(&counts, &design)
                .result(ResultColumn::Spec(&contrast))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();
    let formula_wald = run_wald_replacement(
        builder.clone().model_frame(model_frame.clone()).run(
            &Analysis::new(&counts, &design)
                .result(ResultColumn::Spec(&contrast))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();

    assert_eq!(
        generic_wald.results.rows[0].gene.as_deref(),
        Some("low_count")
    );
    assert_eq!(generic_wald.results.rows[0].cooks_outlier, Some(true));
    assert_eq!(formula_wald.results.rows[0].cooks_outlier, Some(false));
    assert!(formula_wald.results.rows[0].pvalue.is_some());

    let generic_reverse_wald = run_wald_replacement(
        builder.clone().run(
            &Analysis::new(&counts, &design)
                .result(ResultColumn::Numeric(&[0.0, -1.0]))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();
    let formula_reverse_wald = run_wald_replacement(
        builder.clone().model_frame(model_frame.clone()).run(
            &Analysis::new(&counts, &design)
                .result(ResultColumn::Numeric(&[0.0, -1.0]))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();

    assert_eq!(
        generic_reverse_wald.results.rows[0].gene.as_deref(),
        Some("low_count")
    );
    assert_eq!(
        generic_reverse_wald.results.rows[0].cooks_outlier,
        Some(true)
    );
    assert_eq!(
        formula_reverse_wald.results.rows[0].cooks_outlier,
        Some(false)
    );
    assert!(formula_reverse_wald.results.rows[0].pvalue.is_some());

    let generic_lrt = run_replacement(
        builder.clone().run(
            &Analysis::new(&counts, &design)
                .lrt(&reduced_design)
                .result(ResultColumn::Spec(&contrast))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();
    let formula_lrt = run_replacement(
        builder.clone().model_frame(model_frame.clone()).run(
            &Analysis::new(&counts, &design)
                .lrt(&reduced_design)
                .result(ResultColumn::Spec(&contrast))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();
    let CooksReplacementTestOutput::Lrt(generic_lrt) = generic_lrt else {
        panic!("named contrast should route to LRT");
    };
    let CooksReplacementTestOutput::Lrt(formula_lrt) = formula_lrt else {
        panic!("named contrast should route to LRT");
    };

    assert_eq!(
        generic_lrt.results.rows[0].gene.as_deref(),
        Some("low_count")
    );
    assert_eq!(generic_lrt.results.rows[0].cooks_outlier, Some(true));
    assert_eq!(formula_lrt.results.rows[0].cooks_outlier, Some(false));
    assert!(formula_lrt.results.rows[0].pvalue.is_some());

    let generic_reverse_lrt = run_replacement(
        builder.clone().run(
            &Analysis::new(&counts, &design)
                .lrt(&reduced_design)
                .result(ResultColumn::Numeric(&[0.0, -1.0]))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();
    let formula_reverse_lrt = run_replacement(
        builder.model_frame(model_frame).run(
            &Analysis::new(&counts, &design)
                .lrt(&reduced_design)
                .result(ResultColumn::Numeric(&[0.0, -1.0]))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();
    let CooksReplacementTestOutput::Lrt(generic_reverse_lrt) = generic_reverse_lrt else {
        panic!("reverse numeric contrast should route to LRT");
    };
    let CooksReplacementTestOutput::Lrt(formula_reverse_lrt) = formula_reverse_lrt else {
        panic!("reverse numeric contrast should route to LRT");
    };

    assert_eq!(
        generic_reverse_lrt.results.rows[0].gene.as_deref(),
        Some("low_count")
    );
    assert_eq!(
        generic_reverse_lrt.results.rows[0].cooks_outlier,
        Some(true)
    );
    assert_eq!(
        formula_reverse_lrt.results.rows[0].cooks_outlier,
        Some(false)
    );
    assert!(formula_reverse_lrt.results.rows[0].pvalue.is_some());
}

#[test]
fn native_replacement_numeric_contrasts_do_not_infer_low_count_gate_from_multifactor_frame() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        3,
        6,
        vec![
            1, 20, 21, 20, 20, 20, //
            10, 11, 12, 15, 16, 17, //
            30, 32, 31, 50, 52, 51,
        ],
        Some(vec!["low_count".into(), "stable".into(), "up".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        6,
        3,
        vec![
            1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, //
            1.0, 0.0, 0.0, //
            1.0, 1.0, 1.0, //
            1.0, 1.0, 0.0, //
            1.0, 1.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "condition_B_vs_A".into(),
            "batch_Y_vs_X".into(),
        ]),
    )
    .unwrap();
    let model_frame = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "condition".to_string(),
                sample_levels: ["A", "A", "A", "B", "B", "B"]
                    .into_iter()
                    .map(String::from)
                    .collect(),
                levels: Some(vec!["A".to_string(), "B".to_string()]),
                reference: None,
            },
            FormulaFactorColumn {
                name: "batch".to_string(),
                sample_levels: ["X", "Y", "X", "Y", "X", "Y"]
                    .into_iter()
                    .map(String::from)
                    .collect(),
                levels: Some(vec!["X".to_string(), "Y".to_string()]),
                reference: None,
            },
        ],
        numeric_covariates: Vec::new(),
    };
    let replacement_options = CooksReplacementOptions::new(f64::MAX);
    let builder = DeseqBuilder::new()
        .fit_type(FitType::Mean)
        .size_factors(vec![1.0; 6])
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .cooks_cutoff_threshold(0.0)
        .disable_independent_filtering();

    let generic = run_wald_replacement(
        builder.clone().run(
            &Analysis::new(&counts, &design)
                .result(ResultColumn::Numeric(&[0.0, 1.0, 0.0]))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();
    let formula = run_wald_replacement(
        builder.model_frame(model_frame).run(
            &Analysis::new(&counts, &design)
                .result(ResultColumn::Numeric(&[0.0, 1.0, 0.0]))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();

    assert_eq!(generic.results.rows[0].gene.as_deref(), Some("low_count"));
    assert_eq!(formula.results.rows, generic.results.rows);
    assert_eq!(formula.results.metadata, generic.results.metadata);
}
