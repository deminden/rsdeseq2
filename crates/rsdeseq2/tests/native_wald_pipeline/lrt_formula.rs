use super::*;

#[test]
fn top_level_lrt_formula_replacement_refit_retains_formula_derived_metadata() {
    let counts = native_wald_counts_with_zero_row();
    let full_design = two_group_design();
    let reduced_design = DesignMatrix::intercept_only(counts.n_samples()).unwrap();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels,
            levels: None,
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = glm_mu_native_wald_builder()
        .test(TestType::Lrt)
        .reduced_design(reduced_design.clone())
        .model_frame(model_frame);
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![true, false, false, false, false, false, false, false]),
    };

    let expected = builder
        .fit_lrt_with_results_with_cooks_replacement(
            &counts,
            &full_design,
            &reduced_design,
            &replacement_options,
        )
        .unwrap();
    let transformed = builder
        .fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ factor(condition, levels=c('A','B'))",
            "~ 1",
            &ResultsContrast::character("factor.condition..levels...c..A....B...", "B", "A"),
            &replacement_options,
        )
        .unwrap();

    assert!(transformed.refit_plan.should_refit);
    assert_eq!(transformed.refit_plan, expected.refit_plan);
    assert_result_rows_close(
        &transformed.results.rows,
        &expected.results.rows,
        "transformed LRT formula replacement results",
    );
    let transformed_model_frame = transformed.original_fit.current_model_frame().unwrap();
    assert_eq!(
        transformed_model_frame
            .factors
            .iter()
            .map(|factor| factor.name.as_str())
            .collect::<Vec<_>>(),
        vec!["condition", "factor(condition, levels = c(\"A\", \"B\"))"]
    );
    assert_eq!(
        transformed
            .refit_fit
            .as_ref()
            .and_then(DeseqFit::current_model_frame),
        Some(transformed_model_frame)
    );
    assert_eq!(transformed.results.metadata.test_type, Some(TestType::Lrt));
    assert_eq!(
        transformed.results.metadata.result_name.as_deref(),
        Some("factor.condition..levels...c..A....B..._B_vs_A")
    );
}

#[test]
fn top_level_test_output_routes_lrt_cooks_replacement_contrasts() {
    let counts = native_wald_counts_with_zero_row();
    let full_design = two_group_design();
    let reduced_design = DesignMatrix::intercept_only(counts.n_samples()).unwrap();
    let options = CooksReplacementOptions::new(f64::MAX);
    let builder = glm_mu_native_wald_builder()
        .test(TestType::Lrt)
        .reduced_design(reduced_design.clone());
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

    let numeric = builder
        .fit_with_test_results_contrast_with_cooks_replacement(
            &counts,
            &full_design,
            &[0.0, 1.0],
            &options,
        )
        .unwrap();
    let spec = builder
        .fit_with_test_results_contrast_spec_with_cooks_replacement(
            &counts,
            &full_design,
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
            &options,
        )
        .unwrap();
    let factor = builder
        .fit_with_test_results_factor_level_contrast_with_cooks_replacement(
            &counts,
            &full_design,
            FactorLevelContrast::new("condition", "B", "A", &levels),
            &options,
        )
        .unwrap();
    let request = builder
        .fit_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            &full_design,
            &ResultsContrast::character("condition", "B", "A"),
            Some(&levels),
            &options,
        )
        .unwrap();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels.clone(),
            levels: None,
            reference: None,
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "exposure".to_string(),
            values: vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
        }],
    };
    let model_frame_request = builder
        .fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement(
            &counts,
            &full_design,
            &ResultsContrast::character("condition", "B", "A"),
            &model_frame,
            &options,
        )
        .unwrap();
    let expected = builder
        .fit_lrt_with_results_contrast_with_cooks_replacement(
            &counts,
            &full_design,
            &reduced_design,
            &[0.0, 1.0],
            &options,
        )
        .unwrap();

    let CooksReplacementTestOutput::Lrt(numeric) = numeric else {
        panic!("numeric contrast should route to LRT");
    };
    let CooksReplacementTestOutput::Lrt(spec) = spec else {
        panic!("contrast spec should route to LRT");
    };
    let CooksReplacementTestOutput::Lrt(factor) = factor else {
        panic!("factor contrast should route to LRT");
    };
    let CooksReplacementTestOutput::Lrt(request) = request else {
        panic!("results contrast request should route to LRT");
    };
    let CooksReplacementTestOutput::Lrt(model_frame_request) = model_frame_request else {
        panic!("model-frame results contrast request should route to LRT");
    };

    assert_eq!(numeric.results, expected.results);
    assert_eq!(request.results, factor.results);
    assert_eq!(model_frame_request.results, factor.results);
    assert_eq!(
        model_frame_request.original_fit.current_model_frame(),
        Some(&model_frame)
    );
    assert_eq!(
        spec.results.metadata.comparison.as_deref(),
        Some("coefficient condition_B_vs_A")
    );
    assert_eq!(
        factor.results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
}

#[test]
fn top_level_lrt_formula_contrast_uses_formula_derived_factor_metadata() {
    let counts = native_wald_counts_with_zero_row();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels,
            levels: Some(vec!["A".to_string(), "B".to_string()]),
            reference: Some("A".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = glm_mu_native_wald_builder()
        .test(TestType::Lrt)
        .model_frame(model_frame);
    let formula = "~ factor(condition, levels=c('B','A'))";
    let reduced = "~ 1";
    let contrast = ResultsContrast::character("factor.condition..levels...c..B....A...", "A", "B");

    let (fit, mut results) = builder
        .fit_lrt_formula_with_results_contrast_request(&counts, formula, reduced, &contrast)
        .unwrap();
    results.independent_filtering = None;
    let object_results = fit
        .results_contrast_request(&counts, &contrast, counts.gene_names())
        .unwrap();
    let replacement = builder
        .fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
            &counts,
            formula,
            reduced,
            &contrast,
            &CooksReplacementOptions::new(f64::MAX),
        )
        .unwrap();

    assert_eq!(replacement.refit_plan.n_refit, 0);
    assert_eq!(replacement.results.rows, results.rows);
    assert_eq!(replacement.results.metadata, results.metadata);
    assert_eq!(object_results, results);
    assert_eq!(results.metadata.test_type, Some(TestType::Lrt));
    assert_eq!(
        results.metadata.result_name.as_deref(),
        Some("factor.condition..levels...c..B....A..._A_vs_B")
    );
    assert_eq!(
        results.metadata.comparison.as_deref(),
        Some("factor-level contrast: factor(condition, levels = c(\"B\", \"A\")) A vs B")
    );

    let ordered_formula = "~ ordered(condition, levels=c('B','A'))";
    let ordered_contrast = ResultsContrast::numeric(vec![0.0, 1.0]);
    let invalid = builder
        .fit_lrt_formula_with_results_contrast_request(
            &counts,
            ordered_formula,
            reduced,
            &ResultsContrast::character("ordered.condition..levels...c..B....A...", "A", "B"),
        )
        .unwrap_err();
    assert!(invalid.to_string().contains("ordered factor"));
    let (ordered_fit, mut ordered_results) = builder
        .fit_lrt_formula_with_results_contrast_request(
            &counts,
            ordered_formula,
            reduced,
            &ordered_contrast,
        )
        .unwrap();
    ordered_results.independent_filtering = None;
    let ordered_object_results = ordered_fit
        .results_contrast_request(&counts, &ordered_contrast, counts.gene_names())
        .unwrap();
    let ordered_replacement = builder
        .fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
            &counts,
            ordered_formula,
            reduced,
            &ordered_contrast,
            &CooksReplacementOptions::new(f64::MAX),
        )
        .unwrap();
    assert_eq!(ordered_replacement.refit_plan.n_refit, 0);
    assert_eq!(ordered_replacement.results.rows, ordered_results.rows);
    assert_eq!(
        ordered_replacement.results.metadata,
        ordered_results.metadata
    );
    assert_eq!(ordered_object_results, ordered_results);
    assert_eq!(
        ordered_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        ordered_results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );

    let relabeled_formula = "~ factor(condition, levels=c('A','B'), labels=c('control','treated'))";
    let relabeled_contrast = ResultsContrast::character(
        "factor.condition..levels...c..A....B....labels...c..control....treated...",
        "treated",
        "control",
    );
    let (relabeled_fit, mut relabeled_results) = builder
        .fit_lrt_formula_with_results_contrast_request(
            &counts,
            relabeled_formula,
            reduced,
            &relabeled_contrast,
        )
        .unwrap();
    relabeled_results.independent_filtering = None;
    let relabeled_object_results = relabeled_fit
        .results_contrast_request(&counts, &relabeled_contrast, counts.gene_names())
        .unwrap();
    let relabeled_replacement = builder
        .fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
            &counts,
            relabeled_formula,
            reduced,
            &relabeled_contrast,
            &CooksReplacementOptions::new(f64::MAX),
        )
        .unwrap();

    assert_eq!(relabeled_replacement.refit_plan.n_refit, 0);
    assert_eq!(relabeled_replacement.results.rows, relabeled_results.rows);
    assert_eq!(relabeled_object_results, relabeled_results);
    assert_eq!(
        relabeled_results.metadata.result_name.as_deref(),
        Some(
            "factor.condition..levels...c..A....B....labels...c..control....treated..._treated_vs_control"
        )
    );
    assert_eq!(
        relabeled_results.metadata.comparison.as_deref(),
        Some(
            "factor-level contrast: factor(condition, levels = c(\"A\", \"B\"), labels = c(\"control\", \"treated\")) treated vs control"
        )
    );

    let error = builder
        .fit_lrt_formula_with_results_contrast_request(
            &counts,
            "~ ordered(condition, levels=c('A','B'), labels=c('control','treated'))",
            reduced,
            &ResultsContrast::character(
                "ordered.condition..levels...c..A....B....labels...c..control....treated...",
                "treated",
                "control",
            ),
        )
        .unwrap_err();
    assert!(error.to_string().contains("ordered factor"));
}
