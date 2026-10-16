use super::*;

#[test]
fn top_level_fit_contrast_specs_return_default_automatic_wald_fits() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().automatic_dispersion_means();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

    let named_fit = builder
        .fit_contrast_spec(
            &counts,
            &design,
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
        )
        .unwrap();
    let factor_fit = builder
        .fit_factor_level_contrast(
            &counts,
            &design,
            FactorLevelContrast::new("condition", "B", "A", &levels),
        )
        .unwrap();
    let (expected_fit, _expected_results) = builder
        .fit_wald_glm_mu_contrast_spec(
            &counts,
            &design,
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
        )
        .unwrap();

    assert_eq!(named_fit.wald, expected_fit.wald);
    assert_eq!(factor_fit.wald, expected_fit.wald);
    assert_eq!(named_fit.beta_converged, expected_fit.beta_converged);
    assert_eq!(factor_fit.beta_converged, expected_fit.beta_converged);
    assert_wald_fit_state_matches(&named_fit, &expected_fit, "top-level named Wald contrast");
    assert_wald_fit_state_matches(&factor_fit, &expected_fit, "top-level factor Wald contrast");
}

#[test]
fn top_level_fit_with_results_contrast_spec_cooks_replacement_preserves_metadata() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: f64::MAX,
        min_replicates: 3,
        which_samples: None,
    };

    let named_output = builder
        .fit_with_results_contrast_spec_with_cooks_replacement(
            &counts,
            &design,
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
            &replacement_options,
        )
        .unwrap();
    let factor_output = builder
        .fit_with_results_factor_level_contrast_with_cooks_replacement(
            &counts,
            &design,
            FactorLevelContrast::new("condition", "B", "A", &levels),
            &replacement_options,
        )
        .unwrap();
    let request_output = builder
        .fit_with_test_results_contrast_request_with_cooks_replacement(
            &counts,
            &design,
            &ResultsContrast::character("condition", "B", "A"),
            Some(&levels),
            &replacement_options,
        )
        .unwrap();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels.clone(),
            levels: None,
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let model_frame_request_output = builder
        .fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement(
            &counts,
            &design,
            &ResultsContrast::character("condition", "B", "A"),
            &model_frame,
            &replacement_options,
        )
        .unwrap();
    let stored_model_frame_request_output = builder
        .clone()
        .model_frame(model_frame.clone())
        .fit_with_test_results_contrast_request_with_cooks_replacement::<String>(
            &counts,
            &design,
            &ResultsContrast::character("condition", "B", "A"),
            None,
            &replacement_options,
        )
        .unwrap();
    let (expected_fit, _expected_results) = builder
        .fit_wald_glm_mu_contrast_spec(
            &counts,
            &design,
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
        )
        .unwrap();
    let CooksReplacementTestOutput::Wald(request_output) = request_output else {
        panic!("results contrast request should route to Wald");
    };
    let CooksReplacementTestOutput::Wald(model_frame_request_output) = model_frame_request_output
    else {
        panic!("model-frame results contrast request should route to Wald");
    };
    let CooksReplacementTestOutput::Wald(stored_model_frame_request_output) =
        stored_model_frame_request_output
    else {
        panic!("stored model-frame results contrast request should route to Wald");
    };

    assert_eq!(named_output.refit_plan.n_refit, 0);
    assert_eq!(factor_output.refit_plan.n_refit, 0);
    assert_eq!(request_output.refit_plan.n_refit, 0);
    assert_eq!(model_frame_request_output.refit_plan.n_refit, 0);
    assert_eq!(
        model_frame_request_output
            .original_fit
            .current_model_frame(),
        Some(&model_frame)
    );
    assert_wald_fit_state_matches(
        &named_output.original_fit,
        &expected_fit,
        "top-level named Wald contrast no-refit original",
    );
    assert_wald_fit_state_matches(
        &factor_output.original_fit,
        &expected_fit,
        "top-level factor Wald contrast no-refit original",
    );
    assert_eq!(
        named_output.results.metadata.comparison.as_deref(),
        Some("coefficient condition_B_vs_A")
    );
    assert_eq!(
        factor_output.results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
    assert_eq!(request_output.results, factor_output.results);
    assert_eq!(model_frame_request_output.results, factor_output.results);
    assert_eq!(
        stored_model_frame_request_output.results,
        factor_output.results
    );
}

#[test]
fn top_level_fit_keeps_lrt_explicit_until_reduced_design_is_supplied() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let err = glm_mu_native_wald_builder()
        .test(TestType::Lrt)
        .fit(&counts, &design)
        .unwrap_err();

    assert!(matches!(err, DeseqError::UnsupportedFeature { .. }));
}

#[test]
fn top_level_fit_uses_stored_reduced_design_for_lrt() {
    let counts = native_wald_counts_with_zero_row();
    let full_design = two_group_design();
    let reduced_design = DesignMatrix::intercept_only(counts.n_samples()).unwrap();
    let builder = glm_mu_native_wald_builder()
        .test(TestType::Lrt)
        .reduced_design(reduced_design.clone())
        .automatic_dispersion_means();

    let top_level_fit = builder.fit(&counts, &full_design).unwrap();
    let (expected_fit, _expected_results) = builder
        .fit_lrt_glm_mu(&counts, &full_design, &reduced_design, 1)
        .unwrap();

    assert_eq!(builder.current_reduced_design(), Some(&reduced_design));
    assert_eq!(top_level_fit.lrt, expected_fit.lrt);
    assert_eq!(top_level_fit.reduced_design, expected_fit.reduced_design);
    assert_wald_fit_state_matches(&top_level_fit, &expected_fit, "top-level LRT");
    assert_matrix_close_or_nan(
        top_level_fit.reduced_mu.as_ref().unwrap(),
        expected_fit.reduced_mu.as_ref().unwrap(),
        "top-level LRT reduced mu",
    );
}

#[test]
fn top_level_fit_with_results_name_uses_stored_reduced_design_for_lrt() {
    let counts = native_wald_counts_with_zero_row();
    let full_design = two_group_design();
    let reduced_design = DesignMatrix::intercept_only(counts.n_samples()).unwrap();
    let builder = glm_mu_native_wald_builder()
        .test(TestType::Lrt)
        .reduced_design(reduced_design.clone())
        .automatic_dispersion_means();

    let (top_level_fit, top_level_results) = builder
        .fit_with_results_name(&counts, &full_design, "condition_B_vs_A")
        .unwrap();
    let (expected_fit, expected_results) = builder
        .fit_lrt_glm_mu(&counts, &full_design, &reduced_design, 1)
        .unwrap();

    assert_eq!(top_level_fit.lrt, expected_fit.lrt);
    assert_eq!(top_level_results, expected_results);
    assert_eq!(top_level_results.metadata.test_type, Some(TestType::Lrt));
}

#[test]
fn top_level_fit_with_results_contrast_keeps_lrt_explicit() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let err = glm_mu_native_wald_builder()
        .test(TestType::Lrt)
        .fit_with_results_contrast(&counts, &design, &[0.0, 1.0])
        .unwrap_err();
    let spec_err = glm_mu_native_wald_builder()
        .test(TestType::Lrt)
        .fit_with_results_contrast_spec(
            &counts,
            &design,
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
        )
        .unwrap_err();
    let factor_err = glm_mu_native_wald_builder()
        .test(TestType::Lrt)
        .fit_with_results_factor_level_contrast(
            &counts,
            &design,
            FactorLevelContrast::new("condition", "B", "A", &levels),
        )
        .unwrap_err();

    assert!(matches!(err, DeseqError::UnsupportedFeature { .. }));
    assert!(matches!(spec_err, DeseqError::UnsupportedFeature { .. }));
    assert!(matches!(factor_err, DeseqError::UnsupportedFeature { .. }));
}

#[test]
fn top_level_lrt_contrast_methods_use_stored_reduced_design() {
    let counts = native_wald_counts_with_zero_row();
    let full_design = two_group_design();
    let reduced_design = DesignMatrix::intercept_only(counts.n_samples()).unwrap();
    let builder = glm_mu_native_wald_builder()
        .test(TestType::Lrt)
        .reduced_design(reduced_design.clone())
        .automatic_dispersion_means();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

    let (numeric_fit, numeric_results) = builder
        .fit_with_results_contrast(&counts, &full_design, &[0.0, 1.0])
        .unwrap();
    let (spec_fit, spec_results) = builder
        .fit_with_results_contrast_spec(
            &counts,
            &full_design,
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
        )
        .unwrap();
    let (factor_fit, factor_results) = builder
        .fit_with_results_factor_level_contrast(
            &counts,
            &full_design,
            FactorLevelContrast::new("condition", "B", "A", &levels),
        )
        .unwrap();
    let (expected_fit, expected_results) = builder
        .fit_lrt_glm_mu_contrast(&counts, &full_design, &reduced_design, &[0.0, 1.0])
        .unwrap();

    assert_eq!(numeric_fit.lrt, expected_fit.lrt);
    assert_eq!(spec_fit.lrt, expected_fit.lrt);
    assert_eq!(factor_fit.lrt, expected_fit.lrt);
    assert_eq!(numeric_results, expected_results);
    assert_eq!(
        spec_results.metadata.comparison.as_deref(),
        Some("coefficient condition_B_vs_A")
    );
    assert_eq!(
        factor_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
}

#[test]
fn fixed_lrt_model_frame_contrast_retains_formula_derived_metadata() {
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
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; counts.n_samples()])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .model_frame(model_frame);
    let formula_design = builder
        .expanded_formula_design_with_offsets("~ factor(condition, levels=c('A','B'))")
        .unwrap();
    let reduced_design = DesignMatrix::intercept_only(counts.n_samples()).unwrap();
    let dispersions = vec![0.1; counts.n_genes()];
    let contrast = ResultsContrast::character("factor.condition..levels...c..A....B...", "B", "A");

    let (fit, mut results) = builder
        .fit_fixed_dispersion_lrt_results_contrast_from_model_frame(
            &counts,
            &formula_design.design.standard_design,
            &reduced_design,
            &dispersions,
            &contrast,
            &formula_design.model_frame,
        )
        .unwrap();
    results.independent_filtering = None;
    let object_results = fit
        .results_contrast_request(&counts, &contrast, counts.gene_names())
        .unwrap();

    assert_eq!(object_results, results);
    assert_eq!(
        fit.current_model_frame()
            .unwrap()
            .factors
            .iter()
            .map(|factor| factor.name.as_str())
            .collect::<Vec<_>>(),
        vec!["condition", "factor(condition, levels = c(\"A\", \"B\"))"]
    );
    assert_eq!(results.metadata.test_type, Some(TestType::Lrt));
}

#[test]
fn top_level_fit_with_results_cooks_replacement_keeps_lrt_explicit() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let err = glm_mu_native_wald_builder()
        .test(TestType::Lrt)
        .fit_with_results_with_cooks_replacement(
            &counts,
            &design,
            &CooksReplacementOptions::new(f64::MAX),
        )
        .unwrap_err();
    let contrast_err = glm_mu_native_wald_builder()
        .test(TestType::Lrt)
        .fit_with_results_contrast_with_cooks_replacement(
            &counts,
            &design,
            &[0.0, 1.0],
            &CooksReplacementOptions::new(f64::MAX),
        )
        .unwrap_err();
    let spec_err = glm_mu_native_wald_builder()
        .test(TestType::Lrt)
        .fit_with_results_contrast_spec_with_cooks_replacement(
            &counts,
            &design,
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
            &CooksReplacementOptions::new(f64::MAX),
        )
        .unwrap_err();
    let factor_err = glm_mu_native_wald_builder()
        .test(TestType::Lrt)
        .fit_with_results_factor_level_contrast_with_cooks_replacement(
            &counts,
            &design,
            FactorLevelContrast::new("condition", "B", "A", &levels),
            &CooksReplacementOptions::new(f64::MAX),
        )
        .unwrap_err();

    assert!(matches!(err, DeseqError::UnsupportedFeature { .. }));
    assert!(matches!(
        contrast_err,
        DeseqError::UnsupportedFeature { .. }
    ));
    assert!(matches!(spec_err, DeseqError::UnsupportedFeature { .. }));
    assert!(matches!(factor_err, DeseqError::UnsupportedFeature { .. }));
}

#[test]
fn top_level_test_output_routes_lrt_cooks_replacement_with_stored_reduced_design() {
    let counts = native_wald_counts_with_zero_row();
    let full_design = two_group_design();
    let reduced_design = DesignMatrix::intercept_only(counts.n_samples()).unwrap();
    let options = CooksReplacementOptions::new(f64::MAX);
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

    let top_level = builder
        .fit_with_test_results_with_cooks_replacement(&counts, &full_design, &options)
        .unwrap();
    let formula_top_level = builder
        .fit_formula_with_test_results_with_cooks_replacement(&counts, "~ condition", &options)
        .unwrap();
    let formula_named_top_level = builder
        .fit_formula_with_test_results_name_with_cooks_replacement(
            &counts,
            "~ condition",
            "condition_B_vs_A",
            &options,
        )
        .unwrap();
    let expected = builder
        .fit_lrt_with_results_with_cooks_replacement(
            &counts,
            &full_design,
            &reduced_design,
            &options,
        )
        .unwrap();

    assert_eq!(top_level.test_type(), TestType::Lrt);
    assert_eq!(top_level.results(), &expected.results);
    assert_eq!(top_level.original_results(), &expected.original_results);
    assert_eq!(top_level.refit_plan(), &expected.refit_plan);
    assert_eq!(
        top_level.original_fit().current_model_frame(),
        expected.original_fit.current_model_frame()
    );
    assert_eq!(
        top_level.refit_fit().is_some(),
        expected.refit_fit.is_some()
    );
    assert_eq!(top_level.refit_results(), expected.refit_results.as_ref());

    let CooksReplacementTestOutput::Lrt(output) = top_level else {
        panic!("top-level test output should route to LRT");
    };
    let CooksReplacementTestOutput::Lrt(formula_output) = formula_top_level else {
        panic!("formula top-level test output should route to LRT");
    };
    let CooksReplacementTestOutput::Lrt(formula_named_output) = formula_named_top_level else {
        panic!("formula named top-level test output should route to LRT");
    };
    assert_eq!(output.results, expected.results);
    assert_eq!(formula_output.results, expected.results);
    assert_eq!(formula_named_output.results, expected.results);
    assert_eq!(output.refit_plan.n_refit, expected.refit_plan.n_refit);
    assert_eq!(
        formula_output.refit_plan.n_refit,
        expected.refit_plan.n_refit
    );
    assert_eq!(
        formula_named_output.refit_plan.n_refit,
        expected.refit_plan.n_refit
    );
    assert!(output.original_fit.lrt.is_some());
    assert!(formula_output.original_fit.lrt.is_some());
    assert!(formula_named_output.original_fit.lrt.is_some());
    assert_eq!(output.results.metadata.test_type, Some(TestType::Lrt));
}
