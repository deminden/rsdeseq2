use super::*;

#[test]
fn top_level_test_output_routes_wald_cooks_replacement() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
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
        .test(TestType::Wald)
        .model_frame(model_frame);

    let top_level = builder
        .fit_with_test_results_name_with_cooks_replacement(
            &counts,
            &design,
            "condition_B_vs_A",
            &options,
        )
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
        .fit_with_results_name_with_cooks_replacement(
            &counts,
            &design,
            "condition_B_vs_A",
            &options,
        )
        .unwrap();

    assert_eq!(top_level.test_type(), TestType::Wald);
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

    let CooksReplacementTestOutput::Wald(output) = top_level else {
        panic!("top-level test output should route to Wald");
    };
    let CooksReplacementTestOutput::Wald(formula_output) = formula_top_level else {
        panic!("formula top-level test output should route to Wald");
    };
    let CooksReplacementTestOutput::Wald(formula_named_output) = formula_named_top_level else {
        panic!("formula named top-level test output should route to Wald");
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
    assert!(output.original_fit.wald.is_some());
    assert!(formula_output.original_fit.wald.is_some());
    assert!(formula_named_output.original_fit.wald.is_some());
    assert_eq!(output.results.metadata.test_type, Some(TestType::Wald));
}

#[test]
fn native_linear_mu_parametric_wald_validates_coefficient_index() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();

    let err = native_wald_builder()
        .fit_wald_linear_mu_parametric(&counts, &design, 2)
        .unwrap_err();
    assert!(err.to_string().contains("Wald coefficient"));
}

#[test]
fn native_linear_mu_parametric_wald_accepts_normalization_factors() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let normalization_factors = RowMajorMatrix::from_row_major(
        counts.n_genes(),
        counts.n_samples(),
        vec![1.0; counts.n_genes() * counts.n_samples()],
    )
    .unwrap();

    let (fit, results) = native_wald_builder()
        .normalization_factors(normalization_factors.clone())
        .fit_wald_linear_mu_parametric(&counts, &design, 1)
        .unwrap();

    assert_eq!(fit.normalization_factors, Some(normalization_factors));
    assert_eq!(results.rows.len(), counts.n_genes());
    assert!(fit.dispersion.as_ref().unwrap()[1].is_finite());
    assert!(results.rows[1].pvalue.unwrap().is_finite());
}

#[test]
fn native_linear_mu_parametric_wald_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = native_wald_builder();

    let (native_fit, native_results) = builder
        .fit_wald_linear_mu_parametric(&counts, &design, 1)
        .unwrap();
    let final_dispersions = native_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_fit, fixed_results) = builder
        .fit_fixed_dispersion_wald(&counts, &design, &final_dispersions, 1)
        .unwrap();

    for gene in 0..counts.n_genes() {
        let native_disp = native_fit.dispersion.as_ref().unwrap()[gene];
        let fixed_disp = fixed_fit.dispersion.as_ref().unwrap()[gene];
        if native_disp.is_nan() {
            assert!(fixed_disp.is_nan());
        } else {
            assert_relative_eq!(native_disp, fixed_disp, epsilon = 1e-12);
        }
        assert_eq!(
            native_results.rows[gene].pvalue,
            fixed_results.rows[gene].pvalue
        );
        assert_eq!(
            native_results.rows[gene].padj,
            fixed_results.rows[gene].padj
        );
        assert_eq!(
            native_results.rows[gene].log2_fold_change,
            fixed_results.rows[gene].log2_fold_change
        );
    }
}

#[test]
fn native_linear_mu_wald_contrast_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = native_wald_builder().fit_type(FitType::Mean);
    let contrast = [0.0, 1.0];

    let (native_fit, native_results) = builder
        .fit_wald_linear_mu_contrast(&counts, &design, &contrast)
        .unwrap();
    let final_dispersions = native_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_fit, fixed_results) = builder
        .fit_fixed_dispersion_wald_contrast(&counts, &design, &final_dispersions, &contrast)
        .unwrap();

    assert!(native_fit.disp_prior_var.unwrap().is_finite());
    assert_eq!(
        native_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        native_results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );
    assert_eq!(native_fit.wald, fixed_fit.wald);

    for gene in 0..counts.n_genes() {
        let native_disp = native_fit.dispersion.as_ref().unwrap()[gene];
        let fixed_disp = fixed_fit.dispersion.as_ref().unwrap()[gene];
        if native_disp.is_nan() {
            assert!(fixed_disp.is_nan());
        } else {
            assert_relative_eq!(native_disp, fixed_disp, epsilon = 1e-12);
        }
        assert_eq!(
            native_results.rows[gene].pvalue,
            fixed_results.rows[gene].pvalue
        );
        assert_eq!(
            native_results.rows[gene].padj,
            fixed_results.rows[gene].padj
        );
        assert_eq!(
            native_results.rows[gene].log2_fold_change,
            fixed_results.rows[gene].log2_fold_change
        );
    }
}

#[test]
fn native_linear_mu_wald_contrast_specs_set_metadata_and_factor_levels() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = native_wald_builder().fit_type(FitType::Mean);
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

    let (primitive_fit, _primitive_results) = builder
        .clone()
        .fit_wald_linear_mu_contrast(&counts, &design, &[0.0, 1.0])
        .unwrap();
    let (named_fit, named_results) = builder
        .clone()
        .fit_wald_linear_mu_contrast_spec(
            &counts,
            &design,
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
        )
        .unwrap();
    let (factor_fit, factor_results) = builder
        .fit_wald_linear_mu_factor_level_contrast(
            &counts,
            &design,
            FactorLevelContrast::new("condition", "B", "A", &levels),
        )
        .unwrap();

    assert_wald_fit_state_matches(
        &named_fit,
        &primitive_fit,
        "native linear named contrast spec",
    );
    assert_wald_fit_state_matches(
        &factor_fit,
        &primitive_fit,
        "native linear factor-level contrast",
    );
    assert_eq!(
        named_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        named_results.metadata.comparison.as_deref(),
        Some("coefficient condition_B_vs_A")
    );
    assert_eq!(
        factor_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        factor_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
    for (named, factor) in named_results.rows.iter().zip(factor_results.rows.iter()) {
        assert_eq!(named.log2_fold_change, factor.log2_fold_change);
        assert_eq!(named.lfc_se, factor.lfc_se);
        assert_eq!(named.stat, factor.stat);
        assert_eq!(named.pvalue, factor.pvalue);
        assert_eq!(named.padj, factor.padj);
    }
}
