use super::*;

#[test]
fn native_wald_parametric_contrast_helpers_ignore_builder_fit_type() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let mean_linear_builder = native_wald_builder().fit_type(FitType::Mean);
    let parametric_linear_builder = mean_linear_builder.clone().fit_type(FitType::Parametric);
    let mean_glm_builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);
    let parametric_glm_builder = mean_glm_builder.clone().fit_type(FitType::Parametric);
    let spec = ContrastSpec::coefficient_name("condition_B_vs_A");
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let factor_contrast = FactorLevelContrast::new("condition", "B", "A", &levels);

    let (linear_parametric_fit, linear_parametric_results) = mean_linear_builder
        .fit_wald_linear_mu_contrast_parametric(&counts, &design, &[0.0, 1.0])
        .unwrap();
    let (linear_expected_fit, linear_expected_results) = parametric_linear_builder
        .fit_wald_linear_mu_contrast(&counts, &design, &[0.0, 1.0])
        .unwrap();
    assert_wald_fit_state_matches(
        &linear_parametric_fit,
        &linear_expected_fit,
        "linear-mu parametric Wald contrast",
    );
    assert_eq!(linear_parametric_results, linear_expected_results);

    let (linear_named_fit, linear_named_results) = mean_linear_builder
        .fit_wald_linear_mu_contrast_spec_parametric(&counts, &design, &spec)
        .unwrap();
    assert_wald_fit_state_matches(
        &linear_named_fit,
        &linear_expected_fit,
        "linear-mu named parametric Wald contrast",
    );
    assert_eq!(linear_named_results.rows, linear_expected_results.rows);
    assert_eq!(
        linear_named_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );

    let (linear_factor_fit, linear_factor_results) = mean_linear_builder
        .fit_wald_linear_mu_factor_level_contrast_parametric(&counts, &design, factor_contrast)
        .unwrap();
    assert_wald_fit_state_matches(
        &linear_factor_fit,
        &linear_expected_fit,
        "linear-mu factor parametric Wald contrast",
    );
    assert_eq!(linear_factor_results.rows, linear_expected_results.rows);
    assert_eq!(
        linear_factor_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );

    let (glm_parametric_fit, glm_parametric_results) = mean_glm_builder
        .fit_wald_glm_mu_contrast_parametric(&counts, &design, &[0.0, 1.0])
        .unwrap();
    let (glm_expected_fit, glm_expected_results) = parametric_glm_builder
        .fit_wald_glm_mu_contrast(&counts, &design, &[0.0, 1.0])
        .unwrap();
    assert_wald_fit_state_matches(
        &glm_parametric_fit,
        &glm_expected_fit,
        "GLM-mu parametric Wald contrast",
    );
    assert_eq!(glm_parametric_results, glm_expected_results);

    let (glm_named_fit, glm_named_results) = mean_glm_builder
        .fit_wald_glm_mu_contrast_spec_parametric(&counts, &design, &spec)
        .unwrap();
    assert_wald_fit_state_matches(
        &glm_named_fit,
        &glm_expected_fit,
        "GLM-mu named parametric Wald contrast",
    );
    assert_eq!(glm_named_results.rows, glm_expected_results.rows);
    assert_eq!(
        glm_named_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );

    let (glm_factor_fit, glm_factor_results) = mean_glm_builder
        .fit_wald_glm_mu_factor_level_contrast_parametric(&counts, &design, factor_contrast)
        .unwrap();
    assert_wald_fit_state_matches(
        &glm_factor_fit,
        &glm_expected_fit,
        "GLM-mu factor parametric Wald contrast",
    );
    assert_eq!(glm_factor_results.rows, glm_expected_results.rows);
    assert_eq!(
        glm_factor_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
}

#[test]
fn native_linear_mu_generic_mean_wald_runs_through_map_and_glm() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = native_wald_builder().fit_type(FitType::Mean);

    let (fit, results) = builder.fit_wald_linear_mu(&counts, &design, 1).unwrap();

    assert_eq!(results.rows.len(), counts.n_genes());
    assert_eq!(
        fit.all_zero,
        vec![true, false, false, false, false, false, false, false]
    );
    assert!(fit.disp_prior_var.unwrap().is_finite());

    let disp_fit = fit.disp_fit.as_ref().unwrap();
    assert!(disp_fit[0].is_nan());
    let first_non_zero_fit = disp_fit[1];
    assert!(first_non_zero_fit.is_finite());
    for value in disp_fit.iter().copied().skip(1) {
        assert_relative_eq!(value, first_non_zero_fit, epsilon = 1e-12);
    }

    let dispersion = fit.dispersion.as_ref().unwrap();
    assert!(dispersion[0].is_nan());
    for value in dispersion.iter().copied().skip(1) {
        assert!(value.is_finite());
        assert!(value > 0.0);
    }

    assert_eq!(fit.beta.as_ref().unwrap().n_cols(), design.n_coefficients());
    assert_eq!(
        fit.beta_covariance.as_ref().unwrap().n_rows(),
        counts.n_genes()
    );
    assert_eq!(
        fit.beta_covariance.as_ref().unwrap().n_cols(),
        design.n_coefficients() * design.n_coefficients()
    );
    assert_eq!(fit.cooks.as_ref().unwrap().n_rows(), counts.n_genes());
    assert_eq!(fit.wald.as_ref().unwrap().pvalue.len(), counts.n_genes());
    assert_eq!(fit.beta_iter.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.log_like.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(results.rows[0].pvalue, None);
    assert!(results.rows[1].pvalue.unwrap().is_finite());
}

#[test]
fn native_glm_mu_parametric_map_preserves_dispersion_intermediates() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();

    let fit = native_wald_builder()
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_map_dispersions_glm_mu_parametric(&counts, &design)
        .unwrap();

    assert_eq!(fit.design.as_ref().unwrap(), &design);
    assert_eq!(
        fit.all_zero,
        vec![true, false, false, false, false, false, false, false]
    );

    let disp_gene_est = fit.disp_gene_est.as_ref().unwrap();
    let disp_gene_iter = fit.disp_gene_iter.as_ref().unwrap();
    let disp_fit = fit.disp_fit.as_ref().unwrap();
    let disp_map = fit.disp_map.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    let disp_iter = fit.disp_iter.as_ref().unwrap();

    assert_eq!(disp_gene_est.len(), counts.n_genes());
    assert_eq!(disp_gene_iter.len(), counts.n_genes());
    assert_eq!(disp_fit.len(), counts.n_genes());
    assert_eq!(disp_map.len(), counts.n_genes());
    assert_eq!(dispersion.len(), counts.n_genes());
    assert_eq!(disp_iter.len(), counts.n_genes());
    assert!(fit.disp_prior_var.unwrap().is_finite());
    assert!(fit.beta.is_none());
    assert!(fit.wald.is_none());

    assert!(disp_gene_est[0].is_nan());
    assert_eq!(disp_gene_iter[0], 0);
    assert!(disp_fit[0].is_nan());
    assert!(disp_map[0].is_nan());
    assert!(dispersion[0].is_nan());

    for gene in 1..counts.n_genes() {
        assert!(disp_gene_est[gene].is_finite());
        assert!(disp_gene_iter[gene] > 0);
        assert!(disp_fit[gene].is_finite());
        assert!(disp_map[gene].is_finite());
        assert!(dispersion[gene].is_finite());
        assert!(dispersion[gene] > 0.0);
    }
    assert_eq!(fit.mu.as_ref().unwrap().n_cols(), counts.n_samples());
}

#[test]
fn native_glm_mu_generic_mean_map_runs_through_selected_trend() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();

    let fit = native_wald_builder()
        .fit_type(FitType::Mean)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_map_dispersions_glm_mu(&counts, &design)
        .unwrap();

    let disp_fit = fit.disp_fit.as_ref().unwrap();
    assert!(disp_fit[0].is_nan());
    let first_non_zero_fit = disp_fit[1];
    assert!(first_non_zero_fit.is_finite());
    for value in disp_fit.iter().copied().skip(1) {
        assert_relative_eq!(value, first_non_zero_fit, epsilon = 1e-12);
    }

    let dispersion = fit.dispersion.as_ref().unwrap();
    assert!(dispersion[0].is_nan());
    for value in dispersion.iter().copied().skip(1) {
        assert!(value.is_finite());
        assert!(value > 0.0);
    }
}

#[test]
fn native_glm_mu_parametric_wald_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);

    let (native_fit, native_results) = builder
        .fit_wald_glm_mu_parametric(&counts, &design, 1)
        .unwrap();
    let final_dispersions = native_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_fit, fixed_results) = builder
        .fit_fixed_dispersion_wald(&counts, &design, &final_dispersions, 1)
        .unwrap();

    assert!(native_fit.disp_prior_var.unwrap().is_finite());
    assert!(native_fit.beta.as_ref().unwrap().n_cols() == design.n_coefficients());
    assert!(native_fit.wald.is_some());
    assert_eq!(native_results.rows.len(), counts.n_genes());
    assert_eq!(native_results.rows[0].pvalue, None);

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
fn native_glm_mu_wald_contrast_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);
    let contrast = [0.0, 1.0];

    let (native_fit, native_results) = builder
        .fit_wald_glm_mu_contrast(&counts, &design, &contrast)
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
fn native_glm_mu_wald_contrast_specs_set_metadata_and_factor_levels() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

    let (primitive_fit, _primitive_results) = builder
        .clone()
        .fit_wald_glm_mu_contrast(&counts, &design, &[0.0, 1.0])
        .unwrap();
    let (named_fit, named_results) = builder
        .clone()
        .fit_wald_glm_mu_contrast_spec(
            &counts,
            &design,
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
        )
        .unwrap();
    let (factor_fit, factor_results) = builder
        .fit_wald_glm_mu_factor_level_contrast(
            &counts,
            &design,
            FactorLevelContrast::new("condition", "B", "A", &levels),
        )
        .unwrap();

    assert_wald_fit_state_matches(
        &named_fit,
        &primitive_fit,
        "native GLM-mu named contrast spec",
    );
    assert_wald_fit_state_matches(
        &factor_fit,
        &primitive_fit,
        "native GLM-mu factor-level contrast",
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
