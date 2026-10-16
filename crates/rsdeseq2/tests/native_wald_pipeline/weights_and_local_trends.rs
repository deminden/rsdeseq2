use super::*;

#[test]
fn native_glm_mu_unit_weights_match_unweighted_wald() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);
    let unit_weights = unit_weights_for(&counts);

    let (unweighted_fit, unweighted_results) = builder
        .clone()
        .fit_wald_glm_mu_parametric(&counts, &design, 1)
        .unwrap();
    let (weighted_fit, weighted_results) = builder
        .observation_weights(unit_weights.clone())
        .fit_wald_glm_mu_parametric(&counts, &design, 1)
        .unwrap();

    let expected_weights_fail = vec![false; counts.n_genes()];
    assert_eq!(
        weighted_fit.observation_weights.as_ref(),
        Some(&unit_weights)
    );
    assert_eq!(
        weighted_fit.weights_fail.as_ref(),
        Some(&expected_weights_fail)
    );
    assert_eq!(weighted_fit.all_zero, unweighted_fit.all_zero);
    assert_eq!(weighted_results.rows.len(), unweighted_results.rows.len());

    assert_slice_close_or_nan(
        &weighted_fit.size_factors,
        &unweighted_fit.size_factors,
        "size_factors",
    );
    assert_slice_close_or_nan(
        &weighted_fit.base_mean,
        &unweighted_fit.base_mean,
        "base_mean",
    );
    assert_slice_close_or_nan(&weighted_fit.base_var, &unweighted_fit.base_var, "base_var");
    assert_slice_close_or_nan(
        weighted_fit.disp_gene_est.as_ref().unwrap(),
        unweighted_fit.disp_gene_est.as_ref().unwrap(),
        "disp_gene_est",
    );
    assert_slice_close_or_nan(
        weighted_fit.disp_fit.as_ref().unwrap(),
        unweighted_fit.disp_fit.as_ref().unwrap(),
        "disp_fit",
    );
    assert_slice_close_or_nan(
        weighted_fit.disp_map.as_ref().unwrap(),
        unweighted_fit.disp_map.as_ref().unwrap(),
        "disp_map",
    );
    assert_slice_close_or_nan(
        weighted_fit.dispersion.as_ref().unwrap(),
        unweighted_fit.dispersion.as_ref().unwrap(),
        "dispersion",
    );
    assert_matrix_close_or_nan(
        weighted_fit.beta.as_ref().unwrap(),
        unweighted_fit.beta.as_ref().unwrap(),
        "beta",
    );
    assert_matrix_close_or_nan(
        weighted_fit.beta_se.as_ref().unwrap(),
        unweighted_fit.beta_se.as_ref().unwrap(),
        "beta_se",
    );
    assert_matrix_close_or_nan(
        weighted_fit.beta_covariance.as_ref().unwrap(),
        unweighted_fit.beta_covariance.as_ref().unwrap(),
        "beta_covariance",
    );
    assert_option_slice_close(
        &weighted_fit.wald.as_ref().unwrap().pvalue,
        &unweighted_fit.wald.as_ref().unwrap().pvalue,
        "wald_pvalue",
    );

    for (gene, (weighted, unweighted)) in weighted_results
        .rows
        .iter()
        .zip(unweighted_results.rows.iter())
        .enumerate()
    {
        assert_eq!(weighted.gene, unweighted.gene, "gene {gene}");
        assert_float_close_or_nan(weighted.base_mean, unweighted.base_mean, "result baseMean");
        assert_option_close(
            weighted.log2_fold_change,
            unweighted.log2_fold_change,
            "result log2FoldChange",
        );
        assert_option_close(weighted.lfc_se, unweighted.lfc_se, "result lfcSE");
        assert_option_close(weighted.stat, unweighted.stat, "result stat");
        assert_option_close(weighted.pvalue, unweighted.pvalue, "result pvalue");
        assert_option_close(weighted.padj, unweighted.padj, "result padj");
        assert_option_close(
            weighted.dispersion,
            unweighted.dispersion,
            "result dispersion",
        );
        assert_eq!(weighted.converged, unweighted.converged, "gene {gene}");
    }
}

#[test]
fn native_glm_mu_observation_weights_run_through_wald() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder();
    let weights = nonunit_weights_for(&counts);

    let (unweighted_fit, _) = builder
        .clone()
        .fit_wald_glm_mu_parametric(&counts, &design, 1)
        .unwrap();
    let (fit, results) = builder
        .observation_weights(weights)
        .fit_wald_glm_mu_parametric(&counts, &design, 1)
        .unwrap();

    let expected_weights_fail = vec![false; counts.n_genes()];
    assert_eq!(fit.weights_fail.as_ref(), Some(&expected_weights_fail));
    assert!(fit.observation_weights.is_some());
    assert_eq!(
        fit.all_zero,
        vec![true, false, false, false, false, false, false, false]
    );
    assert!(fit.base_mean[2] < unweighted_fit.base_mean[2]);
    assert_eq!(results.rows.len(), counts.n_genes());
    assert_eq!(results.rows[0].pvalue, None);

    let disp_gene_est = fit.disp_gene_est.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    assert!(disp_gene_est[0].is_nan());
    assert!(dispersion[0].is_nan());
    for gene in 1..counts.n_genes() {
        assert!(disp_gene_est[gene].is_finite());
        assert!(dispersion[gene].is_finite());
        assert!(dispersion[gene] > 0.0);
    }
    assert!(
        fit.wald
            .as_ref()
            .unwrap()
            .pvalue
            .iter()
            .skip(1)
            .any(|pvalue| pvalue.is_some_and(f64::is_finite))
    );
    assert!(
        results
            .rows
            .iter()
            .skip(1)
            .any(|row| row.pvalue.is_some_and(f64::is_finite))
    );
}

#[test]
fn native_glm_mu_generic_mean_wald_runs_through_map_and_glm() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = native_wald_builder()
        .fit_type(FitType::Mean)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        });

    let (fit, results) = builder.fit_wald_glm_mu(&counts, &design, 1).unwrap();

    assert_eq!(results.rows.len(), counts.n_genes());
    assert!(fit.disp_prior_var.unwrap().is_finite());
    assert_eq!(
        fit.all_zero,
        vec![true, false, false, false, false, false, false, false]
    );
    assert!(fit.dispersion.as_ref().unwrap()[0].is_nan());
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
    assert_eq!(results.rows[0].pvalue, None);
    assert!(results.rows[1].pvalue.unwrap().is_finite());
}

#[test]
fn native_glm_mu_cooks_replacement_refit_merges_refit_rows() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder();

    let output = builder
        .fit_wald_glm_mu_with_cooks_replacement(
            &counts,
            &design,
            1,
            &CooksReplacementOptions {
                trim: 0.2,
                cooks_cutoff: 0.0,
                min_replicates: 3,
                which_samples: Some(vec![true, false, false, false, false, false, false, false]),
            },
        )
        .unwrap();

    assert!(output.refit_plan.n_refit > 0);
    assert!(output.refit_plan.should_refit);
    assert!(!output.refit_plan.refit_rows.is_empty());
    assert!(output.refit_fit.is_some());
    assert!(output.refit_results.is_some());
    assert_ne!(
        output.refit_plan.replacement.replaced_counts.as_slice(),
        counts.as_slice()
    );
    assert_eq!(
        output.refit_plan.replacement.replaceable_samples,
        vec![true, false, false, false, false, false, false, false]
    );

    let refit_results = output.refit_results.as_ref().unwrap();
    for gene in output.refit_plan.refit_rows.iter().copied() {
        assert_eq!(
            output.results.rows[gene].log2_fold_change,
            refit_results.rows[gene].log2_fold_change
        );
        assert_eq!(
            output.results.rows[gene].pvalue,
            refit_results.rows[gene].pvalue
        );
        assert_eq!(
            output.results.rows[gene].max_cooks,
            output.refit_plan.post_refit_max_cooks[gene]
        );
        assert_relative_eq!(
            output.results.rows[gene].base_mean,
            output.refit_plan.replaced_base_mean[gene],
            epsilon = 1e-12
        );
    }
    for gene in output.refit_plan.new_all_zero_rows.iter().copied() {
        assert_eq!(output.results.rows[gene].pvalue, Some(1.0));
        assert_eq!(output.results.rows[gene].log2_fold_change, Some(0.0));
        assert_eq!(output.results.rows[gene].lfc_se, Some(0.0));
        assert_eq!(output.results.rows[gene].stat, Some(0.0));
        assert_eq!(output.results.rows[gene].dispersion, None);
    }
}

#[test]
fn native_glm_mu_contrast_cooks_replacement_refit_merges_refit_rows() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder();
    let contrast = [0.0, 1.0];

    let output = builder
        .fit_wald_glm_mu_contrast_with_cooks_replacement(
            &counts,
            &design,
            &contrast,
            &CooksReplacementOptions {
                trim: 0.2,
                cooks_cutoff: 0.0,
                min_replicates: 3,
                which_samples: Some(vec![true, false, false, false, false, false, false, false]),
            },
        )
        .unwrap();

    assert!(output.refit_plan.n_refit > 0);
    assert!(output.refit_plan.should_refit);
    assert!(output.refit_fit.is_some());
    assert!(output.refit_results.is_some());
    assert_ne!(
        output.refit_plan.replacement.replaced_counts.as_slice(),
        counts.as_slice()
    );
    assert_eq!(
        output.results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        output.results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );

    let refit_results = output.refit_results.as_ref().unwrap();
    for gene in output.refit_plan.refit_rows.iter().copied() {
        assert_eq!(
            output.results.rows[gene].log2_fold_change,
            refit_results.rows[gene].log2_fold_change
        );
        assert_eq!(
            output.results.rows[gene].pvalue,
            refit_results.rows[gene].pvalue
        );
        assert_eq!(
            output.results.rows[gene].max_cooks,
            output.refit_plan.post_refit_max_cooks[gene]
        );
        assert_relative_eq!(
            output.results.rows[gene].base_mean,
            output.refit_plan.replaced_base_mean[gene],
            epsilon = 1e-12
        );
    }
    for gene in output.refit_plan.new_all_zero_rows.iter().copied() {
        assert_eq!(output.results.rows[gene].pvalue, None);
        assert_eq!(output.results.rows[gene].log2_fold_change, None);
        assert_eq!(output.results.rows[gene].dispersion, None);
    }
}

#[test]
fn native_glm_mu_contrast_spec_cooks_replacement_preserves_metadata() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder();

    let (expected_fit, expected_results) = builder
        .clone()
        .fit_wald_glm_mu_contrast(&counts, &design, &[0.0, 1.0])
        .unwrap();
    let output = builder
        .fit_wald_glm_mu_contrast_spec_with_cooks_replacement(
            &counts,
            &design,
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
            &CooksReplacementOptions {
                trim: 0.2,
                cooks_cutoff: f64::MAX,
                min_replicates: 3,
                which_samples: None,
            },
        )
        .unwrap();

    assert_eq!(output.refit_plan.n_refit, 0);
    assert!(output.refit_fit.is_none());
    assert!(output.refit_results.is_none());
    assert_wald_fit_state_matches(
        &output.original_fit,
        &expected_fit,
        "native GLM-mu named contrast no-refit original",
    );
    assert_eq!(output.original_results.rows, expected_results.rows);
    assert_eq!(output.results.rows, expected_results.rows);
    assert_eq!(
        output.original_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        output.original_results.metadata.comparison.as_deref(),
        Some("coefficient condition_B_vs_A")
    );
    assert_eq!(
        output.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        output.results.metadata.comparison.as_deref(),
        Some("coefficient condition_B_vs_A")
    );
}

#[test]
fn native_glm_mu_factor_level_contrast_cooks_replacement_preserves_metadata() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

    let (expected_fit, expected_results) = builder
        .clone()
        .fit_wald_glm_mu_contrast(&counts, &design, &[0.0, 1.0])
        .unwrap();
    let output = builder
        .fit_wald_glm_mu_factor_level_contrast_with_cooks_replacement(
            &counts,
            &design,
            FactorLevelContrast::new("condition", "B", "A", &levels),
            &CooksReplacementOptions {
                trim: 0.2,
                cooks_cutoff: f64::MAX,
                min_replicates: 3,
                which_samples: None,
            },
        )
        .unwrap();

    assert_eq!(output.refit_plan.n_refit, 0);
    assert!(output.refit_fit.is_none());
    assert!(output.refit_results.is_none());
    assert_wald_fit_state_matches(
        &output.original_fit,
        &expected_fit,
        "native GLM-mu factor-level contrast no-refit original",
    );
    assert_eq!(output.original_results.rows, expected_results.rows);
    assert_eq!(output.results.rows, expected_results.rows);
    assert_eq!(
        output.original_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        output.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        output.results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
}

#[test]
fn native_glm_mu_cooks_replacement_refit_skips_when_no_rows_are_marked() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder();

    let output = builder
        .fit_wald_glm_mu_with_cooks_replacement(
            &counts,
            &design,
            1,
            &CooksReplacementOptions {
                trim: 0.2,
                cooks_cutoff: f64::MAX,
                min_replicates: 3,
                which_samples: None,
            },
        )
        .unwrap();

    assert_eq!(output.refit_plan.n_refit, 0);
    assert!(!output.refit_plan.should_refit);
    assert!(output.refit_fit.is_none());
    assert!(output.refit_results.is_none());
    assert_eq!(output.results, output.original_results);
    assert_eq!(
        output.refit_plan.replacement.replaced_counts.as_slice(),
        counts.as_slice()
    );
}

#[test]
fn native_linear_mu_generic_mean_wald_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = native_wald_builder().fit_type(FitType::Mean);

    let (native_fit, native_results) = builder.fit_wald_linear_mu(&counts, &design, 1).unwrap();
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
fn native_linear_mu_local_wald_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = native_wald_builder().fit_type(FitType::Local);

    let (native_fit, native_results) = builder.fit_wald_linear_mu(&counts, &design, 1).unwrap();
    let final_dispersions = native_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_fit, fixed_results) = builder
        .fit_fixed_dispersion_wald(&counts, &design, &final_dispersions, 1)
        .unwrap();

    assert!(native_fit.disp_prior_var.unwrap().is_finite());
    assert_eq!(
        native_fit.disp_fit.as_ref().unwrap().len(),
        counts.n_genes()
    );

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
fn native_glm_mu_local_wald_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Local);

    let (native_fit, native_results) = builder.fit_wald_glm_mu(&counts, &design, 1).unwrap();
    let final_dispersions = native_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_fit, fixed_results) = builder
        .fit_fixed_dispersion_wald(&counts, &design, &final_dispersions, 1)
        .unwrap();

    assert!(native_fit.disp_prior_var.unwrap().is_finite());
    assert_eq!(
        native_fit.disp_fit.as_ref().unwrap().len(),
        counts.n_genes()
    );

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
fn native_linear_mu_generic_pipeline_accepts_local_and_rejects_glm_gam_poi() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();

    let local_fit = native_wald_builder()
        .fit_type(FitType::Local)
        .fit_map_dispersions_linear_mu(&counts, &design)
        .unwrap();
    assert_eq!(
        local_fit.dispersion.as_ref().unwrap().len(),
        counts.n_genes()
    );

    let glm_gam_poi_err = native_wald_builder()
        .fit_type(FitType::GlmGamPoi)
        .fit_wald_linear_mu(&counts, &design, 1)
        .unwrap_err();
    assert!(
        glm_gam_poi_err
            .to_string()
            .contains("glmGamPoi dispersion trend")
    );
}
