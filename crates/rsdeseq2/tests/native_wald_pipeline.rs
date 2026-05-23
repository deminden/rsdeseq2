use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

fn two_group_design() -> DesignMatrix {
    DesignMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap()
}

fn native_wald_counts_with_zero_row() -> CountMatrix {
    CountMatrix::from_row_major_u32_with_names(
        8,
        8,
        vec![
            0, 0, 0, 0, 0, 0, 0, 0, //
            0, 20, 1, 19, 2, 18, 3, 17, //
            12, 28, 10, 30, 14, 26, 11, 29, //
            30, 50, 25, 55, 35, 45, 28, 52, //
            55, 105, 60, 100, 50, 110, 65, 95, //
            120, 200, 130, 190, 115, 205, 125, 195, //
            240, 400, 250, 390, 230, 410, 260, 380, //
            15, 18, 12, 17, 45, 50, 40, 55,
        ],
        Some(
            [
                "zero", "up", "flat", "variable", "high_up", "stable", "low_up", "broad",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
        ),
        None,
    )
    .unwrap()
}

fn native_wald_builder() -> DeseqBuilder {
    DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            ..GeneWiseDispersionOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
}

fn glm_mu_native_wald_builder() -> DeseqBuilder {
    native_wald_builder().gene_wise_dispersion_options(GeneWiseDispersionOptions {
        use_cox_reid: false,
        fit_method: GeneWiseDispersionFitMethod::Grid,
        niter: 2,
        ..GeneWiseDispersionOptions::default()
    })
}

fn unit_weights_for(counts: &CountMatrix) -> RowMajorMatrix<f64> {
    RowMajorMatrix::from_elem(counts.n_genes(), counts.n_samples(), 1.0).unwrap()
}

fn nonunit_weights_for(counts: &CountMatrix) -> RowMajorMatrix<f64> {
    let sample_pattern = [1.0, 0.85, 1.0, 0.8, 1.0, 0.75, 1.0, 0.7];
    let mut values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for gene in 0..counts.n_genes() {
        let gene_scale = 1.0 - 0.03 * (gene % 3) as f64;
        for weight in sample_pattern.iter().take(counts.n_samples()) {
            values.push(*weight * gene_scale);
        }
    }
    RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values).unwrap()
}

fn assert_float_close_or_nan(actual: f64, expected: f64, label: &str) {
    if expected.is_nan() {
        assert!(actual.is_nan(), "{label}: expected NaN, got {actual}");
        return;
    }
    let diff = (actual - expected).abs();
    let allowed = 1e-9 + 1e-9 * expected.abs().max(1.0);
    assert!(
        diff <= allowed,
        "{label}: actual={actual}, expected={expected}, diff={diff}, allowed={allowed}"
    );
}

fn assert_slice_close_or_nan(actual: &[f64], expected: &[f64], label: &str) {
    assert_eq!(actual.len(), expected.len(), "{label}: length mismatch");
    for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
        assert_float_close_or_nan(*actual, *expected, &format!("{label}[{index}]"));
    }
}

fn assert_matrix_close_or_nan(
    actual: &RowMajorMatrix<f64>,
    expected: &RowMajorMatrix<f64>,
    label: &str,
) {
    assert_eq!(actual.n_rows(), expected.n_rows(), "{label}: row mismatch");
    assert_eq!(
        actual.n_cols(),
        expected.n_cols(),
        "{label}: column mismatch"
    );
    assert_slice_close_or_nan(actual.as_slice(), expected.as_slice(), label);
}

fn assert_option_close(actual: Option<f64>, expected: Option<f64>, label: &str) {
    match (actual, expected) {
        (Some(actual), Some(expected)) => assert_float_close_or_nan(actual, expected, label),
        (None, None) => {}
        _ => panic!("{label}: actual={actual:?}, expected={expected:?}"),
    }
}

fn assert_option_slice_close(actual: &[Option<f64>], expected: &[Option<f64>], label: &str) {
    assert_eq!(actual.len(), expected.len(), "{label}: length mismatch");
    for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
        assert_option_close(*actual, *expected, &format!("{label}[{index}]"));
    }
}

#[test]
fn native_linear_mu_parametric_wald_preserves_dispersion_intermediates() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();

    let (fit, results) = native_wald_builder()
        .fit_wald_linear_mu_parametric(&counts, &design, 1)
        .unwrap();

    assert_eq!(fit.design.as_ref().unwrap(), &design);
    assert_eq!(results.rows.len(), counts.n_genes());
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
    let disp_outlier = fit.disp_outlier.as_ref().unwrap();
    let disp_converged = fit.dispersion_converged.as_ref().unwrap();

    assert_eq!(disp_gene_est.len(), counts.n_genes());
    assert_eq!(disp_gene_iter.len(), counts.n_genes());
    assert_eq!(disp_fit.len(), counts.n_genes());
    assert_eq!(disp_map.len(), counts.n_genes());
    assert_eq!(dispersion.len(), counts.n_genes());
    assert_eq!(disp_iter.len(), counts.n_genes());
    assert_eq!(disp_outlier.len(), counts.n_genes());
    assert_eq!(disp_converged.len(), counts.n_genes());
    assert!(fit.disp_prior_var.unwrap().is_finite());

    assert!(disp_gene_est[0].is_nan());
    assert_eq!(disp_gene_iter[0], 0);
    assert!(disp_fit[0].is_nan());
    assert!(disp_map[0].is_nan());
    assert!(dispersion[0].is_nan());
    assert_eq!(disp_iter[0], 0);
    assert!(!disp_outlier[0]);
    assert!(!disp_converged[0]);

    for gene in 1..counts.n_genes() {
        assert!(disp_gene_est[gene].is_finite());
        assert!(disp_gene_iter[gene] > 0);
        assert!(disp_fit[gene].is_finite());
        assert!(disp_map[gene].is_finite());
        assert!(dispersion[gene].is_finite());
        assert!(dispersion[gene] > 0.0);
    }

    assert_eq!(fit.beta.as_ref().unwrap().n_cols(), design.n_coefficients());
    assert_eq!(fit.beta_se.as_ref().unwrap().n_rows(), counts.n_genes());
    assert_eq!(
        fit.beta_covariance.as_ref().unwrap().n_rows(),
        counts.n_genes()
    );
    assert_eq!(
        fit.beta_covariance.as_ref().unwrap().n_cols(),
        design.n_coefficients() * design.n_coefficients()
    );
    assert!(fit.beta_covariance.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert_eq!(fit.mu.as_ref().unwrap().n_cols(), counts.n_samples());
    assert_eq!(
        fit.hat_diagonal.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert_eq!(fit.cooks.as_ref().unwrap().n_rows(), counts.n_genes());
    assert_eq!(fit.max_cooks.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.wald.as_ref().unwrap().pvalue.len(), counts.n_genes());

    assert_eq!(results.rows[0].gene.as_deref(), Some("zero"));
    assert_eq!(results.rows[0].log2_fold_change, None);
    assert_eq!(results.rows[0].pvalue, None);
    assert_eq!(results.rows[0].padj, None);
    assert_eq!(results.rows[0].dispersion, None);
    assert_eq!(results.rows[0].converged, None);

    assert_eq!(results.rows[1].gene.as_deref(), Some("up"));
    assert!(results.rows[1].log2_fold_change.unwrap().is_finite());
    assert!(results.rows[1].lfc_se.unwrap().is_finite());
    assert!(results.rows[1].stat.unwrap().is_finite());
    assert!(results.rows[1].pvalue.unwrap().is_finite());
    assert!(results.rows[1].padj.unwrap().is_finite());
    assert_eq!(results.rows[1].dispersion, Some(dispersion[1]));
    assert_eq!(results.rows[1].pvalue, fit.wald.as_ref().unwrap().pvalue[1]);
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
    assert!(fit
        .wald
        .as_ref()
        .unwrap()
        .pvalue
        .iter()
        .skip(1)
        .any(|pvalue| pvalue.is_some_and(f64::is_finite)));
    assert!(results
        .rows
        .iter()
        .skip(1)
        .any(|row| row.pvalue.is_some_and(f64::is_finite)));
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
        assert_eq!(output.results.rows[gene].pvalue, None);
        assert_eq!(output.results.rows[gene].log2_fold_change, None);
        assert_eq!(output.results.rows[gene].dispersion, None);
    }
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
fn native_linear_mu_generic_pipeline_rejects_unimplemented_fit_types() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();

    let local_err = native_wald_builder()
        .fit_type(FitType::Local)
        .fit_map_dispersions_linear_mu(&counts, &design)
        .unwrap_err();
    assert!(local_err.to_string().contains("local dispersion trend"));

    let glm_gam_poi_err = native_wald_builder()
        .fit_type(FitType::GlmGamPoi)
        .fit_wald_linear_mu(&counts, &design, 1)
        .unwrap_err();
    assert!(glm_gam_poi_err
        .to_string()
        .contains("glmGamPoi dispersion trend"));
}
