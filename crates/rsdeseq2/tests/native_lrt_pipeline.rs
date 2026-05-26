use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

fn full_design() -> DesignMatrix {
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

fn reduced_design() -> DesignMatrix {
    DesignMatrix::from_row_major(
        8,
        1,
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into()]),
    )
    .unwrap()
}

fn counts_with_zero_row() -> CountMatrix {
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

fn native_lrt_builder() -> DeseqBuilder {
    DeseqBuilder::new()
        .fit_type(FitType::Mean)
        .size_factors(vec![1.0; 8])
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
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
    for (index, (actual, expected)) in actual
        .as_slice()
        .iter()
        .zip(expected.as_slice().iter())
        .enumerate()
    {
        if expected.is_nan() {
            assert!(
                actual.is_nan(),
                "{label}[{index}]: expected NaN, got {actual}"
            );
        } else {
            assert_relative_eq!(*actual, *expected, epsilon = 1e-12);
        }
    }
}

fn assert_float_close_or_nan(actual: f64, expected: f64, label: &str) {
    if expected.is_nan() {
        assert!(actual.is_nan(), "{label}: expected NaN, got {actual}");
        return;
    }
    assert_relative_eq!(actual, expected, epsilon = 1e-12);
}

fn assert_option_close(actual: Option<f64>, expected: Option<f64>, label: &str) {
    match (actual, expected) {
        (Some(actual), Some(expected)) => assert_float_close_or_nan(actual, expected, label),
        (None, None) => {}
        _ => panic!("{label}: actual={actual:?}, expected={expected:?}"),
    }
}

fn assert_slice_close_or_nan(actual: &[f64], expected: &[f64], label: &str) {
    assert_eq!(actual.len(), expected.len(), "{label}: length mismatch");
    for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
        assert_float_close_or_nan(*actual, *expected, &format!("{label}[{index}]"));
    }
}

fn assert_lrt_fit_state_matches(actual: &DeseqFit, expected: &DeseqFit, label: &str) {
    assert_matrix_close_or_nan(
        actual.beta.as_ref().unwrap(),
        expected.beta.as_ref().unwrap(),
        &format!("{label} beta"),
    );
    assert_matrix_close_or_nan(
        actual.beta_se.as_ref().unwrap(),
        expected.beta_se.as_ref().unwrap(),
        &format!("{label} beta_se"),
    );
    assert_matrix_close_or_nan(
        actual.beta_covariance.as_ref().unwrap(),
        expected.beta_covariance.as_ref().unwrap(),
        &format!("{label} beta_covariance"),
    );
    assert_slice_close_or_nan(
        actual.log_like.as_ref().unwrap(),
        expected.log_like.as_ref().unwrap(),
        &format!("{label} full log_like"),
    );
    assert_slice_close_or_nan(
        actual.full_deviance.as_ref().unwrap(),
        expected.full_deviance.as_ref().unwrap(),
        &format!("{label} full_deviance"),
    );
    assert_slice_close_or_nan(
        actual.reduced_log_like.as_ref().unwrap(),
        expected.reduced_log_like.as_ref().unwrap(),
        &format!("{label} reduced log_like"),
    );
    assert_matrix_close_or_nan(
        actual.mu.as_ref().unwrap(),
        expected.mu.as_ref().unwrap(),
        &format!("{label} full mu"),
    );
    assert_matrix_close_or_nan(
        actual.hat_diagonal.as_ref().unwrap(),
        expected.hat_diagonal.as_ref().unwrap(),
        &format!("{label} full hat"),
    );
    assert_matrix_close_or_nan(
        actual.reduced_mu.as_ref().unwrap(),
        expected.reduced_mu.as_ref().unwrap(),
        &format!("{label} reduced mu"),
    );
    assert_matrix_close_or_nan(
        actual.reduced_hat_diagonal.as_ref().unwrap(),
        expected.reduced_hat_diagonal.as_ref().unwrap(),
        &format!("{label} reduced hat"),
    );
}

#[test]
fn native_glm_mu_lrt_preserves_diagnostics_and_all_zero_rows() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();

    let (fit, results) = native_lrt_builder()
        .fit_lrt_glm_mu(&counts, &full, &reduced, 1)
        .unwrap();

    assert_eq!(fit.design.as_ref().unwrap(), &full);
    assert_eq!(fit.reduced_design.as_ref().unwrap(), &reduced);
    assert_eq!(
        fit.all_zero,
        vec![true, false, false, false, false, false, false, false]
    );
    assert!(fit.disp_prior_var.unwrap().is_finite());
    assert_eq!(fit.disp_gene_est.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.disp_gene_iter.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.disp_fit.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.disp_map.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.dispersion.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.disp_gene_iter.as_ref().unwrap()[0], 0);
    assert!(fit.disp_gene_iter.as_ref().unwrap()[1..]
        .iter()
        .all(|iterations| *iterations > 0));

    assert_eq!(fit.beta.as_ref().unwrap().n_cols(), full.n_coefficients());
    assert_eq!(
        fit.reduced_log_like.as_ref().unwrap().len(),
        counts.n_genes()
    );
    assert_eq!(fit.reduced_mu.as_ref().unwrap().n_rows(), counts.n_genes());
    assert_eq!(
        fit.reduced_mu.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert_eq!(
        fit.reduced_hat_diagonal.as_ref().unwrap().n_rows(),
        counts.n_genes()
    );
    assert_eq!(
        fit.reduced_hat_diagonal.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert!(fit.reduced_mu.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert!(fit.reduced_hat_diagonal.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert!(fit.reduced_mu.as_ref().unwrap().row(1).unwrap()[0].is_finite());
    assert!(fit.reduced_hat_diagonal.as_ref().unwrap().row(1).unwrap()[0].is_finite());
    assert_eq!(
        fit.reduced_beta_converged.as_deref(),
        Some(&[false, true, true, true, true, true, true, true][..])
    );
    assert_eq!(fit.lrt.as_ref().unwrap().degrees_of_freedom, 1);
    assert_eq!(fit.lrt.as_ref().unwrap().deviance[0], None);
    assert_eq!(fit.lrt.as_ref().unwrap().pvalue[0], None);
    assert_eq!(fit.cooks.as_ref().unwrap().n_rows(), counts.n_genes());

    assert_eq!(results.rows.len(), counts.n_genes());
    assert_eq!(results.rows[0].gene.as_deref(), Some("zero"));
    assert_eq!(results.rows[0].pvalue, None);
    assert_eq!(results.rows[0].padj, None);
    assert_eq!(results.rows[0].converged, None);
    assert!(results.rows[1].stat.unwrap().is_finite());
    assert!(results.rows[1].pvalue.unwrap().is_finite());
    assert_eq!(results.rows[1].pvalue, fit.lrt.as_ref().unwrap().pvalue[1]);
}

#[test]
fn native_glm_mu_lrt_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();

    let (native_fit, native_results) = builder.fit_lrt_glm_mu(&counts, &full, &reduced, 1).unwrap();
    let final_dispersions = native_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_fit, fixed_results) = builder
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &final_dispersions, 1)
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
            native_fit.lrt.as_ref().unwrap().deviance[gene],
            fixed_fit.lrt.as_ref().unwrap().deviance[gene]
        );
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
        for sample in 0..counts.n_samples() {
            let native_mu = native_fit.reduced_mu.as_ref().unwrap().row(gene).unwrap()[sample];
            let fixed_mu = fixed_fit.reduced_mu.as_ref().unwrap().row(gene).unwrap()[sample];
            if native_mu.is_nan() {
                assert!(fixed_mu.is_nan());
            } else {
                assert_relative_eq!(native_mu, fixed_mu, epsilon = 1e-12);
            }
            let native_hat = native_fit
                .reduced_hat_diagonal
                .as_ref()
                .unwrap()
                .row(gene)
                .unwrap()[sample];
            let fixed_hat = fixed_fit
                .reduced_hat_diagonal
                .as_ref()
                .unwrap()
                .row(gene)
                .unwrap()[sample];
            if native_hat.is_nan() {
                assert!(fixed_hat.is_nan());
            } else {
                assert_relative_eq!(native_hat, fixed_hat, epsilon = 1e-12);
            }
        }
    }
}

#[test]
fn native_glm_mu_lrt_contrast_keeps_lrt_pvalues_and_reports_contrast_effect() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let contrast = [0.0, 1.0];

    let (coefficient_fit, coefficient_results) =
        builder.fit_lrt_glm_mu(&counts, &full, &reduced, 1).unwrap();
    let (contrast_fit, contrast_results) = builder
        .fit_lrt_glm_mu_contrast(&counts, &full, &reduced, &contrast)
        .unwrap();
    let final_dispersions = coefficient_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_contrast_fit, fixed_contrast_results) = builder
        .fit_fixed_dispersion_lrt_contrast(&counts, &full, &reduced, &final_dispersions, &contrast)
        .unwrap();

    assert_lrt_fit_state_matches(&contrast_fit, &coefficient_fit, "LRT contrast");
    assert_lrt_fit_state_matches(
        &fixed_contrast_fit,
        &coefficient_fit,
        "fixed-dispersion LRT contrast",
    );
    assert_eq!(
        contrast_results.metadata.test_type,
        Some(TestType::Lrt),
        "contrast result remains an LRT table"
    );
    assert_eq!(
        contrast_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        contrast_results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );
    for gene in 0..counts.n_genes() {
        assert_eq!(
            contrast_results.rows[gene].stat,
            coefficient_results.rows[gene].stat
        );
        assert_eq!(
            contrast_results.rows[gene].pvalue,
            coefficient_results.rows[gene].pvalue
        );
        assert_eq!(
            contrast_results.rows[gene].padj,
            coefficient_results.rows[gene].padj
        );
        assert_option_close(
            contrast_results.rows[gene].log2_fold_change,
            coefficient_results.rows[gene].log2_fold_change,
            &format!("LRT contrast LFC gene {gene}"),
        );
        assert_option_close(
            contrast_results.rows[gene].lfc_se,
            coefficient_results.rows[gene].lfc_se,
            &format!("LRT contrast SE gene {gene}"),
        );
        assert_eq!(
            fixed_contrast_results.rows[gene], contrast_results.rows[gene],
            "fixed/native LRT contrast row {gene}"
        );
    }
}

#[test]
fn native_linear_mu_lrt_contrast_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let contrast = [0.0, 1.0];

    let (coefficient_fit, coefficient_results) = builder
        .fit_lrt_linear_mu(&counts, &full, &reduced, 1)
        .unwrap();
    let (contrast_fit, contrast_results) = builder
        .fit_lrt_linear_mu_contrast(&counts, &full, &reduced, &contrast)
        .unwrap();
    let final_dispersions = coefficient_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_contrast_fit, fixed_contrast_results) = builder
        .fit_fixed_dispersion_lrt_contrast(&counts, &full, &reduced, &final_dispersions, &contrast)
        .unwrap();

    assert_lrt_fit_state_matches(&contrast_fit, &coefficient_fit, "linear-mu LRT contrast");
    assert_lrt_fit_state_matches(
        &fixed_contrast_fit,
        &coefficient_fit,
        "linear-mu fixed LRT contrast",
    );
    assert_eq!(contrast_results.rows, fixed_contrast_results.rows);
    assert_eq!(
        contrast_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        contrast_results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );
    for gene in 0..counts.n_genes() {
        assert_eq!(
            contrast_results.rows[gene].stat,
            coefficient_results.rows[gene].stat
        );
        assert_eq!(
            contrast_results.rows[gene].pvalue,
            coefficient_results.rows[gene].pvalue
        );
        assert_eq!(
            contrast_results.rows[gene].padj,
            coefficient_results.rows[gene].padj
        );
    }
}

#[test]
fn native_linear_mu_lrt_contrast_specs_set_metadata_and_factor_levels() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let spec = ContrastSpec::coefficient_name("condition_B_vs_A");
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let factor_contrast = FactorLevelContrast {
        factor: "condition",
        numerator: "B",
        denominator: "A",
        reference: None,
        sample_levels: &levels,
    };

    let (numeric_fit, numeric_results) = builder
        .fit_lrt_linear_mu_contrast(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    let (named_fit, named_results) = builder
        .fit_lrt_linear_mu_contrast_spec(&counts, &full, &reduced, &spec)
        .unwrap();
    let (factor_fit, factor_results) = builder
        .fit_lrt_linear_mu_factor_level_contrast(&counts, &full, &reduced, factor_contrast)
        .unwrap();

    assert_lrt_fit_state_matches(&named_fit, &numeric_fit, "linear-mu named LRT contrast");
    assert_lrt_fit_state_matches(
        &factor_fit,
        &numeric_fit,
        "linear-mu factor-level LRT contrast",
    );
    assert_eq!(named_results.rows, numeric_results.rows);
    assert_eq!(factor_results.rows, numeric_results.rows);
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
}

#[test]
fn native_lrt_parametric_contrast_helpers_ignore_builder_fit_type() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let mean_builder = native_lrt_builder().fit_type(FitType::Mean);
    let parametric_builder = mean_builder.clone().fit_type(FitType::Parametric);
    let spec = ContrastSpec::coefficient_name("condition_B_vs_A");
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let factor_contrast = FactorLevelContrast {
        factor: "condition",
        numerator: "B",
        denominator: "A",
        reference: None,
        sample_levels: &levels,
    };

    let (linear_parametric_fit, linear_parametric_results) = mean_builder
        .fit_lrt_linear_mu_contrast_parametric(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    let (linear_expected_fit, linear_expected_results) = parametric_builder
        .fit_lrt_linear_mu_contrast(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    assert_lrt_fit_state_matches(
        &linear_parametric_fit,
        &linear_expected_fit,
        "linear-mu parametric LRT contrast",
    );
    assert_eq!(linear_parametric_results, linear_expected_results);

    let (linear_named_fit, linear_named_results) = mean_builder
        .fit_lrt_linear_mu_contrast_spec_parametric(&counts, &full, &reduced, &spec)
        .unwrap();
    assert_lrt_fit_state_matches(
        &linear_named_fit,
        &linear_expected_fit,
        "linear-mu named parametric LRT contrast",
    );
    assert_eq!(linear_named_results.rows, linear_expected_results.rows);
    assert_eq!(
        linear_named_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );

    let (linear_factor_fit, linear_factor_results) = mean_builder
        .fit_lrt_linear_mu_factor_level_contrast_parametric(
            &counts,
            &full,
            &reduced,
            factor_contrast,
        )
        .unwrap();
    assert_lrt_fit_state_matches(
        &linear_factor_fit,
        &linear_expected_fit,
        "linear-mu factor parametric LRT contrast",
    );
    assert_eq!(linear_factor_results.rows, linear_expected_results.rows);
    assert_eq!(
        linear_factor_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );

    let (glm_parametric_fit, glm_parametric_results) = mean_builder
        .fit_lrt_glm_mu_contrast_parametric(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    let (glm_expected_fit, glm_expected_results) = parametric_builder
        .fit_lrt_glm_mu_contrast(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    assert_lrt_fit_state_matches(
        &glm_parametric_fit,
        &glm_expected_fit,
        "GLM-mu parametric LRT contrast",
    );
    assert_eq!(glm_parametric_results, glm_expected_results);

    let (glm_named_fit, glm_named_results) = mean_builder
        .fit_lrt_glm_mu_contrast_spec_parametric(&counts, &full, &reduced, &spec)
        .unwrap();
    assert_lrt_fit_state_matches(
        &glm_named_fit,
        &glm_expected_fit,
        "GLM-mu named parametric LRT contrast",
    );
    assert_eq!(glm_named_results.rows, glm_expected_results.rows);
    assert_eq!(
        glm_named_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );

    let (glm_factor_fit, glm_factor_results) = mean_builder
        .fit_lrt_glm_mu_factor_level_contrast_parametric(&counts, &full, &reduced, factor_contrast)
        .unwrap();
    assert_lrt_fit_state_matches(
        &glm_factor_fit,
        &glm_expected_fit,
        "GLM-mu factor parametric LRT contrast",
    );
    assert_eq!(glm_factor_results.rows, glm_expected_results.rows);
    assert_eq!(
        glm_factor_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
}

#[test]
fn native_glm_mu_lrt_contrast_all_zero_only_zeroes_lfc() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        5,
        6,
        vec![
            0, 0, 0, 0, 0, 0, //
            20, 22, 0, 0, 0, 0, //
            12, 28, 14, 26, 11, 29, //
            55, 105, 50, 110, 65, 95, //
            120, 200, 115, 205, 125, 195,
        ],
        Some(
            [
                "zero",
                "contrast_groups_zero",
                "variable",
                "high_up",
                "stable",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
        ),
        None,
    )
    .unwrap();
    let full = DesignMatrix::from_row_major(
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
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "condition_D_vs_A".to_string(),
        ]),
    )
    .unwrap();
    let reduced = DesignMatrix::from_row_major(
        6,
        1,
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".to_string()]),
    )
    .unwrap();
    let builder = native_lrt_builder().size_factors(vec![1.0; 6]);
    let contrast = [0.0, -1.0, 1.0];

    let (_coefficient_fit, coefficient_results) =
        builder.fit_lrt_glm_mu(&counts, &full, &reduced, 2).unwrap();
    let (_contrast_fit, contrast_results) = builder
        .fit_lrt_glm_mu_contrast(&counts, &full, &reduced, &contrast)
        .unwrap();

    assert_eq!(
        contrast_results.rows[1].gene.as_deref(),
        Some("contrast_groups_zero")
    );
    assert_eq!(contrast_results.rows[1].log2_fold_change, Some(0.0));
    assert_eq!(
        contrast_results.rows[1].stat,
        coefficient_results.rows[1].stat
    );
    assert_eq!(
        contrast_results.rows[1].pvalue,
        coefficient_results.rows[1].pvalue
    );
    assert_eq!(
        contrast_results.rows[1].padj,
        coefficient_results.rows[1].padj
    );
    assert!(contrast_results.rows[1].stat.unwrap().is_finite());
    assert!(contrast_results.rows[1].pvalue.unwrap().is_finite());
}

#[test]
fn native_glm_mu_lrt_factor_level_contrast_uses_character_all_zero_cleanup() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        5,
        6,
        vec![
            0, 0, 0, 0, 0, 0, //
            20, 22, 0, 0, 0, 0, //
            12, 28, 14, 26, 11, 29, //
            55, 105, 50, 110, 65, 95, //
            120, 200, 115, 205, 125, 195,
        ],
        Some(
            [
                "zero",
                "contrast_groups_zero",
                "variable",
                "high_up",
                "stable",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
        ),
        None,
    )
    .unwrap();
    let full = DesignMatrix::from_row_major(
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
            "Intercept".to_string(),
            "condition_B_vs_A".to_string(),
            "condition_D_vs_A".to_string(),
        ]),
    )
    .unwrap();
    let reduced = DesignMatrix::from_row_major(
        6,
        1,
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".to_string()]),
    )
    .unwrap();
    let levels = ["A", "A", "B", "B", "D", "D"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let builder = native_lrt_builder().size_factors(vec![1.0; 6]);
    let contrast = FactorLevelContrast {
        factor: "condition",
        numerator: "D",
        denominator: "B",
        reference: Some("A"),
        sample_levels: &levels,
    };

    let (_coefficient_fit, coefficient_results) =
        builder.fit_lrt_glm_mu(&counts, &full, &reduced, 2).unwrap();
    let (_contrast_fit, contrast_results) = builder
        .fit_lrt_glm_mu_factor_level_contrast(&counts, &full, &reduced, contrast)
        .unwrap();

    assert_eq!(
        contrast_results.rows[1].gene.as_deref(),
        Some("contrast_groups_zero")
    );
    assert_eq!(contrast_results.rows[1].log2_fold_change, Some(0.0));
    assert_eq!(
        contrast_results.rows[1].stat,
        coefficient_results.rows[1].stat
    );
    assert_eq!(
        contrast_results.rows[1].pvalue,
        coefficient_results.rows[1].pvalue
    );
    assert_eq!(
        contrast_results.metadata.result_name.as_deref(),
        Some("condition_D_vs_B")
    );
    assert_eq!(
        contrast_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition D vs B")
    );
}

#[test]
fn top_level_lrt_contrast_spec_resolves_named_full_model_effect() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let spec = ContrastSpec::coefficient_name("condition_B_vs_A");

    let (numeric_fit, numeric_results) = builder
        .fit_lrt_with_results_contrast(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    let (named_fit, named_results) = builder
        .fit_lrt_with_results_contrast_spec(&counts, &full, &reduced, &spec)
        .unwrap();

    assert_lrt_fit_state_matches(&named_fit, &numeric_fit, "named LRT contrast");
    assert_eq!(named_results.rows, numeric_results.rows);
    assert_eq!(
        named_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        named_results.metadata.comparison.as_deref(),
        Some("coefficient condition_B_vs_A")
    );
}

#[test]
fn top_level_lrt_fit_only_helpers_match_result_helpers() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let spec = ContrastSpec::coefficient_name("condition_B_vs_A");

    let default_fit = builder.fit_lrt(&counts, &full, &reduced).unwrap();
    let (default_result_fit, _default_results) = builder
        .fit_lrt_with_results(&counts, &full, &reduced)
        .unwrap();
    assert_lrt_fit_state_matches(&default_fit, &default_result_fit, "fit-only default LRT");

    let named_fit = builder
        .fit_lrt_name(&counts, &full, &reduced, "condition_B_vs_A")
        .unwrap();
    let (named_result_fit, _named_results) = builder
        .fit_lrt_with_results_name(&counts, &full, &reduced, "condition_B_vs_A")
        .unwrap();
    assert_lrt_fit_state_matches(&named_fit, &named_result_fit, "fit-only named LRT");

    let contrast_fit = builder
        .fit_lrt_contrast(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    let (contrast_result_fit, _contrast_results) = builder
        .fit_lrt_with_results_contrast(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    assert_lrt_fit_state_matches(
        &contrast_fit,
        &contrast_result_fit,
        "fit-only numeric LRT contrast",
    );

    let spec_fit = builder
        .fit_lrt_contrast_spec(&counts, &full, &reduced, &spec)
        .unwrap();
    let (spec_result_fit, _spec_results) = builder
        .fit_lrt_with_results_contrast_spec(&counts, &full, &reduced, &spec)
        .unwrap();
    assert_lrt_fit_state_matches(&spec_fit, &spec_result_fit, "fit-only named LRT contrast");

    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let factor_contrast = FactorLevelContrast {
        factor: "condition",
        numerator: "B",
        denominator: "A",
        reference: None,
        sample_levels: &levels,
    };
    let factor_fit = builder
        .fit_lrt_factor_level_contrast(&counts, &full, &reduced, factor_contrast)
        .unwrap();
    let (factor_result_fit, _factor_results) = builder
        .fit_lrt_with_results_factor_level_contrast(&counts, &full, &reduced, factor_contrast)
        .unwrap();
    assert_lrt_fit_state_matches(
        &factor_fit,
        &factor_result_fit,
        "fit-only factor-level LRT contrast",
    );
}

#[test]
fn native_glm_mu_lrt_cooks_replacement_refit_merges_refit_rows() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();

    let output = builder
        .fit_lrt_glm_mu_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            1,
            &CooksReplacementOptions {
                trim: 0.2,
                cooks_cutoff: 0.0,
                min_replicates: 3,
                which_samples: Some(vec![true, false, false, false, false, false, false, false]),
            },
        )
        .unwrap();

    assert!(output.original_fit.lrt.is_some());
    assert!(output.refit_plan.n_refit > 0);
    assert!(output.refit_plan.should_refit);
    assert!(!output.refit_plan.refit_rows.is_empty());
    assert!(output.refit_fit.as_ref().unwrap().lrt.is_some());
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
            output.results.rows[gene].stat,
            refit_results.rows[gene].stat
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
        assert_eq!(output.results.rows[gene].stat, Some(0.0));
        assert_eq!(output.results.rows[gene].log2_fold_change, Some(0.0));
        assert_eq!(output.results.rows[gene].lfc_se, Some(0.0));
        assert_eq!(output.results.rows[gene].dispersion, None);
    }
}

#[test]
fn native_glm_mu_lrt_contrast_cooks_replacement_refit_merges_refit_rows() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![true, false, false, false, false, false, false, false]),
    };

    let output = builder
        .fit_lrt_glm_mu_contrast_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            &[0.0, 1.0],
            &options,
        )
        .unwrap();
    let named_output = builder
        .fit_lrt_glm_mu_contrast_spec_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
            &options,
        )
        .unwrap();

    assert!(output.original_fit.lrt.is_some());
    assert!(output.refit_plan.should_refit);
    assert!(output.refit_fit.as_ref().unwrap().lrt.is_some());
    let refit_results = output.refit_results.as_ref().unwrap();
    for gene in output.refit_plan.refit_rows.iter().copied() {
        assert_eq!(
            output.results.rows[gene].log2_fold_change,
            refit_results.rows[gene].log2_fold_change
        );
        assert_eq!(
            output.results.rows[gene].stat,
            refit_results.rows[gene].stat
        );
        assert_eq!(
            output.results.rows[gene].pvalue,
            refit_results.rows[gene].pvalue
        );
    }
    assert_eq!(named_output.results.rows, output.results.rows);
    assert_eq!(
        named_output.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        named_output
            .refit_results
            .as_ref()
            .unwrap()
            .metadata
            .comparison
            .as_deref(),
        Some("coefficient condition_B_vs_A")
    );
}

#[test]
fn top_level_fit_lrt_runs_default_glm_mu_lrt_for_last_coefficient() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();

    let top_level_fit = builder.fit_lrt(&counts, &full, &reduced).unwrap();
    let (lrt_fit, _results) = builder.fit_lrt_glm_mu(&counts, &full, &reduced, 1).unwrap();

    assert_eq!(top_level_fit.counts_summary, lrt_fit.counts_summary);
    assert_eq!(top_level_fit.design, lrt_fit.design);
    assert_eq!(top_level_fit.reduced_design, lrt_fit.reduced_design);
    assert_eq!(top_level_fit.all_zero, lrt_fit.all_zero);
    assert_eq!(
        top_level_fit.reduced_beta_converged,
        lrt_fit.reduced_beta_converged
    );
    assert_lrt_fit_state_matches(&top_level_fit, &lrt_fit, "top-level LRT");
    assert_eq!(top_level_fit.lrt, lrt_fit.lrt);
    assert!(top_level_fit.dispersion.is_some());
    assert!(matches!(
        top_level_fit.dispersion_trend.as_ref(),
        Some(DispersionTrendFit::Mean(_))
    ));
    assert!(top_level_fit.beta.is_some());
    assert!(top_level_fit.lrt.is_some());

    let transformed = top_level_fit.vst(&counts).unwrap();
    assert_eq!(transformed.n_rows(), counts.n_genes());
    assert_eq!(transformed.n_cols(), counts.n_samples());
    assert_eq!(
        transformed,
        top_level_fit
            .variance_stabilizing_transform(&counts)
            .unwrap()
    );
    assert!(transformed.as_slice().iter().all(|value| value.is_finite()));
}

#[test]
fn top_level_fit_lrt_with_results_returns_default_lrt_results() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();

    let (top_level_fit, top_level_results) = builder
        .fit_lrt_with_results(&counts, &full, &reduced)
        .unwrap();
    let (lrt_fit, lrt_results) = builder.fit_lrt_glm_mu(&counts, &full, &reduced, 1).unwrap();

    assert_eq!(top_level_fit.lrt, lrt_fit.lrt);
    assert_lrt_fit_state_matches(&top_level_fit, &lrt_fit, "top-level LRT results");
    assert_eq!(top_level_results, lrt_results);
    assert_eq!(top_level_results.metadata.test_type, Some(TestType::Lrt));
    assert_eq!(
        top_level_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
}

#[test]
fn top_level_fit_lrt_with_results_accepts_coefficient_name() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();

    let (named_fit, named_results) = builder
        .fit_lrt_with_results_name(&counts, &full, &reduced, "condition_B_vs_A")
        .unwrap();
    let (indexed_fit, indexed_results) =
        builder.fit_lrt_glm_mu(&counts, &full, &reduced, 1).unwrap();

    assert_eq!(named_fit.lrt, indexed_fit.lrt);
    assert_lrt_fit_state_matches(&named_fit, &indexed_fit, "named LRT");
    assert_eq!(named_results, indexed_results);
    assert!(builder
        .fit_lrt_with_results_name(&counts, &full, &reduced, "missing")
        .is_err());
}

#[test]
fn top_level_fit_lrt_with_results_cooks_replacement_runs_default_glm_mu_lrt() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![true, false, false, false, false, false, false, false]),
    };

    let top_level_output = builder
        .fit_lrt_with_results_with_cooks_replacement(&counts, &full, &reduced, &replacement_options)
        .unwrap();
    let lrt_output = builder
        .fit_lrt_glm_mu_with_cooks_replacement(&counts, &full, &reduced, 1, &replacement_options)
        .unwrap();

    assert_eq!(top_level_output.refit_plan, lrt_output.refit_plan);
    assert_eq!(
        top_level_output.original_results,
        lrt_output.original_results
    );
    assert_lrt_fit_state_matches(
        &top_level_output.original_fit,
        &lrt_output.original_fit,
        "top-level LRT replacement original",
    );
    assert_eq!(top_level_output.refit_results, lrt_output.refit_results);
    assert_lrt_fit_state_matches(
        top_level_output.refit_fit.as_ref().unwrap(),
        lrt_output.refit_fit.as_ref().unwrap(),
        "top-level LRT replacement refit",
    );
    assert_eq!(top_level_output.results, lrt_output.results);
    assert!(top_level_output.refit_plan.should_refit);
    assert!(top_level_output.refit_fit.as_ref().unwrap().lrt.is_some());
    assert_eq!(
        top_level_output.results.metadata.test_type,
        Some(TestType::Lrt)
    );
    assert_eq!(
        top_level_output.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
}

#[test]
fn top_level_fit_lrt_with_results_cooks_replacement_accepts_coefficient_name() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![true, false, false, false, false, false, false, false]),
    };

    let named_output = builder
        .fit_lrt_with_results_name_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            "condition_B_vs_A",
            &replacement_options,
        )
        .unwrap();
    let indexed_output = builder
        .fit_lrt_glm_mu_with_cooks_replacement(&counts, &full, &reduced, 1, &replacement_options)
        .unwrap();

    assert_eq!(named_output.refit_plan, indexed_output.refit_plan);
    assert_eq!(named_output.results, indexed_output.results);
    assert_lrt_fit_state_matches(
        &named_output.original_fit,
        &indexed_output.original_fit,
        "named LRT replacement original",
    );
}

#[test]
fn top_level_fit_lrt_with_results_cooks_replacement_validates_reduced_design() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let invalid_reduced = full_design();
    let err = native_lrt_builder()
        .fit_lrt_with_results_with_cooks_replacement(
            &counts,
            &full,
            &invalid_reduced,
            &CooksReplacementOptions::new(f64::MAX),
        )
        .unwrap_err();

    assert!(matches!(err, DeseqError::InvalidDimensions { .. }));
}

#[test]
fn native_glm_mu_lrt_cooks_replacement_refit_skips_when_no_rows_are_marked() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();

    let (expected_fit, expected_results) = builder
        .clone()
        .fit_lrt_glm_mu(&counts, &full, &reduced, 1)
        .unwrap();
    let output = builder
        .fit_lrt_glm_mu_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
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
    assert_lrt_fit_state_matches(
        &output.original_fit,
        &expected_fit,
        "native GLM-mu LRT no-refit original",
    );
    assert_eq!(output.original_results, expected_results);
    assert_eq!(output.results, expected_results);
    assert_eq!(output.results, output.original_results);
    assert_eq!(
        output.refit_plan.replacement.replaced_counts.as_slice(),
        counts.as_slice()
    );
}

#[test]
fn native_linear_mu_lrt_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();

    let (native_fit, native_results) = builder
        .fit_lrt_linear_mu(&counts, &full, &reduced, 1)
        .unwrap();
    let final_dispersions = native_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_fit, fixed_results) = builder
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &final_dispersions, 1)
        .unwrap();

    assert!(native_fit.disp_prior_var.unwrap().is_finite());
    assert_eq!(native_fit.lrt.as_ref().unwrap().degrees_of_freedom, 1);
    assert_eq!(native_results.rows[0].pvalue, None);

    for gene in 0..counts.n_genes() {
        assert_eq!(
            native_fit.lrt.as_ref().unwrap().pvalue[gene],
            fixed_fit.lrt.as_ref().unwrap().pvalue[gene]
        );
        assert_eq!(
            native_results.rows[gene].pvalue,
            fixed_results.rows[gene].pvalue
        );
        assert_eq!(
            native_results.rows[gene].stat,
            fixed_results.rows[gene].stat
        );
    }
}

#[test]
fn native_linear_mu_local_lrt_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder().fit_type(FitType::Local);

    let (native_fit, native_results) = builder
        .fit_lrt_linear_mu(&counts, &full, &reduced, 1)
        .unwrap();
    let final_dispersions = native_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_fit, fixed_results) = builder
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &final_dispersions, 1)
        .unwrap();

    assert!(native_fit.disp_prior_var.unwrap().is_finite());
    assert_eq!(
        native_fit.disp_fit.as_ref().unwrap().len(),
        counts.n_genes()
    );
    assert_eq!(native_fit.lrt.as_ref().unwrap().degrees_of_freedom, 1);
    assert_eq!(native_results.rows[0].pvalue, None);

    for gene in 0..counts.n_genes() {
        assert_eq!(
            native_fit.lrt.as_ref().unwrap().pvalue[gene],
            fixed_fit.lrt.as_ref().unwrap().pvalue[gene]
        );
        assert_eq!(
            native_results.rows[gene].pvalue,
            fixed_results.rows[gene].pvalue
        );
        assert_eq!(
            native_results.rows[gene].stat,
            fixed_results.rows[gene].stat
        );
    }
}

#[test]
fn native_glm_mu_local_lrt_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder().fit_type(FitType::Local);

    let (native_fit, native_results) = builder.fit_lrt_glm_mu(&counts, &full, &reduced, 1).unwrap();
    let final_dispersions = native_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_fit, fixed_results) = builder
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &final_dispersions, 1)
        .unwrap();

    assert!(native_fit.disp_prior_var.unwrap().is_finite());
    assert_eq!(
        native_fit.disp_fit.as_ref().unwrap().len(),
        counts.n_genes()
    );
    assert_eq!(native_fit.lrt.as_ref().unwrap().degrees_of_freedom, 1);
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
            native_fit.lrt.as_ref().unwrap().pvalue[gene],
            fixed_fit.lrt.as_ref().unwrap().pvalue[gene]
        );
        assert_eq!(
            native_results.rows[gene].pvalue,
            fixed_results.rows[gene].pvalue
        );
        assert_eq!(
            native_results.rows[gene].stat,
            fixed_results.rows[gene].stat
        );
    }
}
