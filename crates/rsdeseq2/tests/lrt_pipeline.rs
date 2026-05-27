use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

fn assert_lrt_likelihood_state(fit: &DeseqFit, counts: &CountMatrix) {
    let log_like = fit.log_like.as_ref().unwrap();
    let full_deviance = fit.full_deviance.as_ref().unwrap();
    let reduced_log_like = fit.reduced_log_like.as_ref().unwrap();
    let lrt = fit.lrt.as_ref().unwrap();
    assert_eq!(log_like.len(), counts.n_genes());
    assert_eq!(full_deviance.len(), counts.n_genes());
    assert_eq!(reduced_log_like.len(), counts.n_genes());
    assert_eq!(lrt.deviance.len(), counts.n_genes());
    for gene in 0..counts.n_genes() {
        if log_like[gene].is_nan() {
            assert!(
                full_deviance[gene].is_nan(),
                "full deviance for gene {gene} should be NaN when log_like is NaN"
            );
            assert_eq!(lrt.deviance[gene], None);
            continue;
        }
        assert_relative_eq!(full_deviance[gene], -2.0 * log_like[gene], epsilon = 1e-12);
        if reduced_log_like[gene].is_nan() {
            assert_eq!(lrt.deviance[gene], None);
        } else {
            assert_relative_eq!(
                lrt.deviance[gene].unwrap(),
                2.0 * (log_like[gene] - reduced_log_like[gene]),
                epsilon = 1e-12
            );
        }
    }
}

#[test]
fn fixed_dispersion_lrt_pipeline_fits_full_and_reduced_models() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        4,
        vec![10, 10, 20, 20],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let full = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let reduced = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &[0.05], 1)
        .unwrap();

    assert_eq!(fit.design.as_ref().unwrap().n_coefficients(), 2);
    assert_eq!(fit.reduced_design.as_ref().unwrap().n_coefficients(), 1);
    assert_eq!(fit.lrt.as_ref().unwrap().degrees_of_freedom, 1);
    assert_eq!(fit.lrt.as_ref().unwrap().reduced_converged, vec![true]);
    assert_eq!(fit.reduced_beta_converged.as_deref(), Some(&[true][..]));
    assert_eq!(fit.reduced_beta_iter.as_ref().unwrap().len(), 1);
    assert!(fit.reduced_beta_iter.as_ref().unwrap()[0] > 0);
    assert_eq!(fit.reduced_mu.as_ref().unwrap().n_rows(), 1);
    assert_eq!(
        fit.reduced_mu.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert_eq!(fit.reduced_hat_diagonal.as_ref().unwrap().n_rows(), 1);
    assert_eq!(
        fit.reduced_hat_diagonal.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert!(fit
        .reduced_mu
        .as_ref()
        .unwrap()
        .as_slice()
        .iter()
        .all(|value| value.is_finite()));
    assert!(fit
        .reduced_hat_diagonal
        .as_ref()
        .unwrap()
        .as_slice()
        .iter()
        .all(|value| value.is_finite()));
    assert_lrt_likelihood_state(&fit, &counts);
    assert_relative_eq!(
        results.rows[0].log2_fold_change.unwrap(),
        2.0_f64.log2(),
        epsilon = 1e-8
    );
    assert!(results.rows[0].stat.unwrap() > 0.0);
    assert!(results.rows[0].pvalue.unwrap() < 1.0);
    assert_eq!(results.rows[0].pvalue, fit.lrt.as_ref().unwrap().pvalue[0]);
    assert_eq!(fit.cooks.as_ref().unwrap().n_cols(), 4);
}

#[test]
fn fixed_dispersion_lrt_pipeline_uses_normalization_factors_for_offsets() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        4,
        vec![10, 20, 20, 40],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let normalization_factors =
        RowMajorMatrix::from_row_major(1, 4, vec![1.0, 2.0, 1.0, 2.0]).unwrap();
    let full = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let reduced = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![100.0, 100.0, 100.0, 100.0])
        .normalization_factors(normalization_factors.clone())
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &[0.05], 1)
        .unwrap();

    assert_eq!(fit.normalization_factors, Some(normalization_factors));
    assert_relative_eq!(fit.base_mean[0], 15.0, epsilon = 1e-12);
    assert_relative_eq!(
        fit.beta.as_ref().unwrap().as_slice()[1],
        2.0_f64.log2(),
        epsilon = 1e-8
    );
    assert_relative_eq!(
        results.rows[0].log2_fold_change.unwrap(),
        2.0_f64.log2(),
        epsilon = 1e-8
    );
    assert_lrt_likelihood_state(&fit, &counts);
    assert_eq!(results.rows[0].pvalue, fit.lrt.as_ref().unwrap().pvalue[0]);
}

#[test]
fn fixed_dispersion_lrt_pipeline_expands_all_zero_rows() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        2,
        4,
        vec![0, 0, 0, 0, 10, 10, 20, 20],
        Some(vec!["zero_gene".into(), "signal_gene".into()]),
        None,
    )
    .unwrap();
    let full =
        DesignMatrix::from_row_major(4, 2, vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0], None)
            .unwrap();
    let reduced = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &[0.1, 0.05], 1)
        .unwrap();

    assert_eq!(fit.all_zero, vec![true, false]);
    assert!(fit.beta.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert_eq!(
        fit.reduced_beta_converged.as_deref(),
        Some(&[false, true][..])
    );
    assert_eq!(fit.reduced_beta_iter.as_ref().unwrap()[0], 0);
    assert!(fit.reduced_beta_iter.as_ref().unwrap()[1] > 0);
    assert!(fit.reduced_mu.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert!(fit.reduced_hat_diagonal.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert!(fit
        .reduced_mu
        .as_ref()
        .unwrap()
        .row(1)
        .unwrap()
        .iter()
        .all(|value| value.is_finite()));
    assert!(fit
        .reduced_hat_diagonal
        .as_ref()
        .unwrap()
        .row(1)
        .unwrap()
        .iter()
        .all(|value| value.is_finite()));
    assert!(fit.full_deviance.as_ref().unwrap()[0].is_nan());
    assert_lrt_likelihood_state(&fit, &counts);
    assert!(fit.reduced_log_like.as_ref().unwrap()[0].is_nan());
    assert!(fit.reduced_log_like.as_ref().unwrap()[1].is_finite());
    assert_eq!(fit.lrt.as_ref().unwrap().deviance[0], None);
    assert_eq!(fit.lrt.as_ref().unwrap().pvalue[0], None);
    assert_eq!(results.rows[0].gene.as_deref(), Some("zero_gene"));
    assert_eq!(results.rows[0].pvalue, None);
    assert_eq!(results.rows[0].padj, None);
    assert_eq!(results.rows[0].converged, None);
    assert!(results.rows[1].stat.unwrap() > 0.0);
    assert!(results.rows[1].pvalue.is_some());
}

#[test]
fn fixed_dispersion_lrt_factor_level_contrast_only_zeroes_lfc() {
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
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 6])
        .disable_cooks_cutoff()
        .disable_independent_filtering();
    let dispersions = vec![0.2; counts.n_genes()];
    let contrast = FactorLevelContrast {
        factor: "condition",
        numerator: "D",
        denominator: "B",
        reference: Some("A"),
        sample_levels: &levels,
    };

    let (_coefficient_fit, coefficient_results) = builder
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &dispersions, 2)
        .unwrap();
    let (_contrast_fit, contrast_results) = builder
        .fit_fixed_dispersion_lrt_factor_level_contrast(
            &counts,
            &full,
            &reduced,
            &dispersions,
            contrast,
        )
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
fn fixed_dispersion_lrt_factor_level_contrast_applies_low_count_cooks_gate() {
    let counts = CountMatrix::from_row_major_u32(1, 6, vec![1, 20, 21, 20, 20, 20]).unwrap();
    let full = DesignMatrix::from_row_major(
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
    let reduced = DesignMatrix::from_row_major(
        6,
        1,
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into()]),
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
        .fit_fixed_dispersion_lrt_factor_level_contrast(
            &counts,
            &full,
            &reduced,
            &[0.1],
            FactorLevelContrast::new("condition", "B", "A", &levels),
        )
        .unwrap();

    assert!(results.rows[0].max_cooks.unwrap() > 0.0);
    assert_eq!(results.rows[0].cooks_outlier, Some(false));
    assert!(results.rows[0].pvalue.is_some());
}

fn replacement_lrt_fixture() -> (
    CountMatrix,
    DesignMatrix,
    DesignMatrix,
    Vec<f64>,
    CooksReplacementOptions,
) {
    let counts = CountMatrix::from_row_major_u32(
        2,
        8,
        vec![
            0, 20, 1, 19, 2, 18, 3, 17, //
            12, 28, 10, 30, 14, 26, 11, 29,
        ],
    )
    .unwrap();
    let full = DesignMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let reduced = DesignMatrix::from_row_major(
        8,
        1,
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into()]),
    )
    .unwrap();
    let options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![true, false, false, false, false, false, false, false]),
    };
    (counts, full, reduced, vec![0.1, 0.1], options)
}

#[test]
fn fixed_dispersion_lrt_cooks_replacement_refits_marked_rows() {
    let (counts, full, reduced, dispersions, options) = replacement_lrt_fixture();

    let output = DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .disable_independent_filtering()
        .fit_fixed_dispersion_lrt_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            &dispersions,
            1,
            &options,
        )
        .unwrap();

    assert!(output.refit_plan.n_refit > 0);
    assert!(output.refit_plan.should_refit);
    assert!(output.refit_fit.is_some());
    assert!(output.refit_results.is_some());

    let refit_results = output.refit_results.as_ref().unwrap();
    for gene in output.refit_plan.refit_rows.iter().copied() {
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
    }
}

#[test]
fn fixed_dispersion_lrt_contrast_replacement_preserves_metadata() {
    let (counts, full, reduced, dispersions, options) = replacement_lrt_fixture();

    let named = DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .disable_independent_filtering()
        .fit_fixed_dispersion_lrt_contrast_spec_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            &dispersions,
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
            &options,
        )
        .unwrap();

    assert!(named.refit_plan.n_refit > 0);
    assert_eq!(
        named.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        named.results.metadata.comparison.as_deref(),
        Some("coefficient condition_B_vs_A")
    );
}

#[test]
fn fixed_dispersion_lrt_factor_level_replacement_preserves_metadata() {
    let (counts, full, reduced, dispersions, options) = replacement_lrt_fixture();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

    let output = DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .disable_independent_filtering()
        .fit_fixed_dispersion_lrt_factor_level_contrast_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            &dispersions,
            FactorLevelContrast::new("condition", "B", "A", &levels),
            &options,
        )
        .unwrap();

    assert!(output.refit_plan.n_refit > 0);
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
fn fixed_dispersion_lrt_replacement_skips_when_no_rows_are_marked() {
    let (counts, full, reduced, dispersions, _options) = replacement_lrt_fixture();

    let output = DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .disable_independent_filtering()
        .fit_fixed_dispersion_lrt_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            &dispersions,
            1,
            &CooksReplacementOptions::new(f64::MAX),
        )
        .unwrap();

    assert_eq!(output.refit_plan.n_refit, 0);
    assert!(!output.refit_plan.should_refit);
    assert!(output.refit_fit.is_none());
    assert!(output.refit_results.is_none());
    assert_eq!(
        output.results.rows.len(),
        output.original_results.rows.len()
    );
    for (final_row, original_row) in output
        .results
        .rows
        .iter()
        .zip(&output.original_results.rows)
    {
        assert_eq!(final_row.log2_fold_change, original_row.log2_fold_change);
        assert_eq!(final_row.stat, original_row.stat);
        assert_eq!(final_row.pvalue, original_row.pvalue);
        assert_eq!(final_row.max_cooks, original_row.max_cooks);
        assert_eq!(final_row.cooks_outlier, Some(false));
    }
    assert_eq!(
        output.refit_plan.replacement.replaced_counts.as_slice(),
        counts.as_slice()
    );
}

#[test]
fn fixed_dispersion_lrt_pipeline_validates_inputs() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
    let full =
        DesignMatrix::from_row_major(4, 2, vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0], None)
            .unwrap();
    let reduced = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();
    let same_rank = DesignMatrix::from_row_major(4, 2, vec![1.0; 8], None).unwrap();
    let rank_deficient_reduced =
        DesignMatrix::from_row_major(4, 1, vec![0.0, 0.0, 0.0, 0.0], None).unwrap();
    let bad_samples = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    assert!(DeseqBuilder::new()
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &[], 1)
        .is_err());
    assert!(DeseqBuilder::new()
        .fit_fixed_dispersion_lrt(&counts, &full, &same_rank, &[0.1], 1)
        .is_err());
    assert!(DeseqBuilder::new()
        .fit_fixed_dispersion_lrt(&counts, &full, &rank_deficient_reduced, &[0.1], 1)
        .is_err());
    assert!(DeseqBuilder::new()
        .fit_fixed_dispersion_lrt(&counts, &full, &bad_samples, &[0.1], 1)
        .is_err());
    assert!(DeseqBuilder::new()
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &[0.1], 2)
        .is_err());
}
