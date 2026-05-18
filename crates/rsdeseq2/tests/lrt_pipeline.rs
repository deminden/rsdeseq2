use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

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
    assert_eq!(fit.full_deviance.as_ref().unwrap().len(), 1);
    assert_relative_eq!(
        fit.full_deviance.as_ref().unwrap()[0],
        -2.0 * fit.log_like.as_ref().unwrap()[0],
        epsilon = 1e-12
    );
    assert_eq!(fit.reduced_log_like.as_ref().unwrap().len(), 1);
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
    assert!(fit.full_deviance.as_ref().unwrap()[0].is_nan());
    assert_relative_eq!(
        fit.full_deviance.as_ref().unwrap()[1],
        -2.0 * fit.log_like.as_ref().unwrap()[1],
        epsilon = 1e-12
    );
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
