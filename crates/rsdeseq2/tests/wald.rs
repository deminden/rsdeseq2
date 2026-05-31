use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

fn toy_fit(beta: Vec<f64>, beta_se: Vec<f64>, n_genes: usize, n_coef: usize) -> NbinomGlmFit {
    let n_samples = 2;
    NbinomGlmFit {
        log_like: vec![0.0; n_genes],
        beta_converged: vec![true; n_genes],
        beta: RowMajorMatrix::from_row_major(n_genes, n_coef, beta).unwrap(),
        beta_se: RowMajorMatrix::from_row_major(n_genes, n_coef, beta_se).unwrap(),
        beta_optim_start: RowMajorMatrix::from_elem(n_genes, n_coef, f64::NAN).unwrap(),
        beta_covariance: None,
        mu: RowMajorMatrix::from_row_major(n_genes, n_samples, vec![1.0; n_genes * n_samples])
            .unwrap(),
        beta_iter: vec![1; n_genes],
        beta_optim_iter: vec![f64::NAN; n_genes],
        beta_optim_start_objective: vec![f64::NAN; n_genes],
        beta_optim_objective: vec![f64::NAN; n_genes],
        beta_optim_gradient_norm: vec![f64::NAN; n_genes],
        model_matrix: DesignMatrix::from_row_major(
            n_samples,
            n_coef,
            vec![1.0; n_samples * n_coef],
            None,
        )
        .unwrap(),
        n_terms: n_coef,
        hat_diagonal: RowMajorMatrix::from_row_major(
            n_genes,
            n_samples,
            vec![0.5; n_genes * n_samples],
        )
        .unwrap(),
    }
}

#[test]
fn wald_stat_and_pvalue_match_default_deseq2_normal_path() {
    let (stat, pvalue) = wald_stat_and_pvalue(2.0, 0.5).unwrap();
    assert_relative_eq!(stat, 4.0, epsilon = 1e-12);
    assert_relative_eq!(pvalue, 6.334248366623996e-5, epsilon = 1e-14);
    assert_relative_eq!(pvalue, two_sided_normal_pvalue(-4.0), epsilon = 1e-15);
}

#[test]
fn two_sided_t_pvalue_matches_r_pt_shape() {
    assert_relative_eq!(
        two_sided_t_pvalue(2.0, 10.0).unwrap(),
        0.07338803477074039,
        epsilon = 1e-11
    );
    assert_eq!(two_sided_t_pvalue(2.0, 0.0), None);
    assert_eq!(two_sided_t_pvalue(2.0, f64::NAN), None);
}

#[test]
fn two_sided_t_pvalue_uses_stable_upper_tail() {
    let pvalue = two_sided_t_pvalue(1.0e6, 10.0).unwrap();

    assert!(pvalue.is_finite());
    assert!(pvalue > 0.0);
    assert!(pvalue < 1.0e-50);
}

#[test]
fn wald_probability_helpers_bound_or_reject_public_pvalues() {
    assert_eq!(two_sided_normal_pvalue(f64::INFINITY), 0.0);
    assert_eq!(two_sided_normal_pvalue(f64::NAN), 1.0);
    assert_eq!(two_sided_t_pvalue(f64::INFINITY, 10.0), None);
}

#[test]
fn wald_test_coefficient_selects_requested_column() {
    let fit = toy_fit(vec![1.0, 2.0, 3.0, 4.0], vec![1.0, 0.5, 1.5, 2.0], 2, 2);
    let wald = wald_test_coefficient(&fit, 1).unwrap();
    assert_relative_eq!(wald.stat[0].unwrap(), 4.0, epsilon = 1e-12);
    assert_relative_eq!(wald.stat[1].unwrap(), 2.0, epsilon = 1e-12);
    assert_relative_eq!(
        wald.pvalue[0].unwrap(),
        two_sided_normal_pvalue(4.0),
        epsilon = 1e-11
    );
    assert_relative_eq!(
        wald.pvalue[1].unwrap(),
        two_sided_normal_pvalue(2.0),
        epsilon = 1e-11
    );
}

#[test]
fn wald_test_wraps_selected_coefficient_with_options() {
    let fit = toy_fit(vec![1.0, 2.0, 3.0, 4.0], vec![1.0, 0.5, 1.5, 2.0], 2, 2);
    let options = WaldTestOptions::normal().with_lfc_threshold(1.0, WaldAlternative::Greater);

    let wrapped = wald_test(&fit, 1, &options).unwrap();
    let direct = wald_test_coefficient_with_options(&fit, 1, &options).unwrap();

    assert_eq!(wrapped, direct);
}

#[test]
fn wald_test_contrast_uses_full_beta_covariance() {
    let mut fit = toy_fit(vec![1.0, 3.0], vec![0.5, 0.7], 1, 2);
    fit.beta_covariance =
        Some(RowMajorMatrix::from_row_major(1, 4, vec![0.25, 0.05, 0.05, 0.49]).unwrap());

    let contrast = wald_test_contrast(&fit, &[-1.0, 1.0]).unwrap();

    assert_relative_eq!(contrast.log2_fold_change[0].unwrap(), 2.0, epsilon = 1e-12);
    assert_relative_eq!(contrast.lfc_se[0].unwrap(), 0.8, epsilon = 1e-12);
    assert_relative_eq!(contrast.wald.stat[0].unwrap(), 2.5, epsilon = 1e-12);
    assert_relative_eq!(
        contrast.wald.pvalue[0].unwrap(),
        two_sided_normal_pvalue(2.5),
        epsilon = 1e-12
    );
}

#[test]
fn wald_test_contrast_masks_nonfinite_estimate_or_variance_accumulation() {
    let mut fit = toy_fit(
        vec![f64::MAX, f64::MAX, 1.0, 3.0],
        vec![0.5, 0.7, 0.5, 0.7],
        2,
        2,
    );
    fit.beta_covariance = Some(
        RowMajorMatrix::from_row_major(
            2,
            4,
            vec![0.25, 0.05, 0.05, 0.49, f64::MAX, 0.0, 0.0, f64::MAX],
        )
        .unwrap(),
    );

    let contrast = wald_test_contrast(&fit, &[1.0, 1.0]).unwrap();

    assert_eq!(contrast.log2_fold_change, vec![None, None]);
    assert_eq!(contrast.lfc_se, vec![None, None]);
    assert_eq!(contrast.wald.stat, vec![None, None]);
    assert_eq!(contrast.wald.pvalue, vec![None, None]);
}

#[test]
fn wald_test_contrast_keeps_large_cancelling_accumulations_finite() {
    let mut fit = toy_fit(
        vec![f64::MAX, f64::MAX, f64::MAX, -f64::MAX],
        vec![0.5, 0.7, 0.5, 0.7],
        2,
        2,
    );
    fit.beta_covariance = Some(
        RowMajorMatrix::from_row_major(2, 4, vec![0.25, 0.05, 0.05, 0.49, 0.25, 0.05, 0.05, 0.49])
            .unwrap(),
    );

    let contrast = wald_test_contrast(&fit, &[1.0, 1.0]).unwrap();

    assert_eq!(contrast.log2_fold_change[0], None);
    assert_eq!(contrast.log2_fold_change[1], Some(0.0));
    assert_eq!(contrast.lfc_se[0], None);
    assert_relative_eq!(
        contrast.lfc_se[1].unwrap(),
        0.916515138991168,
        epsilon = 1e-12
    );
    assert_eq!(contrast.wald.stat[1], Some(0.0));
    assert_eq!(contrast.wald.pvalue[1], Some(1.0));
}

#[test]
fn wald_test_contrast_orders_variance_products_by_magnitude() {
    let mut fit = toy_fit(vec![0.0, 1.0e200], vec![0.5, 0.7], 1, 2);
    fit.beta_covariance =
        Some(RowMajorMatrix::from_row_major(1, 4, vec![0.0, 1.0e200, 1.0e200, 0.0]).unwrap());

    let contrast = wald_test_contrast(&fit, &[1.0e200, 1.0e-200]).unwrap();

    assert_eq!(contrast.log2_fold_change[0], Some(1.0));
    assert_relative_eq!(
        contrast.lfc_se[0].unwrap(),
        (2.0e200_f64).sqrt(),
        max_relative = 1e-12
    );
    assert!(contrast.wald.stat[0].unwrap().is_finite());
    assert!(contrast.wald.pvalue[0].unwrap().is_finite());
}

#[test]
fn wald_test_contrast_validates_inputs() {
    let fit = toy_fit(vec![1.0, 3.0], vec![0.5, 0.7], 1, 2);

    assert!(wald_test_contrast(&fit, &[-1.0, 1.0]).is_err());
    assert!(wald_test_contrast(&fit, &[1.0]).is_err());
    assert!(wald_test_contrast(&fit, &[0.0, 0.0]).is_err());
    assert!(wald_test_contrast(&fit, &[1.0, f64::NAN]).is_err());
}

#[test]
fn wald_test_coefficient_can_use_scalar_t_degrees_of_freedom() {
    let fit = toy_fit(vec![2.0, 1.0], vec![0.5, 1.0], 2, 1);
    let wald =
        wald_test_coefficient_with_options(&fit, 0, &WaldTestOptions::t_degrees_of_freedom(10.0))
            .unwrap();

    assert_eq!(
        wald.degrees_of_freedom.as_ref().unwrap(),
        &vec![Some(10.0), Some(10.0)]
    );
    assert_relative_eq!(wald.stat[0].unwrap(), 4.0, epsilon = 1e-12);
    assert_relative_eq!(
        wald.pvalue[0].unwrap(),
        two_sided_t_pvalue(4.0, 10.0).unwrap(),
        epsilon = 1e-11
    );
    assert!(wald.pvalue[0].unwrap() > two_sided_normal_pvalue(4.0));
}

#[test]
fn original_use_t_vector_degrees_of_freedom_drive_pvalues() {
    let fit = toy_fit(vec![2.0, -2.0], vec![0.5, 0.5], 2, 1);
    let wald = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::t_per_gene_degrees_of_freedom(vec![12.0, 6.0]),
    )
    .unwrap();

    assert_eq!(
        wald.degrees_of_freedom.as_ref().unwrap(),
        &vec![Some(12.0), Some(6.0)]
    );
    assert_relative_eq!(wald.stat[0].unwrap(), 4.0, epsilon = 1e-12);
    assert_relative_eq!(wald.stat[1].unwrap(), -4.0, epsilon = 1e-12);
    assert_relative_eq!(
        wald.pvalue[0].unwrap(),
        two_sided_t_pvalue(4.0, 12.0).unwrap(),
        epsilon = 1e-15
    );
    assert_relative_eq!(
        wald.pvalue[1].unwrap(),
        two_sided_t_pvalue(-4.0, 6.0).unwrap(),
        epsilon = 1e-15
    );
}

#[test]
fn wald_test_coefficient_can_use_greater_abs_lfc_threshold() {
    let fit = toy_fit(vec![2.0], vec![0.5], 1, 1);
    let wald = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::normal().with_lfc_threshold(1.0, WaldAlternative::GreaterAbs),
    )
    .unwrap();

    assert_relative_eq!(wald.stat[0].unwrap(), 4.0, epsilon = 1e-12);
    assert_relative_eq!(
        wald.pvalue[0].unwrap(),
        0.02275013293476686,
        epsilon = 1e-11
    );
}

#[test]
fn wald_test_coefficient_can_use_older_greater_abs_lfc_threshold() {
    let fit = toy_fit(vec![2.0], vec![0.5], 1, 1);
    let wald = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::normal().with_lfc_threshold(1.0, WaldAlternative::GreaterAbs2014),
    )
    .unwrap();

    assert_relative_eq!(wald.stat[0].unwrap(), 2.0, epsilon = 1e-12);
    assert_relative_eq!(
        wald.pvalue[0].unwrap(),
        0.04550026389635842,
        epsilon = 1e-11
    );
}

#[test]
fn wald_test_coefficient_upper_tail_keeps_extreme_threshold_pvalue_finite() {
    let fit = toy_fit(vec![1.0e6], vec![1.0], 1, 1);
    let wald = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::t_degrees_of_freedom(10.0)
            .with_lfc_threshold(1.0, WaldAlternative::Greater),
    )
    .unwrap();

    let pvalue = wald.pvalue[0].unwrap();
    assert!(pvalue.is_finite());
    assert!(pvalue > 0.0);
    assert!(pvalue < 1.0e-50);
}

#[test]
fn original_greater_abs_upshot_matches_greater_abs_at_zero_threshold() {
    let fit = toy_fit(vec![2.0, -1.5], vec![0.5, 0.75], 2, 1);
    let upshot = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::normal().with_lfc_threshold(0.0, WaldAlternative::GreaterAbsUpshot),
    )
    .unwrap();
    let greater_abs = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::normal().with_lfc_threshold(0.0, WaldAlternative::GreaterAbs),
    )
    .unwrap();

    assert_eq!(upshot.stat, greater_abs.stat);
    assert_eq!(upshot.pvalue, greater_abs.pvalue);
    assert_eq!(upshot.pvalue.len(), 2);
}

#[test]
fn wald_test_coefficient_upshot_masks_overflowed_threshold_formula() {
    let fit = toy_fit(vec![f64::MAX], vec![1.0], 1, 1);
    let wald = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::normal().with_lfc_threshold(1.0, WaldAlternative::GreaterAbsUpshot),
    )
    .unwrap();

    assert_eq!(wald.stat[0], Some(f64::MAX));
    assert_eq!(wald.pvalue[0], None);
}

#[test]
fn wald_test_coefficient_can_use_less_abs_lfc_threshold() {
    let fit = toy_fit(vec![0.2], vec![0.5], 1, 1);
    let wald = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::normal().with_lfc_threshold(1.0, WaldAlternative::LessAbs),
    )
    .unwrap();

    assert_relative_eq!(wald.stat[0].unwrap(), 1.6, epsilon = 1e-12);
    assert_relative_eq!(
        wald.pvalue[0].unwrap(),
        0.054799291699557995,
        epsilon = 1e-11
    );
}

#[test]
fn wald_test_coefficient_can_use_one_sided_thresholds() {
    let greater_fit = toy_fit(vec![2.0], vec![0.5], 1, 1);
    let greater = wald_test_coefficient_with_options(
        &greater_fit,
        0,
        &WaldTestOptions::normal().with_lfc_threshold(1.0, WaldAlternative::Greater),
    )
    .unwrap();
    assert_relative_eq!(greater.stat[0].unwrap(), 2.0, epsilon = 1e-12);
    assert_relative_eq!(
        greater.pvalue[0].unwrap(),
        0.02275013194817921,
        epsilon = 1e-11
    );

    let less_fit = toy_fit(vec![-2.0], vec![0.5], 1, 1);
    let less = wald_test_coefficient_with_options(
        &less_fit,
        0,
        &WaldTestOptions::normal().with_lfc_threshold(1.0, WaldAlternative::Less),
    )
    .unwrap();
    assert_relative_eq!(less.stat[0].unwrap(), -2.0, epsilon = 1e-12);
    assert_relative_eq!(
        less.pvalue[0].unwrap(),
        0.02275013194817921,
        epsilon = 1e-11
    );

    let opposite_fit = toy_fit(vec![2.0], vec![0.5], 1, 1);
    let opposite = wald_test_coefficient_with_options(
        &opposite_fit,
        0,
        &WaldTestOptions::normal().with_lfc_threshold(1.0, WaldAlternative::Less),
    )
    .unwrap();
    assert_eq!(opposite.stat[0], Some(0.0));
    assert!(opposite.pvalue[0].unwrap() > 0.999);
}

#[test]
fn wald_test_coefficient_can_use_t_lfc_threshold() {
    let fit = toy_fit(vec![2.0], vec![0.5], 1, 1);
    let wald = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::t_degrees_of_freedom(10.0)
            .with_lfc_threshold(1.0, WaldAlternative::Greater),
    )
    .unwrap();

    assert_relative_eq!(wald.stat[0].unwrap(), 2.0, epsilon = 1e-12);
    assert_relative_eq!(
        wald.pvalue[0].unwrap(),
        two_sided_t_pvalue(2.0, 10.0).unwrap() / 2.0,
        epsilon = 1e-15
    );
}

#[test]
fn original_use_t_threshold_alternatives_use_t_tails() {
    let fit = toy_fit(vec![2.0, -2.0, 0.2, 2.0], vec![0.5, 0.5, 0.5, 0.5], 4, 1);
    let df = 11.0;

    let greater_abs_2014 = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::t_degrees_of_freedom(df)
            .with_lfc_threshold(1.0, WaldAlternative::GreaterAbs2014),
    )
    .unwrap();
    assert_relative_eq!(greater_abs_2014.stat[0].unwrap(), 2.0, epsilon = 1e-12);
    assert_relative_eq!(
        greater_abs_2014.pvalue[0].unwrap(),
        two_sided_t_pvalue(2.0, df).unwrap(),
        epsilon = 1e-15
    );

    let greater = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::t_degrees_of_freedom(df)
            .with_lfc_threshold(1.0, WaldAlternative::Greater),
    )
    .unwrap();
    assert_relative_eq!(greater.stat[0].unwrap(), 2.0, epsilon = 1e-12);
    assert_relative_eq!(
        greater.pvalue[0].unwrap(),
        two_sided_t_pvalue(2.0, df).unwrap() / 2.0,
        epsilon = 1e-15
    );

    let less = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::t_degrees_of_freedom(df).with_lfc_threshold(1.0, WaldAlternative::Less),
    )
    .unwrap();
    assert_relative_eq!(less.stat[1].unwrap(), -2.0, epsilon = 1e-12);
    assert_relative_eq!(
        less.pvalue[1].unwrap(),
        two_sided_t_pvalue(-2.0, df).unwrap() / 2.0,
        epsilon = 1e-15
    );

    let less_abs = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::t_degrees_of_freedom(df)
            .with_lfc_threshold(1.0, WaldAlternative::LessAbs),
    )
    .unwrap();
    assert_relative_eq!(less_abs.stat[2].unwrap(), 1.6, epsilon = 1e-12);
    assert_relative_eq!(
        less_abs.pvalue[2].unwrap(),
        two_sided_t_pvalue(1.6, df).unwrap() / 2.0,
        epsilon = 1e-15
    );
}

#[test]
fn original_use_t_novel_contrast_uses_t_tail() {
    let mut fit = toy_fit(vec![0.0, 1.0, 3.0, 0.0, -1.0, 2.0], vec![0.5; 6], 2, 3);
    fit.beta_covariance = Some(
        RowMajorMatrix::from_row_major(
            2,
            9,
            vec![
                0.25, 0.02, 0.01, 0.02, 0.25, 0.03, 0.01, 0.03, 0.25, //
                0.25, 0.00, 0.00, 0.00, 0.25, 0.00, 0.00, 0.00, 0.25,
            ],
        )
        .unwrap(),
    );

    let contrast = vec![0.0, -1.0, 1.0];
    let wald = wald_test_contrast_with_options(
        &fit,
        &contrast,
        &WaldTestOptions::t_degrees_of_freedom(11.0),
    )
    .unwrap();

    assert_relative_eq!(wald.log2_fold_change[0].unwrap(), 2.0, epsilon = 1e-12);
    assert_relative_eq!(wald.lfc_se[0].unwrap(), (0.44_f64).sqrt(), epsilon = 1e-12);
    assert_relative_eq!(
        wald.wald.pvalue[0].unwrap(),
        two_sided_t_pvalue(wald.wald.stat[0].unwrap(), 11.0).unwrap(),
        epsilon = 1e-15
    );
    assert_relative_eq!(wald.log2_fold_change[1].unwrap(), 3.0, epsilon = 1e-12);
    assert_relative_eq!(
        wald.wald.pvalue[1].unwrap(),
        two_sided_t_pvalue(wald.wald.stat[1].unwrap(), 11.0).unwrap(),
        epsilon = 1e-15
    );
}

#[test]
fn wald_test_coefficient_rejects_invalid_threshold_options() {
    let fit = toy_fit(vec![2.0], vec![0.5], 1, 1);
    assert!(wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::normal().with_lfc_threshold(-1.0, WaldAlternative::Greater)
    )
    .is_err());
    assert!(wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::normal().with_lfc_threshold(0.0, WaldAlternative::LessAbs)
    )
    .is_err());
    assert!(wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::t_degrees_of_freedom(10.0)
            .with_lfc_threshold(1.0, WaldAlternative::GreaterAbsUpshot)
    )
    .is_err());
}

#[test]
fn wald_test_coefficient_can_use_residual_t_degrees_of_freedom() {
    let fit = toy_fit(vec![2.0], vec![1.0], 1, 1);
    let wald = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::t_residual_degrees_of_freedom(),
    )
    .unwrap();

    assert_eq!(wald.degrees_of_freedom.as_ref().unwrap(), &vec![Some(1.0)]);
    assert_relative_eq!(
        wald.pvalue[0].unwrap(),
        two_sided_t_pvalue(2.0, 1.0).unwrap(),
        epsilon = 1e-15
    );
}

#[test]
fn wald_test_coefficient_preserves_statistic_when_t_df_is_missing() {
    let fit = toy_fit(vec![2.0, 2.0, 2.0], vec![1.0, 1.0, 1.0], 3, 1);
    let wald = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::t_per_gene_degrees_of_freedom(vec![10.0, 0.0, f64::NAN]),
    )
    .unwrap();

    assert_eq!(
        wald.degrees_of_freedom.as_ref().unwrap(),
        &vec![Some(10.0), None, None]
    );
    assert_eq!(wald.stat, vec![Some(2.0), Some(2.0), Some(2.0)]);
    assert!(wald.pvalue[0].is_some());
    assert_eq!(wald.pvalue[1], None);
    assert_eq!(wald.pvalue[2], None);
}

#[test]
fn wald_test_coefficient_validates_per_gene_degrees_of_freedom_length() {
    let fit = toy_fit(vec![2.0, 2.0], vec![1.0, 1.0], 2, 1);
    let err = wald_test_coefficient_with_options(
        &fit,
        0,
        &WaldTestOptions::t_per_gene_degrees_of_freedom(vec![10.0]),
    );

    assert!(err.is_err());
}

#[test]
fn wald_test_returns_none_for_missing_or_invalid_se() {
    let fit = toy_fit(
        vec![1.0, 2.0, f64::NAN, 4.0],
        vec![1.0, 0.0, 1.0, f64::NAN],
        2,
        2,
    );
    let wald = wald_test_coefficient(&fit, 1).unwrap();
    assert_eq!(wald.stat[0], None);
    assert_eq!(wald.pvalue[0], None);
    assert_eq!(wald.stat[1], None);
    assert_eq!(wald.pvalue[1], None);
}

#[test]
fn wald_test_validates_coefficient_index() {
    let fit = toy_fit(vec![1.0, 2.0], vec![1.0, 1.0], 1, 2);
    assert!(wald_test_coefficient(&fit, 2).is_err());
}

#[test]
fn wald_pvalues_compose_with_bh_adjustment() {
    let fit = toy_fit(vec![2.0, 1.0, 0.0], vec![0.5, 1.0, 1.0], 3, 1);
    let wald = wald_test_coefficient(&fit, 0).unwrap();
    let padj = bh_adjust(&wald.pvalue);
    assert_eq!(padj.len(), 3);
    assert!(padj[0].unwrap() <= padj[1].unwrap());
    assert!(padj[1].unwrap() <= padj[2].unwrap());
    assert_relative_eq!(wald.pvalue[2].unwrap(), 1.0, epsilon = 1e-15);
}
