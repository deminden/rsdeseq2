use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

fn two_group_design() -> DesignMatrix {
    DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap()
}

#[test]
fn normalized_counts_and_base_mean() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![2, 4, 8, 4, 8, 16]).unwrap();
    let normalized = normalized_counts(&counts, &[0.5, 1.0, 2.0]).unwrap();
    assert_eq!(normalized.as_slice(), &[4.0, 4.0, 4.0, 8.0, 8.0, 8.0]);

    let means = base_mean(&normalized).unwrap();
    assert_relative_eq!(means[0], 4.0, epsilon = 1e-12);
    assert_relative_eq!(means[1], 8.0, epsilon = 1e-12);

    let variances = base_variance(&normalized).unwrap();
    assert_relative_eq!(variances[0], 0.0, epsilon = 1e-12);
    assert_relative_eq!(variances[1], 0.0, epsilon = 1e-12);
}

#[test]
fn normalized_counts_reject_overflowed_tiny_size_factor_division() {
    let counts = CountMatrix::from_row_major_u32(1, 1, vec![u32::MAX]).unwrap();
    let err = normalized_counts(&counts, &[f64::MIN_POSITIVE]).unwrap_err();

    assert!(matches!(
        err,
        DeseqError::NonFiniteValue { context, index, .. }
            if context == "normalized count" && index == Some(0)
    ));
}

#[test]
fn normalized_counts_with_gene_sample_factors_and_base_mean() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![10, 20, 30, 6, 12, 24]).unwrap();
    let normalization_factors =
        RowMajorMatrix::from_row_major(2, 3, vec![1.0, 2.0, 5.0, 2.0, 3.0, 4.0]).unwrap();

    let normalized = normalized_counts_with_factors(&counts, &normalization_factors).unwrap();

    assert_eq!(normalized.as_slice(), &[10.0, 10.0, 6.0, 3.0, 4.0, 6.0]);
    let means = base_mean(&normalized).unwrap();
    assert_relative_eq!(means[0], 26.0 / 3.0, epsilon = 1e-12);
    assert_relative_eq!(means[1], 13.0 / 3.0, epsilon = 1e-12);
}

#[test]
fn normalized_counts_with_factors_reject_overflowed_tiny_factor_division() {
    let counts = CountMatrix::from_row_major_u32(1, 1, vec![u32::MAX]).unwrap();
    let normalization_factors =
        RowMajorMatrix::from_row_major(1, 1, vec![f64::MIN_POSITIVE]).unwrap();
    let err = normalized_counts_with_factors(&counts, &normalization_factors).unwrap_err();

    assert!(matches!(
        err,
        DeseqError::NonFiniteValue { context, index, .. }
            if context == "normalization-factor normalized count" && index == Some(0)
    ));
}

#[test]
fn normalization_factors_validate_dimensions_and_positive_values() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![10, 20, 30]).unwrap();
    let bad_dims = RowMajorMatrix::from_row_major(1, 2, vec![1.0, 2.0]).unwrap();
    let zero_factor = RowMajorMatrix::from_row_major(1, 3, vec![1.0, 0.0, 2.0]).unwrap();
    let negative_factor = RowMajorMatrix::from_row_major(1, 3, vec![1.0, -1.0, 2.0]).unwrap();
    let nan_factor = RowMajorMatrix::from_row_major(1, 3, vec![1.0, f64::NAN, 2.0]).unwrap();

    assert!(matches!(
        normalized_counts_with_factors(&counts, &bad_dims).unwrap_err(),
        DeseqError::InvalidDimensions { .. }
    ));
    assert!(matches!(
        normalized_counts_with_factors(&counts, &zero_factor).unwrap_err(),
        DeseqError::InvalidSizeFactors { .. }
    ));
    assert!(matches!(
        normalized_counts_with_factors(&counts, &negative_factor).unwrap_err(),
        DeseqError::InvalidSizeFactors { .. }
    ));
    assert!(matches!(
        normalized_counts_with_factors(&counts, &nan_factor).unwrap_err(),
        DeseqError::InvalidSizeFactors { .. }
    ));
}

#[test]
fn builder_returns_initial_fit_state() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![2, 4, 8, 4, 8, 16]).unwrap();
    let fit = DeseqBuilder::new()
        .size_factor_method(SizeFactorMethod::Ratio)
        .execution_mode(ExecutionMode::Strict)
        .fit_size_factors_and_base_means(&counts)
        .unwrap();

    assert_eq!(fit.counts_summary.n_genes, 2);
    assert_eq!(fit.counts_summary.n_samples, 3);
    assert!(fit.dispersion.is_none());
    assert!(fit.beta.is_none());
    assert_eq!(fit.base_mean, vec![4.0, 8.0]);
    assert_relative_eq!(fit.base_var[0], 0.0, epsilon = 1e-24);
    assert_relative_eq!(fit.base_var[1], 0.0, epsilon = 1e-24);
    assert_eq!(fit.all_zero, vec![false, false]);
}

#[test]
fn builder_uses_normalization_factors_for_initial_fit_state() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![10, 20, 30, 6, 12, 24]).unwrap();
    let normalization_factors =
        RowMajorMatrix::from_row_major(2, 3, vec![1.0, 2.0, 5.0, 2.0, 3.0, 4.0]).unwrap();

    let fit = DeseqBuilder::new()
        .size_factors(vec![100.0, 100.0, 100.0])
        .normalization_factors(normalization_factors.clone())
        .fit_size_factors_and_base_means(&counts)
        .unwrap();

    assert_eq!(fit.size_factors, vec![100.0, 100.0, 100.0]);
    assert_eq!(fit.normalization_factors, Some(normalization_factors));
    assert_relative_eq!(fit.base_mean[0], 26.0 / 3.0, epsilon = 1e-12);
    assert_relative_eq!(fit.base_mean[1], 13.0 / 3.0, epsilon = 1e-12);
}

#[test]
fn builder_with_normalization_factors_does_not_require_estimable_size_factors() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![0, 0, 0]).unwrap();
    let normalization_factors = RowMajorMatrix::from_row_major(1, 3, vec![1.0, 1.5, 2.0]).unwrap();

    let fit = DeseqBuilder::new()
        .normalization_factors(normalization_factors)
        .fit_size_factors_and_base_means(&counts)
        .unwrap();

    assert_eq!(fit.size_factors, vec![1.0, 1.0, 1.0]);
    assert_eq!(fit.base_mean, vec![0.0]);
    assert_eq!(fit.all_zero, vec![true]);
}

#[test]
fn base_variance_uses_sample_variance_like_row_vars() {
    let normalized =
        RowMajorMatrix::from_row_major(2, 3, vec![1.0, 2.0, 3.0, 2.0, 2.0, 8.0]).unwrap();
    let variances = base_variance(&normalized).unwrap();
    assert_relative_eq!(variances[0], 1.0, epsilon = 1e-12);
    assert_relative_eq!(variances[1], 12.0, epsilon = 1e-12);
}

#[test]
fn base_metadata_rejects_nonfinite_accumulation() {
    let normalized = RowMajorMatrix::from_row_major(1, 2, vec![f64::MAX, f64::MAX]).unwrap();

    let means = base_mean(&normalized).unwrap();
    assert_eq!(means[0], f64::MAX);
    assert_eq!(base_variance(&normalized).unwrap()[0], 0.0);

    let signed = RowMajorMatrix::from_row_major(1, 2, vec![f64::MAX, -f64::MAX]).unwrap();
    assert_eq!(base_mean(&signed).unwrap()[0], 0.0);
    assert!(
        base_variance(&RowMajorMatrix::from_row_major(1, 2, vec![f64::MAX, 0.0]).unwrap())
            .unwrap_err()
            .to_string()
            .contains("baseVar")
    );
}

#[test]
fn base_variance_avoids_finite_sum_squares_overflow() {
    let deviation = (f64::MAX * 0.5).sqrt();
    let normalized =
        RowMajorMatrix::from_row_major(1, 4, vec![0.0, 2.0 * deviation, 0.0, 2.0 * deviation])
            .unwrap();

    let variances = base_variance(&normalized).unwrap();
    assert_relative_eq!(variances[0], f64::MAX / 3.0 * 2.0, max_relative = 1e-12);
}

#[test]
fn base_variance_rejects_overflowed_scaled_deviation() {
    let normalized = RowMajorMatrix::from_row_major(1, 2, vec![f64::MAX, -f64::MAX]).unwrap();

    assert!(base_variance(&normalized)
        .unwrap_err()
        .to_string()
        .contains("baseVar"));
}

#[test]
fn weighted_base_metadata_multiplies_normalized_counts_by_weights() {
    let normalized =
        RowMajorMatrix::from_row_major(2, 3, vec![10.0, 20.0, 30.0, 2.0, 4.0, 8.0]).unwrap();
    let weights =
        RowMajorMatrix::from_row_major(2, 3, vec![1.0, 0.5, 0.0, 0.5, 1.0, 0.25]).unwrap();

    let means = base_mean_with_weights(&normalized, &weights).unwrap();
    let variances = base_variance_with_weights(&normalized, &weights).unwrap();

    assert_relative_eq!(means[0], 20.0 / 3.0, epsilon = 1e-12);
    assert_relative_eq!(variances[0], 100.0 / 3.0, epsilon = 1e-12);
    assert_relative_eq!(means[1], 7.0 / 3.0, epsilon = 1e-12);
    assert_relative_eq!(variances[1], 7.0 / 3.0, epsilon = 1e-12);
}

#[test]
fn weighted_base_metadata_rejects_nonfinite_products() {
    let normalized = RowMajorMatrix::from_row_major(1, 2, vec![f64::MAX, 1.0]).unwrap();
    let weights = RowMajorMatrix::from_row_major(1, 2, vec![2.0, 1.0]).unwrap();

    assert!(base_mean_with_weights(&normalized, &weights)
        .unwrap_err()
        .to_string()
        .contains("weighted baseMean"));
    assert!(base_variance_with_weights(&normalized, &weights)
        .unwrap_err()
        .to_string()
        .contains("weighted baseVar"));
}

#[test]
fn weighted_base_mean_avoids_finite_sum_overflow() {
    let normalized =
        RowMajorMatrix::from_row_major(2, 2, vec![f64::MAX, f64::MAX, f64::MAX, -f64::MAX])
            .unwrap();
    let weights = RowMajorMatrix::from_row_major(2, 2, vec![1.0, 1.0, 1.0, 1.0]).unwrap();

    let means = base_mean_with_weights(&normalized, &weights).unwrap();
    assert_eq!(means[0], f64::MAX);
    assert_eq!(means[1], 0.0);
}

#[test]
fn weighted_base_variance_avoids_finite_sum_squares_overflow() {
    let deviation = (f64::MAX * 0.5).sqrt();
    let normalized =
        RowMajorMatrix::from_row_major(1, 4, vec![0.0, 2.0 * deviation, 0.0, 2.0 * deviation])
            .unwrap();
    let weights = RowMajorMatrix::from_row_major(1, 4, vec![1.0, 1.0, 1.0, 1.0]).unwrap();

    let variances = base_variance_with_weights(&normalized, &weights).unwrap();
    assert_relative_eq!(variances[0], f64::MAX / 3.0 * 2.0, max_relative = 1e-12);
}

#[test]
fn weighted_base_variance_rejects_overflowed_scaled_deviation() {
    let normalized = RowMajorMatrix::from_row_major(1, 2, vec![f64::MAX, -f64::MAX]).unwrap();
    let weights = RowMajorMatrix::from_row_major(1, 2, vec![1.0, 1.0]).unwrap();

    assert!(base_variance_with_weights(&normalized, &weights)
        .unwrap_err()
        .to_string()
        .contains("weighted baseVar"));
}

#[test]
fn weighted_base_variance_is_nan_for_one_sample() {
    let normalized = RowMajorMatrix::from_row_major(1, 1, vec![10.0]).unwrap();
    let weights = RowMajorMatrix::from_row_major(1, 1, vec![0.5]).unwrap();
    let variances = base_variance_with_weights(&normalized, &weights).unwrap();
    assert!(variances[0].is_nan());
}

#[test]
fn weighted_base_metadata_validates_inputs() {
    let normalized = RowMajorMatrix::from_row_major(1, 3, vec![10.0, 20.0, 30.0]).unwrap();
    let bad_dims = RowMajorMatrix::from_row_major(1, 2, vec![1.0, 1.0]).unwrap();
    let bad_weight = RowMajorMatrix::from_row_major(1, 3, vec![1.0, -0.5, 1.0]).unwrap();

    assert!(base_mean_with_weights(&normalized, &bad_dims).is_err());
    assert!(base_variance_with_weights(&normalized, &bad_dims).is_err());
    assert!(base_mean_with_weights(&normalized, &bad_weight).is_err());
    assert!(base_variance_with_weights(&normalized, &bad_weight).is_err());
}

#[test]
fn builder_uses_observation_weights_for_initial_base_metadata() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![10, 20, 30]).unwrap();
    let weights = RowMajorMatrix::from_row_major(1, 3, vec![1.0, 0.5, 0.0]).unwrap();

    let fit = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .observation_weights(weights.clone())
        .fit_size_factors_and_base_means(&counts)
        .unwrap();

    assert_eq!(fit.observation_weights, Some(weights));
    assert_eq!(fit.weights_fail, None);
    assert_eq!(fit.weights_design_rank, None);
    assert_eq!(fit.all_zero, vec![false]);
    assert_relative_eq!(fit.base_mean[0], 20.0 / 3.0, epsilon = 1e-12);
    assert_relative_eq!(fit.base_var[0], 100.0 / 3.0, epsilon = 1e-12);
}

#[test]
fn builder_preprocesses_observation_weights_when_design_is_supplied() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        4,
        vec![
            10, 20, 30, 40, //
            50, 60, 70, 80,
        ],
    )
    .unwrap();
    let weights = RowMajorMatrix::from_row_major(
        2,
        4,
        vec![
            1.0, 2.0, 1.0, 2.0, //
            1.0, 1.0, 0.0, 0.0,
        ],
    )
    .unwrap();

    let fit = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .observation_weights(weights)
        .fit_size_factors_and_base_means_with_design(&counts, &two_group_design())
        .unwrap();

    assert_eq!(fit.weights_design_rank, Some(2));
    assert_eq!(fit.weights_fail, Some(vec![false, true]));
    assert_eq!(fit.all_zero, vec![false, true]);
    assert_eq!(
        fit.observation_weights.as_ref().unwrap().row(0).unwrap(),
        &[0.5, 1.0, 0.5, 1.0]
    );
    assert_eq!(
        fit.observation_weights.as_ref().unwrap().row(1).unwrap(),
        &[1.0, 1.0, 0.0, 0.0]
    );
    assert_relative_eq!(fit.base_mean[0], 40.0, epsilon = 1e-12);
    assert_relative_eq!(fit.base_var[0], 2600.0 / 3.0, epsilon = 1e-12);
    assert_relative_eq!(fit.base_mean[1], 27.5, epsilon = 1e-12);
}

#[test]
fn fixed_dispersion_wald_pipeline_marks_weight_failed_rows_as_skipped() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        4,
        vec![
            10, 20, 30, 40, //
            50, 60, 70, 80,
        ],
    )
    .unwrap();
    let weights = RowMajorMatrix::from_row_major(
        2,
        4,
        vec![
            1.0, 1.0, 1.0, 1.0, //
            1.0, 1.0, 0.0, 0.0,
        ],
    )
    .unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .observation_weights(weights)
        .disable_cooks_cutoff()
        .fit_fixed_dispersion_wald(&counts, &two_group_design(), &[0.1, 0.1], 1)
        .unwrap();

    assert_eq!(fit.weights_fail, Some(vec![false, true]));
    assert_eq!(fit.all_zero, vec![false, true]);
    assert!(results.rows[0].pvalue.is_some());
    assert!(results.rows[1].pvalue.is_none());
    assert!(fit.beta.as_ref().unwrap().row(1).unwrap()[0].is_nan());
}

#[test]
fn base_variance_is_nan_for_one_sample() {
    let normalized = RowMajorMatrix::from_row_major(2, 1, vec![1.0, 2.0]).unwrap();
    let variances = base_variance(&normalized).unwrap();
    assert!(variances[0].is_nan());
    assert!(variances[1].is_nan());
}
