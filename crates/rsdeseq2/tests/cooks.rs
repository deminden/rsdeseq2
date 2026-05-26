use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

#[test]
fn samples_for_cooks_matches_model_matrix_cells() {
    let design = DesignMatrix::from_row_major(
        5,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        None,
    )
    .unwrap();

    assert_eq!(
        samples_for_cooks(&design, 3).unwrap(),
        vec![true, true, true, false, false]
    );
}

#[test]
fn robust_dispersion_uses_trimmed_cell_variance_when_replicated_cells_exist() {
    let normalized = RowMajorMatrix::from_row_major(1, 3, vec![2.0_f64, 4.0_f64, 6.0_f64]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let dispersion = robust_method_of_moments_dispersion(&normalized, &design).unwrap();

    assert_relative_eq!(dispersion[0], 0.26, epsilon = 1e-12);
}

#[test]
fn robust_dispersion_keeps_large_finite_row_mean() {
    let normalized =
        RowMajorMatrix::from_row_major(1, 3, vec![f64::MAX, f64::MAX, f64::MAX]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let dispersion = robust_method_of_moments_dispersion(&normalized, &design).unwrap();

    assert_eq!(dispersion[0], 0.04);
}

#[test]
fn cooks_distance_matches_hand_formula_for_intercept_cell() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 6]).unwrap();
    let normalized = RowMajorMatrix::from_row_major(1, 3, vec![2.0_f64, 4.0_f64, 6.0_f64]).unwrap();
    let mu = RowMajorMatrix::from_row_major(1, 3, vec![4.0_f64, 4.0_f64, 4.0_f64]).unwrap();
    let hat =
        RowMajorMatrix::from_row_major(1, 3, vec![1.0_f64 / 3.0, 1.0_f64 / 3.0, 1.0_f64 / 3.0])
            .unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let output = calculate_cooks_distance(&counts, &normalized, &mu, &hat, &design).unwrap();

    let expected = 4.0 / 8.16 * 0.75;
    assert_relative_eq!(output.robust_dispersion[0], 0.26, epsilon = 1e-12);
    assert_relative_eq!(output.cooks.row(0).unwrap()[0], expected, epsilon = 1e-12);
    assert_relative_eq!(output.cooks.row(0).unwrap()[1], 0.0, epsilon = 1e-12);
    assert_relative_eq!(output.cooks.row(0).unwrap()[2], expected, epsilon = 1e-12);
    assert_relative_eq!(output.max_cooks[0].unwrap(), expected, epsilon = 1e-12);
}

#[test]
fn cooks_distance_masks_overflowed_variance() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 6]).unwrap();
    let normalized = RowMajorMatrix::from_row_major(1, 3, vec![2.0_f64, 4.0_f64, 6.0_f64]).unwrap();
    let mu = RowMajorMatrix::from_row_major(1, 3, vec![2e154, 2e154, 2e154]).unwrap();
    let hat = RowMajorMatrix::from_row_major(1, 3, vec![0.1_f64, 0.1_f64, 0.1_f64]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let output = calculate_cooks_distance(&counts, &normalized, &mu, &hat, &design).unwrap();

    assert!(output.cooks.as_slice().iter().all(|value| value.is_nan()));
    assert_eq!(output.max_cooks, vec![None]);
}

#[test]
fn replace_outlier_counts_uses_trimmed_mean_and_size_factors() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        2,
        4,
        vec![
            10, 30, 20, 50, //
            5, 10, 15, 20,
        ],
        Some(vec!["g1".to_string(), "g2".to_string()]),
        Some(vec![
            "s1".to_string(),
            "s2".to_string(),
            "s3".to_string(),
            "s4".to_string(),
        ]),
    )
    .unwrap();
    let size_factors = vec![1.0, 2.0, 1.0, 2.0];
    let normalized = normalized_counts(&counts, &size_factors).unwrap();
    let cooks = RowMajorMatrix::from_row_major(
        2,
        4,
        vec![
            0.0, 9.0, 0.0, 0.0, //
            0.0, 0.0, 0.0, 0.0,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();

    let output = replace_outlier_counts(
        &counts,
        &normalized,
        &size_factors,
        None,
        &cooks,
        &design,
        &CooksReplacementOptions {
            trim: 0.25,
            cooks_cutoff: 5.0,
            min_replicates: 3,
            which_samples: None,
        },
    )
    .unwrap();

    assert_eq!(output.replaceable_samples, vec![true; 4]);
    assert_eq!(output.replace, vec![Some(true), Some(false)]);
    assert_eq!(
        output.candidate_replacement_counts.as_slice(),
        &[17, 35, 17, 35, 7, 15, 7, 15]
    );
    assert_eq!(
        output.replaced_counts.as_slice(),
        &[10, 35, 20, 50, 5, 10, 15, 20]
    );
    assert_eq!(
        output.outlier_cells.row(0).unwrap(),
        &[false, true, false, false]
    );
    assert_eq!(
        output.replaced_counts.gene_names().unwrap(),
        &["g1".to_string(), "g2".to_string()]
    );
    assert_eq!(
        output.replaced_counts.sample_names().unwrap(),
        &[
            "s1".to_string(),
            "s2".to_string(),
            "s3".to_string(),
            "s4".to_string()
        ]
    );
}

#[test]
fn replace_outlier_counts_rejects_nonfinite_scaled_mean() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![1, 2, 3, 4]).unwrap();
    let normalized = RowMajorMatrix::from_row_major(1, 4, vec![f64::MAX; 4]).unwrap();
    let normalization_factors =
        RowMajorMatrix::from_row_major(1, 4, vec![2.0, 1.0, 1.0, 1.0]).unwrap();
    let cooks = RowMajorMatrix::from_row_major(1, 4, vec![9.0, 0.0, 0.0, 0.0]).unwrap();
    let design = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();

    let err = replace_outlier_counts(
        &counts,
        &normalized,
        &[1.0; 4],
        Some(&normalization_factors),
        &cooks,
        &design,
        &CooksReplacementOptions {
            trim: 0.0,
            cooks_cutoff: 5.0,
            min_replicates: 3,
            which_samples: None,
        },
    )
    .unwrap_err();

    assert!(matches!(
        err,
        DeseqError::NonFiniteValue { context, index, .. }
            if context == "replacement scaled mean" && index == Some(0)
    ));
}

#[test]
fn replace_outlier_counts_uses_normalization_factors_when_supplied() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 20, 5, 80]).unwrap();
    let factors = RowMajorMatrix::from_row_major(1, 4, vec![1.0, 2.0, 4.0, 8.0]).unwrap();
    let normalized = normalized_counts_with_factors(&counts, &factors).unwrap();
    let cooks = RowMajorMatrix::from_row_major(1, 4, vec![0.0, 0.0, 9.0, 0.0]).unwrap();
    let design = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();

    let output = replace_outlier_counts(
        &counts,
        &normalized,
        &[],
        Some(&factors),
        &cooks,
        &design,
        &CooksReplacementOptions {
            trim: 0.25,
            cooks_cutoff: 5.0,
            min_replicates: 3,
            which_samples: None,
        },
    )
    .unwrap();

    assert_eq!(
        output.candidate_replacement_counts.as_slice(),
        &[10, 20, 40, 80]
    );
    assert_eq!(output.replaced_counts.as_slice(), &[10, 20, 40, 80]);
}

#[test]
fn replace_outlier_counts_respects_replaceable_samples_and_missing_cooks() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        4,
        vec![
            10, 30, 20, 50, //
            5, 10, 15, 20,
        ],
    )
    .unwrap();
    let size_factors = vec![1.0, 2.0, 1.0, 2.0];
    let normalized = normalized_counts(&counts, &size_factors).unwrap();
    let cooks = RowMajorMatrix::from_row_major(
        2,
        4,
        vec![
            0.0,
            9.0,
            0.0,
            0.0, //
            f64::NAN,
            0.0,
            0.0,
            0.0,
        ],
    )
    .unwrap();
    let design =
        DesignMatrix::from_row_major(4, 2, vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0], None)
            .unwrap();

    let output = replace_outlier_counts(
        &counts,
        &normalized,
        &size_factors,
        None,
        &cooks,
        &design,
        &CooksReplacementOptions {
            trim: 0.25,
            cooks_cutoff: 5.0,
            min_replicates: 3,
            which_samples: Some(vec![false, false, true, true]),
        },
    )
    .unwrap();

    assert_eq!(output.replaceable_samples, vec![false, false, true, true]);
    assert_eq!(output.replace, vec![Some(true), None]);
    assert_eq!(output.replaced_counts.as_slice(), counts.as_slice());
}

#[test]
fn replace_outlier_counts_validates_options_and_dimensions() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![1, 2, 3]).unwrap();
    let normalized = normalized_counts(&counts, &[1.0, 1.0, 1.0]).unwrap();
    let cooks = RowMajorMatrix::from_row_major(1, 3, vec![0.0, 0.0, 0.0]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    assert!(replace_outlier_counts(
        &counts,
        &normalized,
        &[1.0, 1.0, 1.0],
        None,
        &cooks,
        &design,
        &CooksReplacementOptions {
            trim: 0.2,
            cooks_cutoff: 1.0,
            min_replicates: 2,
            which_samples: None,
        },
    )
    .is_err());
    assert!(replace_outlier_counts(
        &counts,
        &normalized,
        &[1.0, 1.0, 1.0],
        None,
        &cooks,
        &design,
        &CooksReplacementOptions {
            trim: 0.2,
            cooks_cutoff: f64::NAN,
            min_replicates: 3,
            which_samples: None,
        },
    )
    .is_err());
    assert!(replace_outlier_counts(
        &counts,
        &normalized,
        &[1.0, 1.0, 1.0],
        None,
        &cooks,
        &design,
        &CooksReplacementOptions {
            trim: 0.2,
            cooks_cutoff: 1.0,
            min_replicates: 3,
            which_samples: Some(vec![true, false]),
        },
    )
    .is_err());
}

#[test]
fn replace_outlier_counts_skips_when_samples_do_not_exceed_coefficients() {
    let counts = CountMatrix::from_row_major_u32(1, 2, vec![10, 100]).unwrap();
    let normalized = normalized_counts(&counts, &[1.0, 1.0]).unwrap();
    let cooks = RowMajorMatrix::from_row_major(1, 2, vec![10.0, 10.0]).unwrap();
    let saturated = DesignMatrix::from_row_major(2, 2, vec![1.0, 0.0, 0.0, 1.0], None).unwrap();

    let output = replace_outlier_counts(
        &counts,
        &normalized,
        &[1.0, 1.0],
        None,
        &cooks,
        &saturated,
        &CooksReplacementOptions {
            trim: 0.2,
            cooks_cutoff: 1.0,
            min_replicates: 3,
            which_samples: None,
        },
    )
    .unwrap();

    assert_eq!(output.replaced_counts.as_slice(), counts.as_slice());
    assert_eq!(
        output.candidate_replacement_counts.as_slice(),
        counts.as_slice()
    );
    assert_eq!(output.outlier_cells.as_slice(), &[false, false]);
    assert_eq!(output.replace, vec![Some(false)]);
    assert_eq!(output.replaceable_samples, vec![false, false]);
}

#[test]
fn original_infinite_min_replicates_disables_outlier_replacement() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![100_000, 10, 10, 10]).unwrap();
    let normalized = normalized_counts(&counts, &[1.0, 1.0, 1.0, 1.0]).unwrap();
    let cooks = RowMajorMatrix::from_row_major(1, 4, vec![100.0, 0.0, 0.0, 0.0]).unwrap();
    let design = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();

    let output = replace_outlier_counts(
        &counts,
        &normalized,
        &[1.0, 1.0, 1.0, 1.0],
        None,
        &cooks,
        &design,
        &CooksReplacementOptions {
            trim: 0.2,
            cooks_cutoff: 1.0,
            min_replicates: usize::MAX,
            which_samples: None,
        },
    )
    .unwrap();

    assert_eq!(output.replaced_counts.as_slice(), counts.as_slice());
    assert_eq!(
        output.outlier_cells.as_slice(),
        &[true, false, false, false]
    );
    assert_eq!(output.replace, vec![Some(true)]);
    assert_eq!(output.replaceable_samples, vec![false, false, false, false]);
    assert_eq!(
        output.candidate_replacement_counts.as_slice(),
        &[25007, 25007, 25007, 25007]
    );
}

#[test]
fn prepare_cooks_replacement_refit_identifies_refit_rows_and_base_metadata() {
    let counts = CountMatrix::from_row_major_u32(
        3,
        4,
        vec![
            10, 30, 20, 50, //
            0, 0, 0, 0, //
            5, 10, 15, 20,
        ],
    )
    .unwrap();
    let size_factors = vec![1.0, 1.0, 1.0, 1.0];
    let normalized = normalized_counts(&counts, &size_factors).unwrap();
    let cooks = RowMajorMatrix::from_row_major(
        3,
        4,
        vec![
            0.0, 9.0, 0.0, 0.0, //
            8.0, 0.0, 0.0, 0.0, //
            0.0, 0.0, 0.0, 0.0,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();

    let plan = prepare_cooks_replacement_refit(
        &counts,
        &normalized,
        &size_factors,
        None,
        &cooks,
        &design,
        &CooksReplacementOptions {
            trim: 0.25,
            cooks_cutoff: 5.0,
            min_replicates: 3,
            which_samples: None,
        },
    )
    .unwrap();

    assert_eq!(plan.n_refit, 2);
    assert_eq!(plan.refit_rows, vec![0]);
    assert_eq!(plan.new_all_zero_rows, vec![1]);
    assert!(plan.should_refit);
    assert_eq!(
        plan.replacement.replaced_counts.as_slice(),
        &[10, 25, 20, 50, 0, 0, 0, 0, 5, 10, 15, 20]
    );
    assert_relative_eq!(plan.replaced_base_mean[0], 26.25, epsilon = 1e-12);
    assert_relative_eq!(plan.replaced_base_mean[1], 0.0, epsilon = 1e-12);
    assert_eq!(plan.replaced_all_zero, vec![false, true, false]);
    assert_eq!(plan.post_refit_max_cooks, vec![None, None, None]);
}

#[test]
fn max_cooks_after_replacement_refit_zeros_replaceable_columns() {
    let cooks = RowMajorMatrix::from_row_major(1, 6, vec![9.0, 2.0, 3.0, 8.0, 4.0, 5.0]).unwrap();
    let design =
        DesignMatrix::from_row_major(6, 1, vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0], None).unwrap();

    let max_cooks = max_cooks_after_replacement_refit(
        &cooks,
        &design,
        &[true, false, false, true, false, false],
    )
    .unwrap();

    assert_eq!(max_cooks, vec![Some(5.0)]);
    assert!(max_cooks_after_replacement_refit(&cooks, &design, &[true, true]).is_err());
}

#[test]
fn record_max_cooks_returns_missing_without_replicated_cells() {
    let cooks = RowMajorMatrix::from_row_major(1, 2, vec![0.2_f64, 0.3_f64]).unwrap();
    let design = DesignMatrix::from_row_major(2, 1, vec![1.0, 1.0], None).unwrap();
    let samples = samples_for_cooks(&design, 3).unwrap();

    assert_eq!(
        record_max_cooks(&cooks, &design, &samples).unwrap(),
        vec![None]
    );
}

#[test]
fn cooks_distance_validates_dimensions() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 6]).unwrap();
    let normalized = RowMajorMatrix::from_row_major(1, 2, vec![2.0_f64, 4.0_f64]).unwrap();
    let mu = RowMajorMatrix::from_row_major(1, 3, vec![4.0_f64, 4.0_f64, 4.0_f64]).unwrap();
    let hat =
        RowMajorMatrix::from_row_major(1, 3, vec![1.0_f64 / 3.0, 1.0_f64 / 3.0, 1.0_f64 / 3.0])
            .unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    assert!(calculate_cooks_distance(&counts, &normalized, &mu, &hat, &design).is_err());
}
