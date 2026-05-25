use approx::assert_relative_eq;
use rsdeseq2::prelude::{
    fast_vst_eligible_count, fast_vst_subset, fast_vst_subset_indices, fast_vst_subset_matrix_rows,
    fast_vst_subset_normalized_counts, local_vst_inverse_size_factor_mean,
    local_vst_inverse_size_factor_mean_from_normalization_factors, norm_transform,
    norm_transform_value, rlog, vst, vst_local, vst_mean, vst_mean_value, vst_parametric,
    vst_parametric_value, vst_with_dispersion_trend,
    vst_with_dispersion_trend_and_normalization_factors,
    vst_with_dispersion_trend_and_size_factors, CountMatrix, DeseqError, DispersionTrendFit,
    LocalDispersionTrend, LocalDispersionTrendFit, MeanDispersionTrend, MeanDispersionTrendFit,
    ParametricDispersionTrend, ParametricDispersionTrendFit, RowMajorMatrix,
};

fn expected_mean_vst(q: f64, alpha: f64) -> f64 {
    (2.0 * (alpha * q).sqrt().asinh() - alpha.ln() - 4.0_f64.ln()) / std::f64::consts::LN_2
}

fn expected_parametric_vst(q: f64, trend: ParametricDispersionTrend) -> f64 {
    ((1.0
        + trend.extra_pois
        + 2.0 * trend.asympt_disp * q
        + 2.0 * (trend.asympt_disp * q * (1.0 + trend.extra_pois + trend.asympt_disp * q)).sqrt())
        / (4.0 * trend.asympt_disp))
        .ln()
        / std::f64::consts::LN_2
}

#[test]
fn norm_transform_applies_log2_count_plus_one() {
    let normalized =
        RowMajorMatrix::from_row_major(2, 3, vec![0.0, 1.0, 3.0, 7.0, 15.0, 31.0]).unwrap();

    let transformed = norm_transform(&normalized).unwrap();

    assert_eq!(transformed.n_rows(), 2);
    assert_eq!(transformed.n_cols(), 3);
    for (observed, count) in transformed.as_slice().iter().zip(normalized.as_slice()) {
        assert_relative_eq!(*observed, (*count + 1.0).log2(), epsilon = 1e-12);
    }
}

#[test]
fn norm_transform_value_matches_scalar_formula() {
    assert_relative_eq!(norm_transform_value(1023.0, 0).unwrap(), 10.0);
}

#[test]
fn norm_transform_value_keeps_tiny_counts_precise() {
    let tiny = 1.0e-12_f64;

    let transformed = norm_transform_value(tiny, 0).unwrap();

    assert_relative_eq!(transformed, tiny / std::f64::consts::LN_2, epsilon = 1e-24);
}

#[test]
fn norm_transform_rejects_negative_and_non_finite_counts() {
    let negative = RowMajorMatrix::from_row_major(1, 1, vec![-1.0]).unwrap();
    assert!(norm_transform(&negative).is_err());

    let non_finite = RowMajorMatrix::from_row_major(1, 1, vec![f64::NAN]).unwrap();
    assert!(norm_transform(&non_finite).is_err());
}

#[test]
fn rlog_is_explicitly_unsupported_until_regularized_log_parity_lands() {
    let err = rlog().unwrap_err();

    assert!(err.to_string().contains("regularized-log transformation"));
}

#[test]
fn fast_vst_subset_indices_match_deseq2_ordered_base_mean_rule() {
    let base_mean = vec![1.0, 12.0, 6.0, 40.0, 8.0, 100.0, 9.0, 7.0];

    let selected = fast_vst_subset_indices(&base_mean, 4).unwrap();

    assert_eq!(selected, vec![2, 4, 1, 5]);
}

#[test]
fn fast_vst_subset_indices_use_r_round_half_to_even_positions() {
    let base_mean = vec![6.0, 7.0, 8.0, 9.0, 10.0, 11.0];

    let selected = fast_vst_subset_indices(&base_mean, 3).unwrap();

    assert_eq!(selected, vec![0, 3, 5]);
}

#[test]
fn original_fast_vst_single_subset_uses_first_ordered_eligible_row() {
    let base_mean = vec![5.0, 12.0, 6.0, 6.0, 100.0];

    let selected = fast_vst_subset_indices(&base_mean, 1).unwrap();

    assert_eq!(selected, vec![2]);
}

#[test]
fn fast_vst_subset_indices_reject_invalid_inputs() {
    assert!(fast_vst_subset_indices(&[6.0, 7.0], 0).is_err());
    assert!(fast_vst_subset_indices(&[1.0, 2.0, 6.0], 2).is_err());
    assert!(fast_vst_subset_indices(&[6.0, f64::NAN, 8.0], 2).is_err());
}

#[test]
fn fast_vst_eligible_count_matches_base_mean_threshold() {
    let base_mean = vec![5.0, 5.01, 6.0, 1.0, 100.0];

    assert_eq!(fast_vst_eligible_count(&base_mean).unwrap(), 3);
    assert!(fast_vst_eligible_count(&[6.0, f64::NAN]).is_err());
}

#[test]
fn fast_vst_subset_normalized_counts_returns_selected_rows_in_deseq2_order() {
    let normalized = RowMajorMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 10.0, 12.0, 13.0, 6.0, 7.0, 40.0, 41.0, 8.0, 9.0, 100.0, 101.0, 9.0, 10.0, 7.0,
            8.0,
        ],
    )
    .unwrap();
    let base_mean = vec![1.0, 12.0, 6.0, 40.0, 8.0, 100.0, 9.0, 7.0];

    let subset = fast_vst_subset_normalized_counts(&normalized, &base_mean, 4).unwrap();

    assert_eq!(subset.n_rows(), 4);
    assert_eq!(subset.n_cols(), 2);
    assert_eq!(
        subset.as_slice(),
        &[6.0, 7.0, 8.0, 9.0, 12.0, 13.0, 100.0, 101.0]
    );
}

#[test]
fn fast_vst_subset_normalized_counts_validates_base_mean_length() {
    let normalized = RowMajorMatrix::from_row_major(2, 2, vec![6.0, 7.0, 8.0, 9.0]).unwrap();

    assert!(fast_vst_subset_normalized_counts(&normalized, &[6.5], 1).is_err());
}

#[test]
fn fast_vst_subset_matrix_rows_keeps_aligned_factor_rows() {
    let normalization_factors = RowMajorMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5,
        ],
    )
    .unwrap();
    let base_mean = vec![1.0, 12.0, 6.0, 40.0, 8.0, 100.0, 9.0, 7.0];

    let subset = fast_vst_subset_matrix_rows(
        &normalization_factors,
        &base_mean,
        4,
        "fast VST normalization factors",
    )
    .unwrap();

    assert_eq!(subset.n_rows(), 4);
    assert_eq!(subset.n_cols(), 2);
    assert_eq!(subset.as_slice(), &[1.4, 1.5, 1.8, 1.9, 1.2, 1.3, 2.0, 2.1]);
}

#[test]
fn count_matrix_select_rows_preserves_fast_vst_order_and_names() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        8,
        2,
        vec![1, 10, 12, 13, 6, 7, 40, 41, 8, 9, 100, 101, 9, 10, 7, 8],
        Some((0..8).map(|idx| format!("gene{idx}")).collect()),
        Some(vec!["sample0".to_string(), "sample1".to_string()]),
    )
    .unwrap();
    let base_mean = vec![1.0, 12.0, 6.0, 40.0, 8.0, 100.0, 9.0, 7.0];
    let row_indices = fast_vst_subset_indices(&base_mean, 4).unwrap();

    let subset = counts.select_rows(&row_indices).unwrap();

    assert_eq!(subset.n_genes(), 4);
    assert_eq!(subset.n_samples(), 2);
    assert_eq!(subset.as_slice(), &[6, 7, 8, 9, 12, 13, 100, 101]);
    assert_eq!(
        subset.gene_names().unwrap(),
        &[
            "gene2".to_string(),
            "gene4".to_string(),
            "gene1".to_string(),
            "gene5".to_string()
        ]
    );
    assert_eq!(
        subset.sample_names().unwrap(),
        &["sample0".to_string(), "sample1".to_string()]
    );
}

#[test]
fn fast_vst_subset_returns_aligned_counts_normalized_factors_and_weights() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        8,
        2,
        vec![1, 10, 12, 13, 6, 7, 40, 41, 8, 9, 100, 101, 9, 10, 7, 8],
        Some((0..8).map(|idx| format!("gene{idx}")).collect()),
        Some(vec!["sample0".to_string(), "sample1".to_string()]),
    )
    .unwrap();
    let normalized = RowMajorMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 10.0, 12.0, 13.0, 6.0, 7.0, 40.0, 41.0, 8.0, 9.0, 100.0, 101.0, 9.0, 10.0, 7.0,
            8.0,
        ],
    )
    .unwrap();
    let normalization_factors = RowMajorMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5,
        ],
    )
    .unwrap();
    let observation_weights = RowMajorMatrix::from_row_major(
        8,
        2,
        vec![
            0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
        ],
    )
    .unwrap();
    let base_mean = vec![1.0, 12.0, 6.0, 40.0, 8.0, 100.0, 9.0, 7.0];

    let subset = fast_vst_subset(
        &counts,
        &normalized,
        &base_mean,
        4,
        Some(&normalization_factors),
        Some(&observation_weights),
    )
    .unwrap();

    assert_eq!(subset.row_indices, vec![2, 4, 1, 5]);
    let metadata = subset.metadata();
    assert_eq!(metadata.rows, 4);
    assert_eq!(metadata.cols, 2);
    assert_eq!(metadata.row_indices, vec![2, 4, 1, 5]);
    assert!(metadata.has_normalization_factors);
    assert!(metadata.has_observation_weights);
    assert_eq!(subset.counts.as_slice(), &[6, 7, 8, 9, 12, 13, 100, 101]);
    assert_eq!(
        subset.normalized_counts.as_slice(),
        &[6.0, 7.0, 8.0, 9.0, 12.0, 13.0, 100.0, 101.0]
    );
    assert_eq!(
        subset.normalization_factors.unwrap().as_slice(),
        &[1.4, 1.5, 1.8, 1.9, 1.2, 1.3, 2.0, 2.1]
    );
    assert_eq!(
        subset.observation_weights.unwrap().as_slice(),
        &[0.5, 0.6, 0.9, 1.0, 0.3, 0.4, 1.1, 1.2]
    );
    assert_eq!(
        subset.counts.gene_names().unwrap(),
        &[
            "gene2".to_string(),
            "gene4".to_string(),
            "gene1".to_string(),
            "gene5".to_string()
        ]
    );
}

#[test]
fn fast_vst_subset_validates_all_aligned_matrix_shapes() {
    let counts = CountMatrix::from_row_major_u32(2, 2, vec![6, 7, 8, 9]).unwrap();
    let normalized = RowMajorMatrix::from_row_major(2, 2, vec![6.0, 7.0, 8.0, 9.0]).unwrap();
    let bad_matrix = RowMajorMatrix::from_row_major(1, 2, vec![1.0, 1.1]).unwrap();
    let base_mean = vec![6.5, 8.5];

    assert!(fast_vst_subset(&counts, &bad_matrix, &base_mean, 1, None, None).is_err());
    assert!(fast_vst_subset(&counts, &normalized, &base_mean, 1, Some(&bad_matrix), None).is_err());
    assert!(fast_vst_subset(&counts, &normalized, &base_mean, 1, None, Some(&bad_matrix)).is_err());
}

#[test]
fn vst_mean_applies_deseq2_mean_fit_closed_form() {
    let normalized =
        RowMajorMatrix::from_row_major(2, 3, vec![0.0, 1.0, 4.0, 10.0, 100.0, 1000.0]).unwrap();
    let alpha = 0.25;

    let transformed = vst_mean(&normalized, alpha).unwrap();

    assert_eq!(transformed.n_rows(), 2);
    assert_eq!(transformed.n_cols(), 3);
    for (observed, q) in transformed.as_slice().iter().zip(normalized.as_slice()) {
        assert_relative_eq!(*observed, expected_mean_vst(*q, alpha), epsilon = 1e-12);
    }
}

#[test]
fn public_vst_aliases_mean_fit_branch() {
    let normalized = RowMajorMatrix::from_row_major(1, 3, vec![2.0, 8.0, 32.0]).unwrap();

    let direct = vst_mean(&normalized, 0.5).unwrap();
    let aliased = vst(&normalized, 0.5).unwrap();

    assert_eq!(aliased, direct);
}

#[test]
fn vst_parametric_applies_deseq2_parametric_closed_form() {
    let normalized =
        RowMajorMatrix::from_row_major(2, 3, vec![0.0, 1.0, 4.0, 10.0, 100.0, 1000.0]).unwrap();
    let trend = ParametricDispersionTrend {
        asympt_disp: 0.2,
        extra_pois: 1.5,
    };

    let transformed = vst_parametric(&normalized, trend).unwrap();

    assert_eq!(transformed.n_rows(), 2);
    assert_eq!(transformed.n_cols(), 3);
    for (observed, q) in transformed.as_slice().iter().zip(normalized.as_slice()) {
        assert_relative_eq!(
            *observed,
            expected_parametric_vst(*q, trend),
            epsilon = 1e-12
        );
    }
}

#[test]
fn vst_local_applies_numerical_integration_and_log2_scaling() {
    let normalized = RowMajorMatrix::from_row_major(
        4,
        3,
        vec![
            1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0,
        ],
    )
    .unwrap();
    let trend = constant_local_trend(0.25);

    let transformed = vst_local(&normalized, &trend, 1.0).unwrap();

    assert_eq!(transformed.n_rows(), 4);
    assert_eq!(transformed.n_cols(), 3);
    for row in 0..transformed.n_rows() {
        let values = transformed.row(row).unwrap();
        assert!(values.windows(2).all(|pair| pair[1] > pair[0]));
    }
    let high = transformed.as_slice().last().copied().unwrap();
    assert_relative_eq!(high, 2048_f64.log2(), epsilon = 0.01);
}

#[test]
fn local_vst_size_factor_summary_matches_mean_inverse_size_factor() {
    let observed = local_vst_inverse_size_factor_mean(&[1.0, 2.0, 4.0]).unwrap();

    assert_relative_eq!(observed, (1.0 + 0.5 + 0.25) / 3.0, epsilon = 1e-12);
}

#[test]
fn local_vst_normalization_factor_summary_uses_column_geometric_means() {
    let normalization_factors =
        RowMajorMatrix::from_row_major(2, 3, vec![1.0, 2.0, 4.0, 4.0, 8.0, 16.0]).unwrap();

    let observed =
        local_vst_inverse_size_factor_mean_from_normalization_factors(&normalization_factors)
            .unwrap();

    let sf0 = (1.0_f64 * 4.0).sqrt();
    let sf1 = (2.0_f64 * 8.0).sqrt();
    let sf2 = (4.0_f64 * 16.0).sqrt();
    assert_relative_eq!(
        observed,
        (sf0.recip() + sf1.recip() + sf2.recip()) / 3.0,
        epsilon = 1e-12
    );
}

#[test]
fn local_vst_size_factor_summaries_reject_invalid_values() {
    assert!(local_vst_inverse_size_factor_mean(&[]).is_err());
    assert!(local_vst_inverse_size_factor_mean(&[1.0, 0.0]).is_err());

    let bad_factors = RowMajorMatrix::from_row_major(1, 2, vec![1.0, f64::NAN]).unwrap();
    assert!(local_vst_inverse_size_factor_mean_from_normalization_factors(&bad_factors).is_err());
}

#[test]
fn local_vst_size_factor_summaries_keep_large_finite_means() {
    let observed = local_vst_inverse_size_factor_mean(&[
        f64::MIN_POSITIVE,
        f64::MIN_POSITIVE,
        f64::MIN_POSITIVE,
        f64::MIN_POSITIVE,
    ])
    .unwrap();
    assert_relative_eq!(observed, f64::MIN_POSITIVE.recip(), max_relative = 1e-12);

    let factors = RowMajorMatrix::from_row_major(
        1,
        5,
        vec![
            f64::MIN_POSITIVE,
            f64::MIN_POSITIVE,
            f64::MIN_POSITIVE,
            f64::MIN_POSITIVE,
            f64::MIN_POSITIVE,
        ],
    )
    .unwrap();
    let observed = local_vst_inverse_size_factor_mean_from_normalization_factors(&factors).unwrap();
    assert_relative_eq!(observed, f64::MIN_POSITIVE.recip(), max_relative = 1e-12);
}

#[test]
fn vst_with_dispersion_trend_dispatches_parametric_mean_and_local_branches() {
    let normalized = RowMajorMatrix::from_row_major(
        4,
        3,
        vec![
            1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0,
        ],
    )
    .unwrap();

    let parametric_trend = ParametricDispersionTrend {
        asympt_disp: 0.2,
        extra_pois: 1.5,
    };
    let parametric = DispersionTrendFit::Parametric(ParametricDispersionTrendFit {
        trend: parametric_trend,
        disp_fit: vec![0.0; normalized.n_rows()],
        use_for_fit: vec![true; normalized.n_rows()],
        genes_used: normalized.n_rows(),
        outer_iterations: 1,
        irls_iterations: 1,
        converged: true,
    });
    assert_eq!(
        vst_with_dispersion_trend(&normalized, &parametric, 1.0).unwrap(),
        vst_parametric(&normalized, parametric_trend).unwrap()
    );

    let mean = DispersionTrendFit::Mean(MeanDispersionTrendFit {
        trend: MeanDispersionTrend { mean_disp: 0.25 },
        disp_fit: vec![0.25; normalized.n_rows()],
        use_for_fit: vec![true; normalized.n_rows()],
        use_for_mean: vec![true; normalized.n_rows()],
        genes_used: normalized.n_rows(),
        genes_used_for_mean: normalized.n_rows(),
    });
    assert_eq!(
        vst_with_dispersion_trend(&normalized, &mean, 1.0).unwrap(),
        vst_mean(&normalized, 0.25).unwrap()
    );

    let local_trend = constant_local_trend(0.25);
    let local = DispersionTrendFit::Local(LocalDispersionTrendFit {
        trend: local_trend.clone(),
        disp_fit: vec![0.25; normalized.n_rows()],
        use_for_fit: vec![true; normalized.n_rows()],
        genes_used: normalized.n_rows(),
        used_min_disp_floor: true,
    });
    assert_eq!(
        vst_with_dispersion_trend(&normalized, &local, 1.0).unwrap(),
        vst_local(&normalized, &local_trend, 1.0).unwrap()
    );
}

#[test]
fn vst_with_dispersion_trend_and_size_factors_computes_local_variance_term() {
    let normalized = RowMajorMatrix::from_row_major(
        4,
        3,
        vec![
            1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0,
        ],
    )
    .unwrap();
    let local_trend = constant_local_trend(0.25);
    let local = DispersionTrendFit::Local(LocalDispersionTrendFit {
        trend: local_trend.clone(),
        disp_fit: vec![0.25; normalized.n_rows()],
        use_for_fit: vec![true; normalized.n_rows()],
        genes_used: normalized.n_rows(),
        used_min_disp_floor: true,
    });
    let size_factors = [1.0, 2.0, 4.0];
    let inverse_size_factor_mean = local_vst_inverse_size_factor_mean(&size_factors).unwrap();

    assert_eq!(
        vst_with_dispersion_trend_and_size_factors(&normalized, &local, &size_factors).unwrap(),
        vst_local(&normalized, &local_trend, inverse_size_factor_mean).unwrap()
    );
}

#[test]
fn vst_with_dispersion_trend_and_normalization_factors_computes_local_variance_term() {
    let normalized = RowMajorMatrix::from_row_major(
        4,
        3,
        vec![
            1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0,
        ],
    )
    .unwrap();
    let normalization_factors = RowMajorMatrix::from_row_major(
        4,
        3,
        vec![
            1.0, 2.0, 4.0, 4.0, 8.0, 16.0, 2.0, 4.0, 8.0, 8.0, 16.0, 32.0,
        ],
    )
    .unwrap();
    let local_trend = constant_local_trend(0.25);
    let local = DispersionTrendFit::Local(LocalDispersionTrendFit {
        trend: local_trend.clone(),
        disp_fit: vec![0.25; normalized.n_rows()],
        use_for_fit: vec![true; normalized.n_rows()],
        genes_used: normalized.n_rows(),
        used_min_disp_floor: true,
    });
    let inverse_size_factor_mean =
        local_vst_inverse_size_factor_mean_from_normalization_factors(&normalization_factors)
            .unwrap();

    assert_eq!(
        vst_with_dispersion_trend_and_normalization_factors(
            &normalized,
            &local,
            &normalization_factors,
        )
        .unwrap(),
        vst_local(&normalized, &local_trend, inverse_size_factor_mean).unwrap()
    );
}

#[test]
fn vst_with_dispersion_trend_factor_helpers_reject_bad_factor_inputs() {
    let normalized =
        RowMajorMatrix::from_row_major(2, 3, vec![1.0, 2.0, 4.0, 8.0, 16.0, 32.0]).unwrap();
    let mean = DispersionTrendFit::Mean(MeanDispersionTrendFit {
        trend: MeanDispersionTrend { mean_disp: 0.25 },
        disp_fit: vec![0.25; normalized.n_rows()],
        use_for_fit: vec![true; normalized.n_rows()],
        use_for_mean: vec![true; normalized.n_rows()],
        genes_used: normalized.n_rows(),
        genes_used_for_mean: normalized.n_rows(),
    });

    assert!(vst_with_dispersion_trend_and_size_factors(&normalized, &mean, &[1.0, 2.0]).is_err());
    assert!(vst_with_dispersion_trend_and_size_factors(&normalized, &mean, &[1.0, 0.0]).is_err());

    let wrong_shape = RowMajorMatrix::from_row_major(1, 3, vec![1.0, 2.0, 4.0]).unwrap();
    assert!(
        vst_with_dispersion_trend_and_normalization_factors(&normalized, &mean, &wrong_shape)
            .is_err()
    );
}

#[test]
fn vst_mean_value_is_log2_like_for_large_counts() {
    let q = 1e10_f64;
    let alpha = 0.25_f64;

    let transformed = vst_mean_value(q, alpha, 0).unwrap();

    assert_relative_eq!(transformed, q.log2(), epsilon = 1e-9);
}

#[test]
fn vst_mean_value_stays_finite_when_dispersion_count_overflows() {
    let q = f64::MAX;
    let transformed = vst_mean_value(q, 2.0, 0).unwrap();

    assert!(transformed.is_finite());
    assert_relative_eq!(transformed, q.log2(), max_relative = 1e-12);
}

#[test]
fn vst_parametric_value_is_log2_like_for_large_counts() {
    let q = 1e10_f64;
    let trend = ParametricDispersionTrend {
        asympt_disp: 0.2,
        extra_pois: 1.5,
    };

    let transformed = vst_parametric_value(q, trend, 0).unwrap();

    assert_relative_eq!(transformed, q.log2(), epsilon = 1e-9);
}

#[test]
fn vst_parametric_value_stays_finite_for_extreme_counts() {
    let q = 1e308_f64;
    let trend = ParametricDispersionTrend {
        asympt_disp: 0.2,
        extra_pois: 1.5,
    };

    let transformed = vst_parametric_value(q, trend, 0).unwrap();

    assert!(transformed.is_finite());
    assert_relative_eq!(transformed, q.log2(), max_relative = 1e-12);
}

#[test]
fn vst_parametric_value_stays_finite_when_dispersion_count_overflows() {
    let q = f64::MAX;
    let trend = ParametricDispersionTrend {
        asympt_disp: 2.0,
        extra_pois: 1.5,
    };

    let transformed = vst_parametric_value(q, trend, 0).unwrap();

    assert!(transformed.is_finite());
    assert_relative_eq!(transformed, q.log2(), max_relative = 1e-12);
}

#[test]
fn vst_mean_rejects_negative_and_non_finite_counts() {
    let negative = RowMajorMatrix::from_row_major(1, 1, vec![-1.0]).unwrap();
    assert!(vst_mean(&negative, 0.25).is_err());

    let non_finite = RowMajorMatrix::from_row_major(1, 1, vec![f64::INFINITY]).unwrap();
    assert!(vst_mean(&non_finite, 0.25).is_err());
}

#[test]
fn vst_parametric_rejects_negative_and_non_finite_counts() {
    let trend = ParametricDispersionTrend {
        asympt_disp: 0.2,
        extra_pois: 1.5,
    };

    let negative = RowMajorMatrix::from_row_major(1, 1, vec![-1.0]).unwrap();
    assert!(vst_parametric(&negative, trend).is_err());

    let non_finite = RowMajorMatrix::from_row_major(1, 1, vec![f64::INFINITY]).unwrap();
    assert!(vst_parametric(&non_finite, trend).is_err());
}

#[test]
fn vst_local_rejects_negative_counts_and_bad_size_factor_summary() {
    let trend = constant_local_trend(0.25);
    let negative = RowMajorMatrix::from_row_major(1, 1, vec![-1.0]).unwrap();
    assert!(vst_local(&negative, &trend, 1.0).is_err());

    let normalized = RowMajorMatrix::from_row_major(2, 2, vec![1.0, 2.0, 4.0, 8.0]).unwrap();
    assert!(vst_local(&normalized, &trend, 0.0).is_err());
    assert!(vst_local(&normalized, &trend, f64::NAN).is_err());
}

#[test]
fn vst_local_rejects_overflowed_variance_curve() {
    let trend = constant_local_trend(0.25);
    let huge = RowMajorMatrix::from_row_major(1, 1, vec![f64::MAX]).unwrap();
    let err = vst_local(&huge, &trend, 1.0).unwrap_err();
    assert!(matches!(
        err,
        DeseqError::NonFiniteValue { .. } | DeseqError::InvalidDispersion { .. }
    ));
}

#[test]
fn vst_mean_rejects_non_positive_and_non_finite_dispersion() {
    let normalized = RowMajorMatrix::from_row_major(1, 1, vec![1.0]).unwrap();

    assert!(vst_mean(&normalized, 0.0).is_err());
    assert!(vst_mean(&normalized, -0.1).is_err());
    assert!(vst_mean(&normalized, f64::NAN).is_err());
}

#[test]
fn vst_parametric_rejects_invalid_trend_coefficients() {
    let normalized = RowMajorMatrix::from_row_major(1, 1, vec![1.0]).unwrap();

    assert!(vst_parametric(
        &normalized,
        ParametricDispersionTrend {
            asympt_disp: 0.0,
            extra_pois: 1.0,
        },
    )
    .is_err());
    assert!(vst_parametric(
        &normalized,
        ParametricDispersionTrend {
            asympt_disp: 0.2,
            extra_pois: -0.1,
        },
    )
    .is_err());
}

fn constant_local_trend(dispersion: f64) -> LocalDispersionTrend {
    LocalDispersionTrend {
        min_disp: dispersion,
        span: 0.7,
        degree: 2,
        log_means: Vec::new(),
        log_disps: Vec::new(),
        weights: Vec::new(),
    }
}
