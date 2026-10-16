use super::*;

#[test]
fn fitted_state_builds_fast_vst_subset_from_base_mean_and_aligned_inputs() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let weights = RowMajorMatrix::from_row_major(
        counts.n_genes(),
        counts.n_samples(),
        (0..counts.n_genes() * counts.n_samples())
            .map(|idx| 0.5 + idx as f64 / 100.0)
            .collect(),
    )
    .unwrap();
    let fit = glm_mu_native_wald_builder()
        .observation_weights(weights)
        .fit_size_factors_and_base_means_with_design(&counts, &design)
        .unwrap();

    let subset = fit.fast_vst_subset(&counts, 2).unwrap();

    assert_eq!(fit.fast_vst_eligible_count().unwrap(), 7);
    assert_eq!(subset.row_indices.len(), 2);
    assert_eq!(subset.counts.n_genes(), 2);
    assert_eq!(subset.counts.n_samples(), counts.n_samples());
    assert_eq!(subset.normalized_counts.n_rows(), 2);
    assert_eq!(subset.normalized_counts.n_cols(), counts.n_samples());
    assert_eq!(
        subset.observation_weights.as_ref().unwrap().n_rows(),
        subset.row_indices.len()
    );
    for (subset_row, original_row) in subset.row_indices.iter().copied().enumerate() {
        assert_eq!(
            subset.counts.row(subset_row).unwrap(),
            counts.row(original_row).unwrap()
        );
        assert_eq!(
            subset
                .observation_weights
                .as_ref()
                .unwrap()
                .row(subset_row)
                .unwrap(),
            fit.observation_weights
                .as_ref()
                .unwrap()
                .row(original_row)
                .unwrap()
        );
    }
}

#[test]
fn builder_fits_glm_mu_dispersion_trend_on_fast_vst_subset() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);

    let (subset_fit, subset) = builder
        .fit_fast_vst_dispersion_trend_glm_mu(&counts, &design, 4)
        .unwrap();

    assert_eq!(subset.row_indices, vec![1, 7, 4, 6]);
    assert_eq!(subset.counts.n_genes(), 4);
    assert_eq!(subset.counts.n_samples(), counts.n_samples());
    assert_eq!(subset_fit.counts_summary.n_genes, 4);
    assert_eq!(subset_fit.counts_summary.n_samples, counts.n_samples());
    assert_eq!(subset_fit.size_factors, vec![1.0; counts.n_samples()]);
    assert!(matches!(
        subset_fit.dispersion_trend.as_ref(),
        Some(DispersionTrendFit::Mean(_))
    ));
    assert!(
        subset_fit
            .disp_fit
            .as_ref()
            .unwrap()
            .iter()
            .all(|value| value.is_finite())
    );
}

#[test]
fn builder_fast_vst_glm_mu_applies_subset_trend_to_full_counts() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);

    let output = builder.fast_vst_glm_mu(&counts, &design, 4).unwrap();
    let normalized = normalized_counts(&counts, &output.subset_fit.size_factors).unwrap();
    let expected = vst_with_dispersion_trend_and_size_factors(
        &normalized,
        output.subset_fit.dispersion_trend.as_ref().unwrap(),
        &output.subset_fit.size_factors,
    )
    .unwrap();

    assert_eq!(output.subset.row_indices, vec![1, 7, 4, 6]);
    assert_eq!(output.transformed.n_rows(), counts.n_genes());
    assert_eq!(output.transformed.n_cols(), counts.n_samples());
    let metadata = output.metadata();
    assert_eq!(metadata.transformed_rows, counts.n_genes());
    assert_eq!(metadata.transformed_cols, counts.n_samples());
    assert_eq!(metadata.fast_subset_rows, 4);
    assert_eq!(metadata.fast_subset_cols, counts.n_samples());
    assert_eq!(metadata.fast_subset_indices, vec![1, 7, 4, 6]);
    assert_eq!(metadata.trend_fit_rows, 4);
    assert_eq!(metadata.trend_fit_cols, counts.n_samples());
    assert_eq!(metadata.trend_fit_type, Some("mean"));
    assert_eq!(output.transformed, expected);
    assert!(
        output
            .transformed
            .as_slice()
            .iter()
            .all(|value| value.is_finite())
    );
}

#[test]
fn builder_default_fast_vst_uses_deseq2_nsub_default() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);

    let err = builder
        .default_fast_vst_glm_mu(&counts, &design)
        .unwrap_err();

    assert_eq!(DEFAULT_FAST_VST_NSUB, 1000);
    assert!(err.to_string().contains("1000"));
}

#[test]
fn builder_auto_vst_uses_full_trend_when_default_fast_subset_is_too_large() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);

    let output = builder.default_vst_glm_mu_auto(&counts, &design).unwrap();
    let full_fit = builder
        .fit_dispersion_trend_glm_mu(&counts, &design)
        .unwrap();
    let expected = full_fit.variance_stabilizing_transform(&counts).unwrap();

    assert!(output.fast_subset.is_none());
    assert_eq!(
        output.trend_source,
        VstTrendSource::FullData {
            nsub: DEFAULT_FAST_VST_NSUB,
            eligible_rows: 7,
            reason: VstFullDataReason::InsufficientEligibleRows,
        }
    );
    assert_eq!(output.trend_source.nsub(), DEFAULT_FAST_VST_NSUB);
    assert_eq!(output.trend_source.eligible_rows(), 7);
    assert!(!output.trend_source.used_fast_subset());
    assert_eq!(output.trend_source.as_str(), "fullData");
    assert_eq!(
        output.trend_source.full_data_reason(),
        Some(VstFullDataReason::InsufficientEligibleRows)
    );
    assert_eq!(
        output.trend_source.full_data_reason().unwrap().as_str(),
        "insufficientEligibleRows"
    );
    let metadata = output.metadata();
    assert_eq!(metadata.trend_source, "fullData");
    assert_eq!(metadata.nsub, DEFAULT_FAST_VST_NSUB);
    assert_eq!(metadata.eligible_rows, 7);
    assert!(!metadata.used_fast_subset);
    assert_eq!(metadata.full_data_reason, Some("insufficientEligibleRows"));
    assert_eq!(metadata.transformed_rows, counts.n_genes());
    assert_eq!(metadata.transformed_cols, counts.n_samples());
    assert_eq!(metadata.trend_fit_rows, counts.n_genes());
    assert_eq!(metadata.trend_fit_cols, counts.n_samples());
    assert_eq!(metadata.trend_fit_type, Some("mean"));
    assert_eq!(metadata.fast_subset_rows, None);
    assert_eq!(metadata.fast_subset_indices, None);
    assert_eq!(output.trend_fit.counts_summary.n_genes, counts.n_genes());
    assert_eq!(output.transformed, expected);
}

#[test]
fn builder_auto_vst_uses_fast_subset_when_enough_rows_are_eligible() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);

    assert!(builder.vst_glm_mu_auto(&counts, &design, 0).is_err());

    let output = builder.vst_glm_mu_auto(&counts, &design, 4).unwrap();
    let fast = builder.fast_vst_glm_mu(&counts, &design, 4).unwrap();

    assert_eq!(
        output.fast_subset.as_ref().unwrap().row_indices,
        vec![1, 7, 4, 6]
    );
    assert_eq!(
        output.trend_source,
        VstTrendSource::FastSubset {
            nsub: 4,
            eligible_rows: 7,
        }
    );
    assert_eq!(output.trend_source.nsub(), 4);
    assert_eq!(output.trend_source.eligible_rows(), 7);
    assert!(output.trend_source.used_fast_subset());
    assert_eq!(output.trend_source.as_str(), "fastSubset");
    assert_eq!(output.trend_source.full_data_reason(), None);
    let metadata = output.metadata();
    assert_eq!(metadata.trend_source, "fastSubset");
    assert_eq!(metadata.nsub, 4);
    assert_eq!(metadata.eligible_rows, 7);
    assert!(metadata.used_fast_subset);
    assert_eq!(metadata.full_data_reason, None);
    assert_eq!(metadata.transformed_rows, counts.n_genes());
    assert_eq!(metadata.transformed_cols, counts.n_samples());
    assert_eq!(metadata.trend_fit_rows, 4);
    assert_eq!(metadata.trend_fit_cols, counts.n_samples());
    assert_eq!(metadata.trend_fit_type, Some("mean"));
    assert_eq!(metadata.fast_subset_rows, Some(4));
    assert_eq!(metadata.fast_subset_indices, Some(vec![1, 7, 4, 6]));
    assert_eq!(output.transformed, fast.transformed);
    assert_eq!(output.trend_fit, fast.subset_fit);
}

#[test]
fn builder_auto_vst_uses_fast_subset_when_observation_weights_are_present() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let weights = nonunit_weights_for(&counts);
    let builder = glm_mu_native_wald_builder()
        .fit_type(FitType::Mean)
        .observation_weights(weights.clone());

    let output = builder.vst_glm_mu_auto(&counts, &design, 4).unwrap();
    let fast = builder.fast_vst_glm_mu(&counts, &design, 4).unwrap();
    let base_fit = builder
        .fit_size_factors_and_base_means_with_design(&counts, &design)
        .unwrap();

    assert_eq!(
        output.fast_subset.as_ref().unwrap().row_indices,
        vec![1, 7, 4, 6]
    );
    assert_eq!(
        output.trend_source,
        VstTrendSource::FastSubset {
            nsub: 4,
            eligible_rows: 7,
        }
    );
    assert!(output.trend_source.used_fast_subset());
    assert_eq!(output.trend_source.as_str(), "fastSubset");
    assert_eq!(output.trend_source.full_data_reason(), None);
    let metadata = output.metadata();
    assert_eq!(metadata.trend_source, "fastSubset");
    assert_eq!(metadata.nsub, 4);
    assert_eq!(metadata.eligible_rows, 7);
    assert!(metadata.used_fast_subset);
    assert_eq!(metadata.full_data_reason, None);
    assert_eq!(metadata.transformed_rows, counts.n_genes());
    assert_eq!(metadata.transformed_cols, counts.n_samples());
    assert_eq!(metadata.trend_fit_rows, 4);
    assert_eq!(metadata.trend_fit_cols, counts.n_samples());
    assert_eq!(metadata.trend_fit_type, Some("mean"));
    assert_eq!(metadata.fast_subset_rows, Some(4));
    assert_eq!(metadata.fast_subset_indices, Some(vec![1, 7, 4, 6]));
    assert_eq!(output.transformed, fast.transformed);
    assert_eq!(output.trend_fit, fast.subset_fit);
    assert!(output.trend_fit.dispersion_trend.is_some());
    assert!(output.trend_fit.observation_weights.is_some());
    assert!(
        output
            .fast_subset
            .as_ref()
            .unwrap()
            .observation_weights
            .is_some()
    );
    for (subset_row, original_row) in output
        .fast_subset
        .as_ref()
        .unwrap()
        .row_indices
        .iter()
        .copied()
        .enumerate()
    {
        assert_eq!(
            output
                .fast_subset
                .as_ref()
                .unwrap()
                .observation_weights
                .as_ref()
                .unwrap()
                .row(subset_row)
                .unwrap(),
            base_fit
                .observation_weights
                .as_ref()
                .unwrap()
                .row(original_row)
                .unwrap()
        );
    }
}

#[test]
fn builder_blind_auto_vst_uses_intercept_only_design() {
    let counts = native_wald_counts_with_zero_row();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);
    let blind_design = DesignMatrix::intercept_only(counts.n_samples()).unwrap();

    let output = builder.blind_vst_glm_mu_auto(&counts, 4).unwrap();
    let expected = builder.vst_glm_mu_auto(&counts, &blind_design, 4).unwrap();

    assert_eq!(
        output
            .trend_fit
            .design
            .as_ref()
            .unwrap()
            .coefficient_names()
            .unwrap(),
        &["Intercept".to_string()]
    );
    assert_eq!(
        output.fast_subset.as_ref().unwrap().row_indices,
        vec![1, 7, 4, 6]
    );
    assert_eq!(output.metadata(), expected.metadata());
    assert_eq!(output.transformed, expected.transformed);
    assert_eq!(output.trend_fit, expected.trend_fit);
}

#[test]
fn fitted_state_can_apply_external_subset_trend_to_full_counts() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);
    let full_fit = builder
        .fit_size_factors_and_base_means_with_design(&counts, &design)
        .unwrap();
    let (subset_fit, _subset) = builder
        .fit_fast_vst_dispersion_trend_glm_mu(&counts, &design, 4)
        .unwrap();

    let transformed = full_fit
        .variance_stabilizing_transform_with_trend(
            &counts,
            subset_fit.dispersion_trend.as_ref().unwrap(),
        )
        .unwrap();
    let normalized = full_fit.normalized_counts(&counts).unwrap();
    let expected = vst_with_dispersion_trend_and_size_factors(
        &normalized,
        subset_fit.dispersion_trend.as_ref().unwrap(),
        &full_fit.size_factors,
    )
    .unwrap();

    assert_eq!(transformed, expected);
    assert_eq!(transformed.n_rows(), counts.n_genes());
    assert_eq!(transformed.n_cols(), counts.n_samples());
}

#[test]
fn builder_fast_vst_glm_mu_preserves_normalization_factors() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let normalization_factors =
        RowMajorMatrix::from_elem(counts.n_genes(), counts.n_samples(), 1.1).unwrap();
    let builder = glm_mu_native_wald_builder()
        .fit_type(FitType::Mean)
        .normalization_factors(normalization_factors.clone());

    let output = builder.fast_vst_glm_mu(&counts, &design, 4).unwrap();
    let normalized = normalized_counts_with_factors(&counts, &normalization_factors).unwrap();
    let expected = vst_with_dispersion_trend_and_normalization_factors(
        &normalized,
        output.subset_fit.dispersion_trend.as_ref().unwrap(),
        &normalization_factors,
    )
    .unwrap();

    assert_eq!(output.subset.row_indices, vec![1, 7, 4, 6]);
    assert_eq!(
        output
            .subset
            .normalization_factors
            .as_ref()
            .unwrap()
            .n_rows(),
        output.subset.row_indices.len()
    );
    for (subset_row, original_row) in output.subset.row_indices.iter().copied().enumerate() {
        assert_eq!(
            output
                .subset
                .normalization_factors
                .as_ref()
                .unwrap()
                .row(subset_row)
                .unwrap(),
            normalization_factors.row(original_row).unwrap()
        );
    }
    assert_eq!(
        output.subset_fit.normalization_factors,
        output.subset.normalization_factors.clone()
    );
    assert_eq!(output.transformed, expected);
}

#[test]
fn builder_fast_vst_subset_trend_preserves_observation_weights() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let n_samples = counts.n_samples();
    let weights = RowMajorMatrix::from_row_major(
        counts.n_genes(),
        n_samples,
        (0..counts.n_genes())
            .flat_map(|gene| {
                (0..n_samples).map(move |sample| if sample == gene % n_samples { 1.0 } else { 0.2 })
            })
            .collect(),
    )
    .unwrap();
    let builder = glm_mu_native_wald_builder().observation_weights(weights.clone());
    let base_fit = builder
        .fit_size_factors_and_base_means_with_design(&counts, &design)
        .unwrap();

    let invalid_nsub = builder
        .fit_fast_vst_dispersion_trend_glm_mu(&counts, &design, 0)
        .unwrap_err();
    assert!(matches!(invalid_nsub, DeseqError::InvalidOptions { .. }));
    assert!(
        invalid_nsub
            .to_string()
            .contains("fast VST subset size must be positive")
    );

    let (subset_fit, subset) = builder
        .fit_fast_vst_dispersion_trend_glm_mu(&counts, &design, 4)
        .unwrap();

    let normalized = base_fit.normalized_counts(&counts).unwrap();
    let ordinary_means = base_mean(&normalized).unwrap();
    let expected_rows = fast_vst_subset_indices(&ordinary_means, 4).unwrap();
    assert_eq!(subset.row_indices, expected_rows);
    assert_ne!(
        expected_rows,
        fast_vst_subset_indices(&base_fit.base_mean, 4).unwrap(),
        "fixture must distinguish ordinary from weighted selection means"
    );

    assert!(subset.observation_weights.is_some());
    assert!(subset_fit.observation_weights.is_some());
    assert_eq!(
        subset_fit.observation_weights,
        subset.observation_weights.clone()
    );
    assert_eq!(subset_fit.counts_summary.n_genes, subset.row_indices.len());
    for (subset_row, original_row) in subset.row_indices.iter().copied().enumerate() {
        assert_eq!(
            subset
                .observation_weights
                .as_ref()
                .unwrap()
                .row(subset_row)
                .unwrap(),
            base_fit
                .observation_weights
                .as_ref()
                .unwrap()
                .row(original_row)
                .unwrap()
        );
    }
}

#[test]
fn fitted_vst_requires_stored_dispersion_trend_and_matching_counts() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let fit = glm_mu_native_wald_builder()
        .fit_size_factors_and_base_means_with_design(&counts, &design)
        .unwrap();

    assert!(fit.variance_stabilizing_transform(&counts).is_err());

    let wrong_counts = CountMatrix::from_row_major_u32(1, 8, vec![1; 8]).unwrap();
    assert!(fit.normalized_counts(&wrong_counts).is_err());
}
