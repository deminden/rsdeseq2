use super::*;

#[test]
fn fitted_trend_state_drives_norm_and_vst_transforms() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let fit = glm_mu_native_wald_builder()
        .fit_type(FitType::Local)
        .fit_dispersion_trend_glm_mu(&counts, &design)
        .unwrap();

    assert!(matches!(
        fit.dispersion_trend.as_ref(),
        Some(DispersionTrendFit::Local(_))
    ));

    let normalized = fit.normalized_counts(&counts).unwrap();
    let norm = fit.norm_transform(&counts).unwrap();
    let vst = fit.variance_stabilizing_transform(&counts).unwrap();

    assert_eq!(normalized.n_rows(), counts.n_genes());
    assert_eq!(normalized.n_cols(), counts.n_samples());
    assert_eq!(norm.n_rows(), counts.n_genes());
    assert_eq!(vst.n_rows(), counts.n_genes());
    assert!(vst.as_slice().iter().all(|value| value.is_finite()));
    assert_relative_eq!(norm.as_slice()[9], (20.0_f64 + 1.0).log2(), epsilon = 1e-12);
}

#[test]
fn fitted_trend_state_drives_rlog_transform() {
    let full_counts = native_wald_counts_with_zero_row();
    let counts = full_counts.select_rows(&[1, 2, 3, 4, 5, 6, 7]).unwrap();
    let design = two_group_design();
    let fit = glm_mu_native_wald_builder()
        .fit_type(FitType::Mean)
        .fit_map_dispersions_glm_mu(&counts, &design)
        .unwrap();

    let output = fit.rlog(&counts).unwrap();
    let expected = rlog_with_estimated_prior_and_size_factors(
        &counts,
        &fit.size_factors,
        &fit.base_mean,
        fit.disp_fit.as_ref().unwrap(),
        fit.disp_fit.as_ref().unwrap(),
        IrlsOptions::default(),
    )
    .unwrap();

    assert_eq!(output.transformed.n_rows(), counts.n_genes());
    assert_eq!(output.transformed.n_cols(), counts.n_samples());
    assert!(output.sample_prior_variance.is_finite());
    assert!(output.sample_prior_variance > 0.0);
    assert_eq!(output.offset_mode, RlogOffsetMode::SizeFactors);
    assert_eq!(output.metadata().transformed_rows, counts.n_genes());
    assert_eq!(output.metadata().transformed_cols, counts.n_samples());
    assert_eq!(output.metadata().intercept_len, counts.n_genes());
    assert_eq!(output.metadata().offset_mode, "sizeFactors");
    assert_relative_eq!(
        output.sample_prior_variance,
        expected.sample_prior_variance,
        epsilon = 1e-12
    );
    assert_eq!(
        output.transformed.as_slice(),
        expected.transformed.as_slice()
    );
    assert_eq!(output.intercept, expected.intercept);
}

#[test]
fn fitted_trend_state_drives_rlog_transform_with_normalization_factors() {
    let full_counts = native_wald_counts_with_zero_row();
    let counts = full_counts.select_rows(&[1, 2, 3, 4, 5, 6, 7]).unwrap();
    let design = two_group_design();
    let mut factor_values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for gene in 0..counts.n_genes() {
        for sample in 0..counts.n_samples() {
            factor_values.push(0.8 + 0.05 * gene as f64 + 0.02 * sample as f64);
        }
    }
    let normalization_factors =
        RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), factor_values)
            .unwrap();
    let fit = glm_mu_native_wald_builder()
        .fit_type(FitType::Mean)
        .normalization_factors(normalization_factors.clone())
        .fit_map_dispersions_glm_mu(&counts, &design)
        .unwrap();

    let output = fit.rlog(&counts).unwrap();
    let expected = rlog_with_estimated_prior_and_normalization_factors(
        &counts,
        fit.normalization_factors.as_ref().unwrap(),
        &fit.base_mean,
        fit.disp_fit.as_ref().unwrap(),
        fit.disp_fit.as_ref().unwrap(),
        IrlsOptions::default(),
    )
    .unwrap();

    assert_eq!(fit.normalization_factors, Some(normalization_factors));
    assert!(output.sample_prior_variance.is_finite());
    assert!(output.sample_prior_variance > 0.0);
    assert_eq!(output.offset_mode, RlogOffsetMode::NormalizationFactors);
    assert_eq!(output.metadata().transformed_rows, counts.n_genes());
    assert_eq!(output.metadata().transformed_cols, counts.n_samples());
    assert_eq!(output.metadata().intercept_len, counts.n_genes());
    assert_eq!(output.metadata().offset_mode, "normalizationFactors");
    assert_relative_eq!(
        output.sample_prior_variance,
        expected.sample_prior_variance,
        epsilon = 1e-12
    );
    assert_eq!(
        output.transformed.as_slice(),
        expected.transformed.as_slice()
    );
    assert_eq!(output.intercept, expected.intercept);
}

#[test]
fn fitted_trend_state_drives_frozen_rlog_transform() {
    let full_counts = native_wald_counts_with_zero_row();
    let counts = full_counts.select_rows(&[1, 2, 3, 4, 5, 6, 7]).unwrap();
    let design = two_group_design();
    let fit = glm_mu_native_wald_builder()
        .fit_type(FitType::Mean)
        .fit_map_dispersions_glm_mu(&counts, &design)
        .unwrap();
    let learned = fit.rlog(&counts).unwrap();

    let output = fit
        .frozen_rlog(&counts, &learned.intercept, learned.sample_prior_variance)
        .unwrap();
    let expected = rsdeseq2::prelude::rlog_frozen_with_size_factors(
        &counts,
        &fit.size_factors,
        fit.disp_fit.as_ref().unwrap(),
        learned.sample_prior_variance,
        &learned.intercept,
        IrlsOptions::default(),
    )
    .unwrap();

    assert_eq!(output.offset_mode, RlogOffsetMode::SizeFactors);
    assert_eq!(output.metadata().intercept_len, counts.n_genes());
    assert_eq!(output.intercept, learned.intercept);
    assert_eq!(
        output.transformed.as_slice(),
        expected.transformed.as_slice()
    );
}

#[test]
fn fitted_trend_state_drives_frozen_rlog_transform_with_normalization_factors() {
    let full_counts = native_wald_counts_with_zero_row();
    let counts = full_counts.select_rows(&[1, 2, 3, 4, 5, 6, 7]).unwrap();
    let design = two_group_design();
    let mut factor_values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for gene in 0..counts.n_genes() {
        for sample in 0..counts.n_samples() {
            factor_values.push(0.8 + 0.05 * gene as f64 + 0.02 * sample as f64);
        }
    }
    let normalization_factors =
        RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), factor_values)
            .unwrap();
    let fit = glm_mu_native_wald_builder()
        .fit_type(FitType::Mean)
        .normalization_factors(normalization_factors)
        .fit_map_dispersions_glm_mu(&counts, &design)
        .unwrap();
    let learned = fit.rlog(&counts).unwrap();

    let output = fit
        .frozen_rlog(&counts, &learned.intercept, learned.sample_prior_variance)
        .unwrap();
    let expected = rsdeseq2::prelude::rlog_frozen_with_normalization_factors(
        &counts,
        fit.normalization_factors.as_ref().unwrap(),
        fit.disp_fit.as_ref().unwrap(),
        learned.sample_prior_variance,
        &learned.intercept,
        IrlsOptions::default(),
    )
    .unwrap();

    assert_eq!(output.offset_mode, RlogOffsetMode::NormalizationFactors);
    assert_eq!(output.intercept, learned.intercept);
    assert_eq!(
        output.transformed.as_slice(),
        expected.transformed.as_slice()
    );
}

#[test]
fn fitted_rlog_uses_trend_dispersions_and_rejects_mismatched_counts() {
    let counts = native_wald_counts_with_zero_row()
        .select_rows(&[1, 2, 3, 4, 5, 6, 7])
        .unwrap();
    let design = two_group_design();
    let fit = glm_mu_native_wald_builder()
        .fit_type(FitType::Mean)
        .fit_dispersion_trend_glm_mu(&counts, &design)
        .unwrap();
    let output = fit.rlog(&counts).unwrap();
    assert_eq!(output.transformed.n_rows(), counts.n_genes());

    let wrong_counts = CountMatrix::from_row_major_u32(1, counts.n_samples(), vec![1; 8]).unwrap();
    assert!(fit.rlog(&wrong_counts).is_err());
}

#[test]
fn fitted_frozen_rlog_uses_trend_dispersions_and_validates_inputs() {
    let counts = native_wald_counts_with_zero_row()
        .select_rows(&[1, 2, 3, 4, 5, 6, 7])
        .unwrap();
    let design = two_group_design();
    let fit = glm_mu_native_wald_builder()
        .fit_type(FitType::Mean)
        .fit_dispersion_trend_glm_mu(&counts, &design)
        .unwrap();
    let frozen_intercept = vec![1.0; counts.n_genes()];
    let output = fit.frozen_rlog(&counts, &frozen_intercept, 1.0).unwrap();
    assert_eq!(output.transformed.n_rows(), counts.n_genes());

    let map_fit = glm_mu_native_wald_builder()
        .fit_type(FitType::Mean)
        .fit_map_dispersions_glm_mu(&counts, &design)
        .unwrap();
    assert!(map_fit.frozen_rlog(&counts, &[1.0; 2], 1.0).is_err());
    assert!(
        map_fit
            .frozen_rlog(&counts, &frozen_intercept, 0.0)
            .is_err()
    );
}

#[test]
fn fitted_rlog_expands_all_zero_rows_as_zero_transform_rows() {
    let counts = native_wald_counts_with_zero_row();
    let compact_counts = counts.select_rows(&[1, 2, 3, 4, 5, 6, 7]).unwrap();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);
    let fit = builder
        .fit_map_dispersions_glm_mu(&counts, &design)
        .unwrap();
    let compact_fit = builder
        .fit_map_dispersions_glm_mu(&compact_counts, &design)
        .unwrap();

    let output = fit.rlog(&counts).unwrap();
    let compact_output = compact_fit.rlog(&compact_counts).unwrap();

    assert_eq!(output.transformed.n_rows(), counts.n_genes());
    assert_eq!(output.transformed.n_cols(), counts.n_samples());
    assert_eq!(output.intercept.len(), counts.n_genes());
    assert_eq!(output.intercept[0], 0.0);
    assert!(
        output
            .transformed
            .row(0)
            .unwrap()
            .iter()
            .all(|value| *value == 0.0)
    );
    assert_eq!(
        output.sample_prior_variance,
        compact_output.sample_prior_variance
    );
    for compact_gene in 0..compact_counts.n_genes() {
        assert_eq!(
            output.transformed.row(compact_gene + 1).unwrap(),
            compact_output.transformed.row(compact_gene).unwrap()
        );
        assert_eq!(
            output.intercept[compact_gene + 1],
            compact_output.intercept[compact_gene]
        );
    }
}

#[test]
fn fitted_frozen_rlog_expands_all_zero_rows_with_frozen_intercepts() {
    let counts = native_wald_counts_with_zero_row();
    let compact_counts = counts.select_rows(&[1, 2, 3, 4, 5, 6, 7]).unwrap();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);
    let fit = builder
        .fit_map_dispersions_glm_mu(&counts, &design)
        .unwrap();
    let compact_fit = builder
        .fit_map_dispersions_glm_mu(&compact_counts, &design)
        .unwrap();
    let frozen_intercept = vec![0.75, 3.1, 2.9, 3.2, 4.0, 2.7, 3.4, 3.8];
    let compact_intercept = frozen_intercept[1..].to_vec();

    let output = fit.frozen_rlog(&counts, &frozen_intercept, 2.0).unwrap();
    let compact_output = compact_fit
        .frozen_rlog(&compact_counts, &compact_intercept, 2.0)
        .unwrap();

    assert_eq!(output.intercept, frozen_intercept);
    assert!(
        output
            .transformed
            .row(0)
            .unwrap()
            .iter()
            .all(|value| *value == frozen_intercept[0])
    );
    for compact_gene in 0..compact_counts.n_genes() {
        assert_eq!(
            output.transformed.row(compact_gene + 1).unwrap(),
            compact_output.transformed.row(compact_gene).unwrap()
        );
    }
}

#[test]
fn builder_rlog_glm_mu_matches_manual_fit_state_transform() {
    let full_counts = native_wald_counts_with_zero_row();
    let counts = full_counts.select_rows(&[1, 2, 3, 4, 5, 6, 7]).unwrap();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);

    let output = builder.rlog_glm_mu(&counts, &design).unwrap();
    let diagnostic_output = builder.rlog_glm_mu_with_fit(&counts, &design).unwrap();
    let fit = builder
        .fit_map_dispersions_glm_mu(&counts, &design)
        .unwrap();
    let expected = fit.rlog(&counts).unwrap();

    assert_eq!(output.metadata(), expected.metadata());
    assert_eq!(
        output.transformed.as_slice(),
        expected.transformed.as_slice()
    );
    assert_eq!(diagnostic_output.rlog.metadata(), expected.metadata());
    assert_eq!(diagnostic_output.design_mode, RlogDesignMode::DesignAware);
    assert_eq!(diagnostic_output.metadata().design_mode, "designAware");
    assert_eq!(diagnostic_output.metadata().fit_rows, counts.n_genes());
    assert_eq!(diagnostic_output.metadata().fit_cols, counts.n_samples());
    assert_eq!(diagnostic_output.metadata().trend_fit_type, Some("mean"));
}

#[test]
fn builder_frozen_rlog_glm_mu_matches_manual_fit_state_reuse() {
    let full_counts = native_wald_counts_with_zero_row();
    let counts = full_counts.select_rows(&[1, 2, 3, 4, 5, 6, 7]).unwrap();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);

    let output = builder
        .frozen_rlog_glm_mu_with_fit(&counts, &design)
        .unwrap();
    let source = builder.rlog_glm_mu_with_fit(&counts, &design).unwrap();
    let expected = source
        .fit
        .frozen_rlog(
            &counts,
            &source.rlog.intercept,
            source.rlog.sample_prior_variance,
        )
        .unwrap();

    assert_eq!(output.source_rlog.metadata(), source.rlog.metadata());
    assert_eq!(
        output.frozen_rlog.transformed.as_slice(),
        expected.transformed.as_slice()
    );
    assert_eq!(output.frozen_rlog.intercept, source.rlog.intercept);
    assert_eq!(output.design_mode, RlogDesignMode::DesignAware);
    assert_eq!(output.metadata().design_mode, "designAware");
    assert_eq!(output.metadata().fit_rows, counts.n_genes());
    assert_eq!(output.metadata().fit_cols, counts.n_samples());
    assert_eq!(output.metadata().trend_fit_type, Some("mean"));
}

#[test]
fn builder_blind_rlog_glm_mu_uses_intercept_only_design() {
    let full_counts = native_wald_counts_with_zero_row();
    let counts = full_counts.select_rows(&[1, 2, 3, 4, 5, 6, 7]).unwrap();
    let blind_design = DesignMatrix::intercept_only(counts.n_samples()).unwrap();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);

    let output = builder.blind_rlog_glm_mu(&counts).unwrap();
    let diagnostic_output = builder.blind_rlog_glm_mu_with_fit(&counts).unwrap();
    let expected = builder.rlog_glm_mu(&counts, &blind_design).unwrap();

    assert_eq!(output.metadata(), expected.metadata());
    assert_eq!(
        output.transformed.as_slice(),
        expected.transformed.as_slice()
    );
    assert_eq!(diagnostic_output.rlog.metadata(), expected.metadata());
    assert_eq!(diagnostic_output.design_mode, RlogDesignMode::Blind);
    assert_eq!(diagnostic_output.metadata().design_mode, "blind");
    assert_eq!(diagnostic_output.fit.design, Some(blind_design));
}

#[test]
fn builder_blind_frozen_rlog_glm_mu_uses_intercept_only_design() {
    let full_counts = native_wald_counts_with_zero_row();
    let counts = full_counts.select_rows(&[1, 2, 3, 4, 5, 6, 7]).unwrap();
    let blind_design = DesignMatrix::intercept_only(counts.n_samples()).unwrap();
    let builder = glm_mu_native_wald_builder().fit_type(FitType::Mean);

    let output = builder.blind_frozen_rlog_glm_mu_with_fit(&counts).unwrap();
    let expected = builder
        .frozen_rlog_glm_mu_with_fit(&counts, &blind_design)
        .unwrap();

    assert_eq!(
        output.source_rlog.metadata(),
        expected.source_rlog.metadata()
    );
    assert_eq!(
        output.frozen_rlog.transformed.as_slice(),
        expected.frozen_rlog.transformed.as_slice()
    );
    assert_eq!(output.design_mode, RlogDesignMode::Blind);
    assert_eq!(output.metadata().design_mode, "blind");
    assert_eq!(output.fit.design, Some(blind_design));
}
