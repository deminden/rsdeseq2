use super::*;

#[test]
fn beta_prior_refit_matches_manual_ridge_options() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 80, 80]).unwrap();
    let design = DesignMatrix::from_row_major(
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
    .unwrap();
    let beta_prior_variance = vec![1e6, 0.25];
    let ridge = beta_prior_variance_to_ridge_lambda(&beta_prior_variance).unwrap();

    let refit = fit_glms_with_beta_prior_variance(
        &counts,
        &design,
        &[1.0; 4],
        &[0.05],
        &beta_prior_variance,
        no_ridge_options(),
    )
    .unwrap();
    let manual = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0; 4],
        &[0.05],
        no_ridge_options().ridge_lambda_by_coefficient(ridge),
    )
    .unwrap();
    let mle = fit_fixed_dispersion_irls(&counts, &design, &[1.0; 4], &[0.05], no_ridge_options())
        .unwrap();

    assert_eq!(refit, manual);
    assert!(refit.beta.as_slice()[1].abs() < mle.beta.as_slice()[1].abs());
}

#[test]
fn beta_prior_refit_matches_deseq2_reference() {
    let rows = read_reference_tsv("beta_prior_refit_reference.tsv");
    let mu_rows = read_reference_tsv("beta_prior_refit_mu.tsv");
    let hat_rows = read_reference_tsv("beta_prior_refit_hat.tsv");
    let prior_rows = read_reference_tsv("beta_prior_variance_reference.tsv");
    let fixed_dispersions = read_fixed_dispersions();
    let size_factors = read_size_factors("size_factors_ratio.tsv");

    let beta_prior_variance = prior_rows
        .iter()
        .map(|row| parse_required_f64(row, "weighted"))
        .collect::<Vec<_>>();
    let fit = fit_glms_with_beta_prior_variance(
        &reference_counts(),
        &reference_full_design(),
        &size_factors,
        &fixed_dispersions,
        &beta_prior_variance,
        no_ridge_options(),
    )
    .unwrap();
    let samples = reference_sample_names();

    assert_eq!(rows.len(), fit.beta.n_rows());
    for (gene, row) in rows.iter().enumerate() {
        assert_float_close(
            fit.beta.row(gene).unwrap()[0],
            parse_required_f64(row, "beta_intercept"),
            1e-5,
            1e-5,
            &format!("beta prior refit beta intercept gene {gene}"),
        );
        assert_float_close(
            fit.beta.row(gene).unwrap()[1],
            parse_required_f64(row, "beta_conditionB"),
            1e-5,
            1e-5,
            &format!("beta prior refit beta conditionB gene {gene}"),
        );
        assert_float_close(
            fit.beta_se.row(gene).unwrap()[0],
            parse_required_f64(row, "beta_se_intercept"),
            1e-5,
            1e-5,
            &format!("beta prior refit beta SE intercept gene {gene}"),
        );
        assert_float_close(
            fit.beta_se.row(gene).unwrap()[1],
            parse_required_f64(row, "beta_se_conditionB"),
            1e-5,
            1e-5,
            &format!("beta prior refit beta SE conditionB gene {gene}"),
        );
        assert_float_close(
            fit.log_like[gene],
            parse_required_f64(row, "log_like"),
            1e-5,
            1e-5,
            &format!("beta prior refit log-likelihood gene {gene}"),
        );
        assert_eq!(
            fit.beta_converged[gene],
            parse_required_bool(row, "converged"),
            "beta prior refit convergence gene {gene}"
        );
        assert_eq!(
            fit.beta_iter[gene],
            parse_required_f64(row, "iterations") as usize,
            "beta prior refit iterations gene {gene}"
        );
    }

    for (gene, (mu_row, hat_row)) in mu_rows.iter().zip(hat_rows.iter()).enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_float_close(
                fit.mu.row(gene).unwrap()[sample],
                parse_required_f64(mu_row, sample_name),
                1e-5,
                1e-5,
                &format!("beta prior refit mu gene {gene} sample {sample}"),
            );
            assert_float_close(
                fit.hat_diagonal.row(gene).unwrap()[sample],
                parse_required_f64(hat_row, sample_name),
                1e-5,
                1e-5,
                &format!("beta prior refit hat gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn beta_prior_estimated_refit_matches_deseq2_reference() {
    let prior_rows = read_reference_tsv("beta_prior_variance_reference.tsv");
    let refit_rows = read_reference_tsv("beta_prior_refit_reference.tsv");
    let fixed_dispersions = read_fixed_dispersions();
    let size_factors = read_size_factors("size_factors_ratio.tsv");
    let base_metadata = read_reference_tsv("base_metadata_ratio.tsv");

    let base_mean = base_metadata
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let fit = fit_glms_with_estimated_beta_prior_variance(
        &reference_counts(),
        &reference_full_design(),
        &size_factors,
        &fixed_dispersions,
        &base_mean,
        &fixed_dispersions,
        BetaPriorRefitOptions {
            fit_options: no_ridge_options(),
            variance_options: BetaPriorVarianceOptions {
                method: BetaPriorVarianceMethod::Weighted,
                upper_quantile: 0.05,
                ..BetaPriorVarianceOptions::default()
            },
        },
    )
    .unwrap();

    assert_eq!(prior_rows.len(), fit.beta_prior_variance.len());
    for (idx, row) in prior_rows.iter().enumerate() {
        assert_float_close(
            fit.beta_prior_variance[idx],
            parse_required_f64(row, "weighted"),
            1e-5,
            1e-5,
            &format!("estimated beta prior variance coefficient {idx}"),
        );
    }

    assert_eq!(refit_rows.len(), fit.prior_fit.beta.n_rows());
    for (gene, row) in refit_rows.iter().enumerate() {
        assert_float_close(
            fit.prior_fit.beta.row(gene).unwrap()[0],
            parse_required_f64(row, "beta_intercept"),
            1e-5,
            1e-5,
            &format!("estimated beta prior refit beta intercept gene {gene}"),
        );
        assert_float_close(
            fit.prior_fit.beta.row(gene).unwrap()[1],
            parse_required_f64(row, "beta_conditionB"),
            1e-5,
            1e-5,
            &format!("estimated beta prior refit beta conditionB gene {gene}"),
        );
        assert_float_close(
            fit.prior_fit.beta_se.row(gene).unwrap()[1],
            parse_required_f64(row, "beta_se_conditionB"),
            1e-5,
            1e-5,
            &format!("estimated beta prior refit beta SE conditionB gene {gene}"),
        );
        assert_float_close(
            fit.prior_fit.log_like[gene],
            parse_required_f64(row, "log_like"),
            1e-5,
            1e-5,
            &format!("estimated beta prior refit log-likelihood gene {gene}"),
        );
    }
}

#[test]
fn beta_prior_refit_with_weights_matches_normalization_factor_path() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 80, 80]).unwrap();
    let design = DesignMatrix::from_row_major(
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
    .unwrap();
    let normalization_factors = RowMajorMatrix::from_row_major(1, 4, vec![1.0; 4]).unwrap();
    let weights = RowMajorMatrix::from_row_major(1, 4, vec![1.0; 4]).unwrap();
    let beta_prior_variance = vec![1e6, 0.25];

    let weighted = fit_glms_with_beta_prior_variance_and_weights(
        &counts,
        &design,
        &[1.0; 4],
        &[0.05],
        Some(&weights),
        &beta_prior_variance,
        no_ridge_options(),
    )
    .unwrap();
    let normalization = fit_glms_with_beta_prior_variance_and_normalization_factors_and_weights(
        &counts,
        &design,
        &normalization_factors,
        &[0.05],
        Some(&weights),
        &beta_prior_variance,
        no_ridge_options(),
    )
    .unwrap();

    assert_eq!(weighted, normalization);
}

#[test]
fn beta_prior_refit_weight_helpers_validate_inputs() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 80, 80]).unwrap();
    let design = DesignMatrix::from_row_major(
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
    .unwrap();
    let beta_prior_variance = vec![1e6, 0.25];
    let bad_weights = RowMajorMatrix::from_row_major(1, 4, vec![1.0, 1.0, -0.5, 1.0]).unwrap();
    let bad_normalization_factors =
        RowMajorMatrix::from_row_major(1, 4, vec![1.0, 1.0, f64::NAN, 1.0]).unwrap();

    assert!(
        fit_glms_with_beta_prior_variance_and_weights(
            &counts,
            &design,
            &[1.0; 4],
            &[0.05],
            Some(&bad_weights),
            &beta_prior_variance,
            no_ridge_options(),
        )
        .is_err()
    );
    assert!(
        fit_glms_with_beta_prior_variance_and_normalization_factors_and_weights(
            &counts,
            &design,
            &bad_normalization_factors,
            &[0.05],
            None,
            &beta_prior_variance,
            no_ridge_options(),
        )
        .is_err()
    );
}

#[test]
fn beta_prior_estimated_refit_runs_mle_then_shrunk_fit() {
    let counts = CountMatrix::from_row_major_u32(
        3,
        4,
        vec![
            10, 10, 20, 20, //
            20, 20, 80, 80, //
            12, 12, 12, 12,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
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
    .unwrap();
    let fit = fit_glms_with_estimated_beta_prior_variance(
        &counts,
        &design,
        &[1.0; 4],
        &[0.05; 3],
        &[15.0, 50.0, 12.0],
        &[0.1; 3],
        BetaPriorRefitOptions {
            fit_options: no_ridge_options(),
            variance_options: BetaPriorVarianceOptions {
                method: BetaPriorVarianceMethod::Quantile,
                upper_quantile: 0.5,
                ..BetaPriorVarianceOptions::default()
            },
        },
    )
    .unwrap();

    assert_eq!(fit.beta_prior_variance.len(), 2);
    assert_relative_eq!(fit.beta_prior_variance[0], 1e6, epsilon = 1e-12);
    assert_relative_eq!(
        fit.beta_prior_variance[1],
        normal_matched_variance(1.0),
        epsilon = 1e-8,
        max_relative = 1e-8
    );
    assert!(fit.prior_fit.beta.as_slice()[3].abs() < fit.mle_fit.beta.as_slice()[3].abs());
}

#[test]
fn beta_prior_estimated_refit_accepts_normalization_factors_and_weights() {
    let counts = CountMatrix::from_row_major_u32(
        3,
        4,
        vec![
            10, 10, 20, 20, //
            20, 20, 80, 80, //
            12, 12, 12, 12,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
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
    .unwrap();
    let normalization_factors = RowMajorMatrix::from_row_major(3, 4, vec![1.0; 12]).unwrap();
    let weights = RowMajorMatrix::from_row_major(3, 4, vec![1.0; 12]).unwrap();
    let options = BetaPriorRefitOptions {
        fit_options: no_ridge_options(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };

    let size_factor_fit = fit_glms_with_estimated_beta_prior_variance(
        &counts,
        &design,
        &[1.0; 4],
        &[0.05; 3],
        &[15.0, 50.0, 12.0],
        &[0.1; 3],
        options.clone(),
    )
    .unwrap();
    let normalization_fit = fit_glms_with_estimated_beta_prior_variance_and_normalization_factors(
        &counts,
        &design,
        &normalization_factors,
        &[0.05; 3],
        &[15.0, 50.0, 12.0],
        &[0.1; 3],
        options.clone(),
    )
    .unwrap();
    let weighted_fit = fit_glms_with_estimated_beta_prior_variance_and_weights(
        &counts,
        &design,
        BetaPriorSizeFactorWeightInput {
            size_factors: &[1.0; 4],
            weights: Some(&weights),
        },
        &[0.05; 3],
        &[15.0, 50.0, 12.0],
        &[0.1; 3],
        options.clone(),
    )
    .unwrap();
    let weighted_normalization_fit =
        fit_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
            &counts,
            &design,
            BetaPriorNormalizationFactorWeightInput {
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
            },
            &[0.05; 3],
            &[15.0, 50.0, 12.0],
            &[0.1; 3],
            options,
        )
        .unwrap();

    assert_eq!(
        normalization_fit.beta_prior_variance,
        size_factor_fit.beta_prior_variance
    );
    assert_eq!(normalization_fit.mle_fit, size_factor_fit.mle_fit);
    assert_eq!(normalization_fit.prior_fit, size_factor_fit.prior_fit);
    assert_eq!(
        weighted_fit.beta_prior_variance,
        size_factor_fit.beta_prior_variance
    );
    assert_eq!(weighted_fit.mle_fit, size_factor_fit.mle_fit);
    assert_eq!(weighted_fit.prior_fit, size_factor_fit.prior_fit);
    assert_eq!(
        weighted_normalization_fit.beta_prior_variance,
        weighted_fit.beta_prior_variance
    );
    assert_eq!(weighted_normalization_fit.mle_fit, weighted_fit.mle_fit);
    assert_eq!(weighted_normalization_fit.prior_fit, weighted_fit.prior_fit);
}

#[test]
fn beta_prior_estimated_refit_weights_influence_mle_and_prior_fit() {
    let counts = CountMatrix::from_row_major_u32(
        3,
        4,
        vec![
            10, 10, 20, 20, //
            20, 20, 80, 80, //
            12, 12, 12, 12,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
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
    .unwrap();
    let weights = RowMajorMatrix::from_row_major(
        3,
        4,
        vec![
            1.0, 1.0, 0.25, 1.0, //
            1.0, 0.25, 1.0, 1.0, //
            1.0, 1.0, 1.0, 1.0,
        ],
    )
    .unwrap();
    let options = BetaPriorRefitOptions {
        fit_options: no_ridge_options(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };

    let unweighted = fit_glms_with_estimated_beta_prior_variance(
        &counts,
        &design,
        &[1.0; 4],
        &[0.05; 3],
        &[15.0, 50.0, 12.0],
        &[0.1; 3],
        options.clone(),
    )
    .unwrap();
    let weighted = fit_glms_with_estimated_beta_prior_variance_and_weights(
        &counts,
        &design,
        BetaPriorSizeFactorWeightInput {
            size_factors: &[1.0; 4],
            weights: Some(&weights),
        },
        &[0.05; 3],
        &[15.0, 50.0, 12.0],
        &[0.1; 3],
        options,
    )
    .unwrap();

    assert_ne!(weighted.mle_fit.beta, unweighted.mle_fit.beta);
    assert_ne!(weighted.prior_fit.beta, unweighted.prior_fit.beta);
    assert!(
        weighted
            .mle_fit
            .beta
            .as_slice()
            .iter()
            .all(|value| value.is_finite())
    );
    assert!(
        weighted
            .prior_fit
            .beta
            .as_slice()
            .iter()
            .all(|value| value.is_finite())
    );
}

#[test]
fn beta_prior_estimated_refit_weight_helpers_validate_inputs() {
    let counts = CountMatrix::from_row_major_u32(
        3,
        4,
        vec![
            10, 10, 20, 20, //
            20, 20, 80, 80, //
            12, 12, 12, 12,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
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
    .unwrap();
    let bad_weights = RowMajorMatrix::from_row_major(2, 4, vec![1.0; 8]).unwrap();
    let bad_normalization_factors = RowMajorMatrix::from_row_major(
        3,
        4,
        vec![
            1.0, 1.0, 1.0, 1.0, //
            1.0, 0.0, 1.0, 1.0, //
            1.0, 1.0, 1.0, 1.0,
        ],
    )
    .unwrap();
    let options = BetaPriorRefitOptions {
        fit_options: no_ridge_options(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };

    assert!(
        fit_glms_with_estimated_beta_prior_variance_and_weights(
            &counts,
            &design,
            BetaPriorSizeFactorWeightInput {
                size_factors: &[1.0; 4],
                weights: Some(&bad_weights),
            },
            &[0.05; 3],
            &[15.0, 50.0, 12.0],
            &[0.1; 3],
            options.clone(),
        )
        .is_err()
    );
    assert!(
        fit_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
            &counts,
            &design,
            BetaPriorNormalizationFactorWeightInput {
                normalization_factors: &bad_normalization_factors,
                weights: None,
            },
            &[0.05; 3],
            &[15.0, 50.0, 12.0],
            &[0.1; 3],
            options,
        )
        .is_err()
    );
}
