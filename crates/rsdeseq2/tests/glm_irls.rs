use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

mod common;
use common::*;

fn no_ridge_options() -> IrlsOptions {
    IrlsOptions {
        ridge_lambda: 0.0,
        ..IrlsOptions::default()
    }
}

fn no_ridge_qr_options() -> IrlsOptions {
    IrlsOptions {
        solver: IrlsSolver::Qr,
        ..no_ridge_options()
    }
}

#[test]
fn irls_intercept_only_agrees_with_intercept_shortcut_for_equal_size_factors() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![4, 4, 8, 8]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        1,
        vec![1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into()]),
    )
    .unwrap();

    let irls = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.1],
        no_ridge_options(),
    )
    .unwrap();
    let shortcut =
        fit_intercept_only_fixed_dispersion(&counts, &[1.0, 1.0, 1.0, 1.0], &[0.1]).unwrap();

    assert!(irls.beta_converged[0]);
    assert_relative_eq!(
        irls.beta.as_slice()[0],
        shortcut.beta.as_slice()[0],
        epsilon = 1e-8
    );
    for (actual, expected) in irls.mu.as_slice().iter().zip(shortcut.mu.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-7);
    }
}

#[test]
fn fit_irls_dispatches_intercept_only_shortcut() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![6, 6, 12, 12]).unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        1,
        vec![1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into()]),
    )
    .unwrap();

    let wrapped = fit_irls(
        &counts,
        &design,
        &[1.0, 1.0, 2.0, 2.0],
        &[0.2],
        IrlsOptions::default(),
    )
    .unwrap();
    let shortcut =
        fit_intercept_only_fixed_dispersion(&counts, &[1.0, 1.0, 2.0, 2.0], &[0.2]).unwrap();

    assert_eq!(wrapped, shortcut);
}

#[test]
fn irls_two_group_design_recovers_known_log2_fold_change() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
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

    let fit = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        no_ridge_options(),
    )
    .unwrap();

    assert!(fit.beta_converged[0]);
    assert_relative_eq!(fit.beta.as_slice()[0], 10.0_f64.log2(), epsilon = 1e-8);
    assert_relative_eq!(fit.beta.as_slice()[1], 2.0_f64.log2(), epsilon = 1e-8);
    assert_eq!(fit.n_terms, 2);
    assert_eq!(
        fit.model_matrix.coefficient_names().unwrap()[1],
        "condition_B_vs_A"
    );
    assert_relative_eq!(fit.mu.as_slice()[0], 10.0, epsilon = 1e-7);
    assert_relative_eq!(fit.mu.as_slice()[2], 20.0, epsilon = 1e-7);
    assert!(fit.beta_iter[0] < no_ridge_options().maxit);
    let beta_covariance = fit.beta_covariance.as_ref().unwrap();
    assert_eq!(beta_covariance.n_cols(), 4);
    assert_relative_eq!(
        beta_covariance.row(0).unwrap()[0].sqrt(),
        fit.beta_se.row(0).unwrap()[0],
        epsilon = 1e-12
    );
}

#[test]
fn fit_irls_dispatches_general_irls_for_multi_coefficient_design() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
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
    let options = no_ridge_options();

    let wrapped = fit_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        options.clone(),
    )
    .unwrap();
    let direct =
        fit_fixed_dispersion_irls(&counts, &design, &[1.0, 1.0, 1.0, 1.0], &[0.05], options)
            .unwrap();

    assert_eq!(wrapped, direct);
}

#[test]
fn estimate_beta_wraps_fixed_dispersion_beta_dispatch() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
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
    let options = no_ridge_options();

    let beta_fit = estimate_beta(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        options.clone(),
    )
    .unwrap();
    let irls_fit = fit_irls(&counts, &design, &[1.0, 1.0, 1.0, 1.0], &[0.05], options).unwrap();

    assert_eq!(beta_fit, irls_fit);
}

fn normal_75_quantile() -> f64 {
    0.674_489_750_196_081_7
}

fn normal_matched_variance(abs_quantile: f64) -> f64 {
    (abs_quantile / normal_75_quantile()).powi(2)
}

#[test]
fn beta_prior_quantile_matches_type7_upper_abs_quantile() {
    let variance = match_upper_quantile_for_variance(&[1.0, -2.0, 3.0], 0.5).unwrap();

    assert_relative_eq!(
        variance,
        normal_matched_variance(2.0),
        epsilon = 1e-12,
        max_relative = 1e-12
    );
}

#[test]
fn beta_prior_weighted_quantile_uses_row_weights() {
    let variance =
        match_weighted_upper_quantile_for_variance(&[1.0, 10.0, 20.0], &[98.0, 1.0, 1.0], 0.5)
            .unwrap();

    assert_relative_eq!(
        variance,
        normal_matched_variance(1.0),
        epsilon = 1e-12,
        max_relative = 1e-12
    );
}

#[test]
fn beta_prior_weighted_quantile_rejects_overflowed_weight_sum() {
    let err = match_weighted_upper_quantile_for_variance(
        &[1.0, 10.0, 20.0],
        &[f64::MAX, f64::MAX, 1.0],
        0.5,
    )
    .unwrap_err();

    assert!(err.to_string().contains("finite total weight"));
}

#[test]
fn beta_prior_variance_handles_extreme_mean_dispersion_weights() {
    let betas = RowMajorMatrix::from_row_major(3, 1, vec![1.0, 2.0, 3.0]).unwrap();
    let variance = estimate_beta_prior_variance(
        &betas,
        &[1.0e200, 2.0, 3.0],
        &[1.0e200, 0.5, 0.25],
        Some(&["conditionB".to_string()]),
        BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Weighted,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    )
    .unwrap();

    assert_eq!(variance.len(), 1);
    assert!(variance[0].is_finite());
}

#[test]
fn beta_prior_variance_filters_large_betas_and_sets_intercept_wide() {
    let betas = RowMajorMatrix::from_row_major(
        4,
        2,
        vec![
            2.0,
            1.0, //
            2.0,
            2.0, //
            2.0,
            100.0, //
            2.0,
            f64::INFINITY,
        ],
    )
    .unwrap();
    let names = vec!["Intercept".to_string(), "condition_B_vs_A".to_string()];
    let options = BetaPriorVarianceOptions {
        method: BetaPriorVarianceMethod::Quantile,
        upper_quantile: 0.5,
        ..BetaPriorVarianceOptions::default()
    };

    let variance = estimate_beta_prior_variance(
        &betas,
        &[10.0, 20.0, 30.0, 40.0],
        &[0.1; 4],
        Some(&names),
        options,
    )
    .unwrap();

    assert_relative_eq!(variance[0], 1e6, epsilon = 1e-12);
    assert_relative_eq!(
        variance[1],
        normal_matched_variance(1.5),
        epsilon = 1e-12,
        max_relative = 1e-12
    );
}

#[test]
fn beta_prior_variance_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("beta_prior_variance_reference.tsv") else {
        return;
    };
    let Some(fixed_dispersions) = read_fixed_dispersions() else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };
    let Some(base_metadata) = read_optional_tsv("base_metadata_ratio.tsv") else {
        return;
    };

    let fit = fit_fixed_dispersion_irls(
        &reference_counts(),
        &reference_full_design(),
        &size_factors,
        &fixed_dispersions,
        no_ridge_options(),
    )
    .unwrap();
    let base_mean = base_metadata
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let names = vec!["Intercept".to_string(), "conditionB".to_string()];

    let weighted = estimate_beta_prior_variance(
        &fit.beta,
        &base_mean,
        &fixed_dispersions,
        Some(&names),
        BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Weighted,
            upper_quantile: 0.05,
            ..BetaPriorVarianceOptions::default()
        },
    )
    .unwrap();
    let quantile = estimate_beta_prior_variance(
        &fit.beta,
        &base_mean,
        &fixed_dispersions,
        Some(&names),
        BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.05,
            ..BetaPriorVarianceOptions::default()
        },
    )
    .unwrap();

    assert_eq!(rows.len(), names.len());
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(
            row.get("coefficient").map(String::as_str),
            Some(names[idx].as_str())
        );
        assert_float_close(
            weighted[idx],
            parse_required_f64(row, "weighted"),
            1e-5,
            1e-5,
            &format!("beta prior weighted variance coefficient {idx}"),
        );
        assert_float_close(
            quantile[idx],
            parse_required_f64(row, "quantile"),
            1e-5,
            1e-5,
            &format!("beta prior quantile variance coefficient {idx}"),
        );
    }
}

#[test]
fn beta_prior_variance_validates_dimensions_and_options() {
    let betas = RowMajorMatrix::from_row_major(2, 1, vec![1.0, 2.0]).unwrap();

    assert!(estimate_beta_prior_variance(
        &betas,
        &[10.0],
        &[0.1, 0.2],
        None,
        BetaPriorVarianceOptions::default(),
    )
    .is_err());

    assert!(estimate_beta_prior_variance(
        &betas,
        &[10.0, 20.0],
        &[0.1, 0.2],
        None,
        BetaPriorVarianceOptions {
            upper_quantile: 1.0,
            ..BetaPriorVarianceOptions::default()
        },
    )
    .is_err());
}

#[test]
fn beta_prior_variance_to_ridge_lambda_matches_deseq2_scale() {
    let ridge = beta_prior_variance_to_ridge_lambda(&[1e6, 0.25]).unwrap();

    assert_relative_eq!(
        ridge[0],
        1.0 / (1e6 * std::f64::consts::LN_2.powi(2)),
        epsilon = 1e-18,
        max_relative = 1e-12
    );
    assert_relative_eq!(
        ridge[1],
        1.0 / (0.25 * std::f64::consts::LN_2.powi(2)),
        epsilon = 1e-14,
        max_relative = 1e-12
    );
    assert!(beta_prior_variance_to_ridge_lambda(&[0.0]).is_err());
}

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
fn beta_prior_refit_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("beta_prior_refit_reference.tsv") else {
        return;
    };
    let Some(mu_rows) = read_optional_tsv("beta_prior_refit_mu.tsv") else {
        return;
    };
    let Some(hat_rows) = read_optional_tsv("beta_prior_refit_hat.tsv") else {
        return;
    };
    let Some(prior_rows) = read_optional_tsv("beta_prior_variance_reference.tsv") else {
        return;
    };
    let Some(fixed_dispersions) = read_fixed_dispersions() else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };

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
fn beta_prior_estimated_refit_matches_optional_deseq2_reference() {
    let Some(prior_rows) = read_optional_tsv("beta_prior_variance_reference.tsv") else {
        return;
    };
    let Some(refit_rows) = read_optional_tsv("beta_prior_refit_reference.tsv") else {
        return;
    };
    let Some(fixed_dispersions) = read_fixed_dispersions() else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };
    let Some(base_metadata) = read_optional_tsv("base_metadata_ratio.tsv") else {
        return;
    };

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

    assert!(fit_glms_with_beta_prior_variance_and_weights(
        &counts,
        &design,
        &[1.0; 4],
        &[0.05],
        Some(&bad_weights),
        &beta_prior_variance,
        no_ridge_options(),
    )
    .is_err());
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
    assert!(weighted
        .mle_fit
        .beta
        .as_slice()
        .iter()
        .all(|value| value.is_finite()));
    assert!(weighted
        .prior_fit
        .beta
        .as_slice()
        .iter()
        .all(|value| value.is_finite()));
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

    assert!(fit_glms_with_estimated_beta_prior_variance_and_weights(
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
    .is_err());
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

#[test]
fn fit_with_dispersion_wraps_fixed_dispersion_irls() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
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
    let options = no_ridge_options();

    let wrapped = fit_with_dispersion(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        options.clone(),
    )
    .unwrap();
    let direct =
        fit_fixed_dispersion_irls(&counts, &design, &[1.0, 1.0, 1.0, 1.0], &[0.05], options)
            .unwrap();

    assert_eq!(wrapped, direct);
}

#[test]
fn irls_optim_fallback_refits_nonconverged_rows_when_enabled() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
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

    let without_optim = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            maxit: 1,
            ridge_lambda: 0.0,
            use_optim: false,
            ..IrlsOptions::default()
        },
    )
    .unwrap();
    let with_optim = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            maxit: 1,
            ridge_lambda: 0.0,
            use_optim: true,
            ..IrlsOptions::default()
        },
    )
    .unwrap();

    assert!(!without_optim.beta_converged[0]);
    assert!(with_optim.beta_converged[0]);
    assert_relative_eq!(
        with_optim.beta.as_slice()[0],
        10.0_f64.log2(),
        epsilon = 1e-5
    );
    assert_relative_eq!(
        with_optim.beta.as_slice()[1],
        2.0_f64.log2(),
        epsilon = 1e-5
    );
    assert!(with_optim.log_like[0] >= without_optim.log_like[0] - 1e-8);
}

#[test]
fn irls_optim_fallback_refreshes_hat_diagonal_for_refit_rows() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
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

    let expected = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        no_ridge_options(),
    )
    .unwrap();
    let with_optim = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            maxit: 1,
            ridge_lambda: 0.0,
            use_optim: true,
            ..IrlsOptions::default()
        },
    )
    .unwrap();

    assert!(with_optim.beta_converged[0]);
    for (actual, expected) in with_optim
        .hat_diagonal
        .as_slice()
        .iter()
        .zip(expected.hat_diagonal.as_slice())
    {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
}

#[test]
fn irls_optim_fallback_log_likelihood_uses_min_mu_floored_means() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![1, 1, 20, 20]).unwrap();
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

    let fit = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            maxit: 1,
            ridge_lambda: 0.0,
            use_optim: true,
            ..IrlsOptions::default()
        },
    )
    .unwrap();
    let stored_mu = fit.mu.row(0).unwrap();
    let inference_mu = stored_mu
        .iter()
        .copied()
        .map(|value| value.max(IrlsOptions::default().min_mu))
        .collect::<Vec<_>>();
    let expected_log_like =
        nbinom_log_likelihood(counts.row(0).unwrap(), &inference_mu, 0.05).unwrap();

    assert!(fit.beta_converged[0]);
    assert_relative_eq!(fit.log_like[0], expected_log_like, epsilon = 1e-12);
}

#[test]
fn irls_qr_solver_matches_normal_equations_for_two_group_design() {
    let counts = CountMatrix::from_row_major_u32(2, 4, vec![10, 10, 20, 20, 8, 14, 7, 15]).unwrap();
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

    let normal = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05, 0.2],
        no_ridge_options(),
    )
    .unwrap();
    let qr = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05, 0.2],
        no_ridge_qr_options(),
    )
    .unwrap();

    assert_eq!(normal.beta_converged, qr.beta_converged);
    for (actual, expected) in qr.beta.as_slice().iter().zip(normal.beta.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
    for (actual, expected) in qr.beta_se.as_slice().iter().zip(normal.beta_se.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
    for (actual, expected) in qr.mu.as_slice().iter().zip(normal.mu.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
}

#[test]
fn original_qr_results_match_non_qr_results_on_stable_rows() {
    let counts = CountMatrix::from_row_major_u32(2, 4, vec![10, 10, 20, 20, 8, 14, 7, 15]).unwrap();
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
    let normal = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05, 0.2],
        no_ridge_options(),
    )
    .unwrap();
    let qr = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05, 0.2],
        no_ridge_qr_options(),
    )
    .unwrap();

    let normal_wald = wald_test_coefficient(&normal, 1).unwrap();
    let qr_wald = wald_test_coefficient(&qr, 1).unwrap();

    for gene in 0..counts.n_genes() {
        assert_relative_eq!(
            qr.beta.row(gene).unwrap()[1],
            normal.beta.row(gene).unwrap()[1],
            epsilon = 1e-8,
            max_relative = 1e-8
        );
        assert_relative_eq!(
            qr.beta_se.row(gene).unwrap()[1],
            normal.beta_se.row(gene).unwrap()[1],
            epsilon = 1e-8,
            max_relative = 1e-8
        );
    }
    for (actual, expected) in qr_wald.stat.iter().zip(normal_wald.stat.iter()) {
        assert_relative_eq!(
            actual.unwrap(),
            expected.unwrap(),
            epsilon = 1e-8,
            max_relative = 1e-8
        );
    }
    for (actual, expected) in qr_wald.pvalue.iter().zip(normal_wald.pvalue.iter()) {
        assert_relative_eq!(
            actual.unwrap(),
            expected.unwrap(),
            epsilon = 1e-8,
            max_relative = 1e-8
        );
    }
}

#[test]
fn irls_qr_solver_handles_default_ridge() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
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

    let fit = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            solver: IrlsSolver::Qr,
            ..IrlsOptions::default()
        },
    )
    .unwrap();

    assert!(fit.beta_converged[0]);
    assert_relative_eq!(fit.beta.as_slice()[1], 2.0_f64.log2(), epsilon = 1e-5);
    assert!(fit
        .hat_diagonal
        .as_slice()
        .iter()
        .all(|value| value.is_finite()));
}

#[test]
fn irls_vector_ridge_matches_scalar_ridge_when_values_match() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
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

    let scalar = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            ridge_lambda: 0.25,
            ..IrlsOptions::default()
        },
    )
    .unwrap();
    let vector = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        }
        .ridge_lambda_by_coefficient(vec![0.25, 0.25]),
    )
    .unwrap();

    for (actual, expected) in vector.beta.as_slice().iter().zip(scalar.beta.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-10, max_relative = 1e-10);
    }
    for (actual, expected) in vector
        .beta_covariance
        .as_ref()
        .unwrap()
        .as_slice()
        .iter()
        .zip(scalar.beta_covariance.as_ref().unwrap().as_slice())
    {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-10, max_relative = 1e-10);
    }
}

#[test]
fn irls_vector_ridge_can_penalize_one_coefficient() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
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

    let unpenalized = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        no_ridge_options(),
    )
    .unwrap();
    let penalized = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        }
        .ridge_lambda_by_coefficient(vec![0.0, 25.0]),
    )
    .unwrap();

    assert!(penalized.beta_converged[0]);
    assert!(penalized.beta.as_slice()[1].abs() < unpenalized.beta.as_slice()[1].abs());
    assert_relative_eq!(
        penalized.beta.as_slice()[0],
        15.0_f64.log2(),
        epsilon = 2e-1
    );
}

#[test]
fn irls_qr_solver_matches_normal_equations_with_vector_ridge() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
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

    let normal = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        }
        .ridge_lambda_by_coefficient(vec![0.0, 25.0]),
    )
    .unwrap();
    let qr = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        IrlsOptions {
            ridge_lambda: 0.0,
            solver: IrlsSolver::Qr,
            ..IrlsOptions::default()
        }
        .ridge_lambda_by_coefficient(vec![0.0, 25.0]),
    )
    .unwrap();

    for (actual, expected) in qr.beta.as_slice().iter().zip(normal.beta.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
    for (actual, expected) in qr.mu.as_slice().iter().zip(normal.mu.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
}

#[test]
fn irls_all_one_weights_match_unweighted_fit() {
    let counts = CountMatrix::from_row_major_u32(2, 4, vec![10, 10, 20, 20, 8, 14, 7, 15]).unwrap();
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
    let weights = RowMajorMatrix::from_row_major(2, 4, vec![1.0; 8]).unwrap();

    let unweighted = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05, 0.2],
        no_ridge_options(),
    )
    .unwrap();
    let weighted = fit_fixed_dispersion_irls_with_weights(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05, 0.2],
        Some(&weights),
        no_ridge_options(),
    )
    .unwrap();

    assert_eq!(weighted.beta_converged, unweighted.beta_converged);
    for (actual, expected) in weighted
        .beta
        .as_slice()
        .iter()
        .zip(unweighted.beta.as_slice())
    {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
    for (actual, expected) in weighted.mu.as_slice().iter().zip(unweighted.mu.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
    for (actual, expected) in weighted.log_like.iter().zip(unweighted.log_like.iter()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
}

#[test]
fn irls_zero_weight_sample_does_not_influence_fit() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 1_000, 20, 20]).unwrap();
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
    let weights = RowMajorMatrix::from_row_major(1, 4, vec![1.0, 0.0, 1.0, 1.0]).unwrap();

    let fit = fit_fixed_dispersion_irls_with_weights(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        Some(&weights),
        no_ridge_options(),
    )
    .unwrap();

    assert!(fit.beta_converged[0]);
    assert_relative_eq!(fit.beta.as_slice()[0], 10.0_f64.log2(), epsilon = 1e-8);
    assert_relative_eq!(fit.beta.as_slice()[1], 2.0_f64.log2(), epsilon = 1e-8);
    assert_relative_eq!(fit.mu.as_slice()[1], 10.0, epsilon = 1e-7);
}

#[test]
fn irls_qr_solver_matches_weighted_normal_equations() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 1_000, 20, 20]).unwrap();
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
    let weights = RowMajorMatrix::from_row_major(1, 4, vec![1.0, 0.0, 1.0, 1.0]).unwrap();

    let normal = fit_fixed_dispersion_irls_with_weights(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        Some(&weights),
        no_ridge_options(),
    )
    .unwrap();
    let qr = fit_fixed_dispersion_irls_with_weights(
        &counts,
        &design,
        &[1.0, 1.0, 1.0, 1.0],
        &[0.05],
        Some(&weights),
        no_ridge_qr_options(),
    )
    .unwrap();

    assert_eq!(qr.beta_converged, normal.beta_converged);
    for (actual, expected) in qr.beta.as_slice().iter().zip(normal.beta.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
    for (actual, expected) in qr.mu.as_slice().iter().zip(normal.mu.as_slice()) {
        assert_relative_eq!(*actual, *expected, epsilon = 1e-8, max_relative = 1e-8);
    }
}

#[test]
fn irls_validates_inputs() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 8]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    assert!(
        fit_fixed_dispersion_irls(&counts, &design, &[1.0], &[0.1], no_ridge_options()).is_err()
    );
    assert!(
        fit_fixed_dispersion_irls(&counts, &design, &[1.0, 1.0, 1.0], &[], no_ridge_options())
            .is_err()
    );

    let bad_design = DesignMatrix::from_row_major(2, 1, vec![1.0, 1.0], None).unwrap();
    assert!(fit_fixed_dispersion_irls(
        &counts,
        &bad_design,
        &[1.0, 1.0, 1.0],
        &[0.1],
        no_ridge_options()
    )
    .is_err());

    let bad_weights = RowMajorMatrix::from_row_major(1, 3, vec![1.0, -1.0, 1.0]).unwrap();
    assert!(fit_fixed_dispersion_irls_with_weights(
        &counts,
        &design,
        &[1.0, 1.0, 1.0],
        &[0.1],
        Some(&bad_weights),
        no_ridge_options()
    )
    .is_err());

    assert!(fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0],
        &[0.1],
        no_ridge_options().ridge_lambda_by_coefficient(vec![0.0, 0.0])
    )
    .is_err());
    assert!(fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0],
        &[0.1],
        no_ridge_options().ridge_lambda_by_coefficient(vec![-1.0])
    )
    .is_err());
}

#[test]
fn irls_rejects_all_zero_rows() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![0, 0, 0]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();
    let err = fit_fixed_dispersion_irls(
        &counts,
        &design,
        &[1.0, 1.0, 1.0],
        &[0.1],
        no_ridge_options(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("all zero"));
}
