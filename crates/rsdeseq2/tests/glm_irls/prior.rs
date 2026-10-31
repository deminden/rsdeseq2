use super::*;

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
fn beta_prior_quantile_rejects_overflowed_variance() {
    let err = match_upper_quantile_for_variance(&[f64::MAX], 0.5).unwrap_err();

    assert!(
        err.to_string()
            .contains("beta prior variance quantile produced non-finite variance")
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
fn beta_prior_weighted_quantile_rejects_overflowed_variance() {
    let err = match_weighted_upper_quantile_for_variance(&[f64::MAX], &[1.0], 0.5).unwrap_err();

    assert!(
        err.to_string()
            .contains("weighted beta prior variance quantile produced non-finite variance")
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
fn beta_prior_variance_matches_deseq2_reference() {
    let rows = read_reference_tsv("beta_prior_variance_reference.tsv");
    let fixed_dispersions = read_fixed_dispersions();
    let size_factors = read_size_factors("size_factors_ratio.tsv");
    let base_metadata = read_reference_tsv("base_metadata_ratio.tsv");

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

    assert!(
        estimate_beta_prior_variance(
            &betas,
            &[10.0],
            &[0.1, 0.2],
            None,
            BetaPriorVarianceOptions::default(),
        )
        .is_err()
    );

    assert!(
        estimate_beta_prior_variance(
            &betas,
            &[10.0, 20.0],
            &[0.1, 0.2],
            None,
            BetaPriorVarianceOptions {
                upper_quantile: 1.0,
                ..BetaPriorVarianceOptions::default()
            },
        )
        .is_err()
    );
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
