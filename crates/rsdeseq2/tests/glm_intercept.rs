use approx::assert_relative_eq;
use rsdeseq2::glm::nb::nbinom_log_likelihood;
use rsdeseq2::prelude::*;

#[test]
fn intercept_only_unweighted_fit_matches_deseq2_shortcut_formulas() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![2, 4, 6, 10, 10, 10]).unwrap();
    let fit = fit_intercept_only_fixed_dispersion(&counts, &[1.0, 1.0, 1.0], &[0.1, 0.2]).unwrap();

    assert_eq!(fit.n_terms, 1);
    assert_eq!(fit.beta_converged, vec![true, true]);
    assert_eq!(fit.beta_iter, vec![1, 1]);
    assert_eq!(fit.model_matrix.matrix().as_slice(), &[1.0, 1.0, 1.0]);
    assert_eq!(
        fit.model_matrix.coefficient_names().unwrap(),
        &["Intercept".to_string()]
    );

    assert_relative_eq!(fit.beta.as_slice()[0], 4.0_f64.log2(), epsilon = 1e-12);
    assert_relative_eq!(fit.beta.as_slice()[1], 10.0_f64.log2(), epsilon = 1e-12);
    assert_relative_eq!(fit.mu.as_slice()[0], 4.0, epsilon = 1e-12);
    assert_relative_eq!(fit.mu.as_slice()[1], 4.0, epsilon = 1e-12);
    assert_relative_eq!(fit.mu.as_slice()[2], 4.0, epsilon = 1e-12);

    let working_weight = (4.0_f64.recip() + 0.1).recip();
    let sigma = (3.0 * working_weight).recip();
    assert_relative_eq!(
        fit.beta_se.as_slice()[0],
        std::f64::consts::LOG2_E * sigma.sqrt(),
        epsilon = 1e-12
    );
    for value in fit.hat_diagonal.row(0).unwrap() {
        assert_relative_eq!(*value, 1.0 / 3.0, epsilon = 1e-12);
    }

    assert_relative_eq!(
        fit.log_like[0],
        nbinom_log_likelihood(&[2, 4, 6], &[4.0, 4.0, 4.0], 0.1).unwrap(),
        epsilon = 1e-12
    );
}

#[test]
fn intercept_only_fit_uses_size_factors_to_reconstruct_mu() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 8]).unwrap();
    let fit = fit_intercept_only_fixed_dispersion(&counts, &[1.0, 2.0, 4.0], &[0.1]).unwrap();

    assert_relative_eq!(fit.beta.as_slice()[0], 2.0_f64.log2(), epsilon = 1e-12);
    assert_eq!(fit.mu.as_slice(), &[2.0, 4.0, 8.0]);
}

#[test]
fn intercept_only_fit_keeps_large_finite_normalized_mean() {
    let counts = CountMatrix::from_row_major_u32(1, 5, vec![u32::MAX; 5]).unwrap();
    let normalization_factors = RowMajorMatrix::from_row_major(1, 5, vec![1e-298; 5]).unwrap();

    let fit = fit_intercept_only_fixed_dispersion_with_normalization_factors(
        &counts,
        &normalization_factors,
        &[0.1],
        None,
    )
    .unwrap();

    assert!(fit.beta.as_slice()[0].is_finite());
    for value in fit.mu.as_slice() {
        assert_relative_eq!(*value, f64::from(u32::MAX), max_relative = 1e-12);
    }
}

#[test]
fn intercept_only_weighted_fit_matches_deseq2_shortcut_formulas() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 100]).unwrap();
    let weights = RowMajorMatrix::from_row_major(1, 3, vec![1.0, 1.0, 0.0]).unwrap();
    let fit = fit_intercept_only_fixed_dispersion_with_weights(
        &counts,
        &[1.0, 1.0, 1.0],
        &[0.1],
        Some(&weights),
    )
    .unwrap();

    assert_relative_eq!(fit.beta.as_slice()[0], 3.0_f64.log2(), epsilon = 1e-12);
    assert_eq!(fit.mu.as_slice(), &[3.0, 3.0, 3.0]);
    assert_relative_eq!(fit.hat_diagonal.as_slice()[0], 0.5, epsilon = 1e-12);
    assert_relative_eq!(fit.hat_diagonal.as_slice()[1], 0.5, epsilon = 1e-12);
    assert_relative_eq!(fit.hat_diagonal.as_slice()[2], 0.0, epsilon = 1e-12);

    let expected_log_like =
        nbinom_log_likelihood_weighted(&[2, 4, 100], &[3.0, 3.0, 3.0], 0.1, Some(&[1.0, 1.0, 0.0]))
            .unwrap();
    assert_relative_eq!(fit.log_like[0], expected_log_like, epsilon = 1e-12);
}

#[test]
fn intercept_only_weighted_fit_keeps_large_finite_normalized_mean() {
    let counts = CountMatrix::from_row_major_u32(1, 5, vec![u32::MAX; 5]).unwrap();
    let normalization_factors = RowMajorMatrix::from_row_major(1, 5, vec![1e-298; 5]).unwrap();
    let weights = RowMajorMatrix::from_row_major(1, 5, vec![1.0; 5]).unwrap();

    let fit = fit_intercept_only_fixed_dispersion_with_normalization_factors(
        &counts,
        &normalization_factors,
        &[0.1],
        Some(&weights),
    )
    .unwrap();

    assert!(fit.beta.as_slice()[0].is_finite());
    for value in fit.mu.as_slice() {
        assert_relative_eq!(*value, f64::from(u32::MAX), max_relative = 1e-12);
    }
}

#[test]
fn intercept_only_fit_rejects_all_zero_rows() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![0, 0, 0]).unwrap();
    let err = fit_intercept_only_fixed_dispersion(&counts, &[1.0, 1.0, 1.0], &[0.1]).unwrap_err();
    assert!(err.to_string().contains("all zero"));
}

#[test]
fn intercept_only_fit_validates_inputs() {
    let counts = CountMatrix::from_row_major_u32(1, 3, vec![2, 4, 8]).unwrap();
    assert!(fit_intercept_only_fixed_dispersion(&counts, &[1.0], &[0.1]).is_err());
    assert!(fit_intercept_only_fixed_dispersion(&counts, &[1.0, 2.0, 4.0], &[]).is_err());
    assert!(fit_intercept_only_fixed_dispersion(&counts, &[1.0, 2.0, 4.0], &[0.0]).is_err());

    let weights = RowMajorMatrix::from_row_major(1, 2, vec![1.0, 1.0]).unwrap();
    assert!(
        fit_intercept_only_fixed_dispersion_with_weights(
            &counts,
            &[1.0, 2.0, 4.0],
            &[0.1],
            Some(&weights),
        )
        .is_err()
    );
}

#[test]
fn intercept_only_fit_rejects_nonfinite_normalized_mean() {
    let counts = CountMatrix::from_row_major_u32(1, 2, vec![u32::MAX, u32::MAX]).unwrap();
    let normalization_factors =
        RowMajorMatrix::from_row_major(1, 2, vec![f64::MIN_POSITIVE, f64::MIN_POSITIVE]).unwrap();

    let err = fit_intercept_only_fixed_dispersion_with_normalization_factors(
        &counts,
        &normalization_factors,
        &[0.1],
        None,
    )
    .unwrap_err();

    assert!(
        err.to_string()
            .contains("non-finite normalized intercept mean")
    );
}

#[test]
fn intercept_only_fit_rejects_nonfinite_reconstructed_mu() {
    let counts = CountMatrix::from_row_major_u32(1, 2, vec![4, 1]).unwrap();
    let normalization_factors = RowMajorMatrix::from_row_major(1, 2, vec![1.0, f64::MAX]).unwrap();

    let err = fit_intercept_only_fixed_dispersion_with_normalization_factors(
        &counts,
        &normalization_factors,
        &[0.1],
        None,
    )
    .unwrap_err();

    assert!(err.to_string().contains("non-finite fitted intercept mean"));
}

#[test]
fn intercept_only_fit_rejects_nonfinite_output_covariance() {
    let counts = CountMatrix::from_row_major_u32(1, 2, vec![1, 1]).unwrap();

    let err = fit_intercept_only_fixed_dispersion(&counts, &[1.0, 1.0], &[f64::MAX]).unwrap_err();

    assert!(
        err.to_string()
            .contains("non-finite intercept beta covariance")
    );
}
