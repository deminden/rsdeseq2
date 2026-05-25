use approx::assert_relative_eq;
use rsdeseq2::glm::nb::{
    nbinom_log_likelihood, nbinom_log_likelihood_matrix, nbinom_log_likelihood_weighted,
    nbinom_log_pmf, nbinom_negative_twice_log_likelihood,
};
use rsdeseq2::prelude::*;
use statrs::function::gamma::ln_gamma;

#[test]
fn nb_log_pmf_matches_hand_formula() {
    let count = 3_u32;
    let mu = 2.5_f64;
    let dispersion = 0.2_f64;
    let size = 1.0 / dispersion;
    let expected = ln_gamma(f64::from(count) + size) - ln_gamma(size) - ln_gamma(4.0)
        + size * (size / (size + mu)).ln()
        + f64::from(count) * (mu / (size + mu)).ln();

    let actual = nbinom_log_pmf(count, mu, dispersion).unwrap();
    assert_relative_eq!(actual, expected, epsilon = 1e-12);
}

#[test]
fn row_log_likelihood_sums_sample_log_pmfs() {
    let counts = [0_u32, 2, 5];
    let mu = [1.0, 2.5, 4.0];
    let dispersion = 0.3;

    let expected = counts
        .iter()
        .copied()
        .zip(mu)
        .map(|(count, mu)| nbinom_log_pmf(count, mu, dispersion).unwrap())
        .sum::<f64>();
    let actual = nbinom_log_likelihood(&counts, &mu, dispersion).unwrap();
    assert_relative_eq!(actual, expected, epsilon = 1e-12);
}

#[test]
fn weighted_log_likelihood_matches_deseq2_row_sum_shape() {
    let counts = [0_u32, 2, 5];
    let mu = [1.0, 2.5, 4.0];
    let weights = [0.5, 1.0, 2.0];
    let dispersion = 0.3;

    let expected = counts
        .iter()
        .copied()
        .zip(mu)
        .zip(weights)
        .map(|((count, mu), weight)| weight * nbinom_log_pmf(count, mu, dispersion).unwrap())
        .sum::<f64>();
    let actual = nbinom_log_likelihood_weighted(&counts, &mu, dispersion, Some(&weights)).unwrap();
    assert_relative_eq!(actual, expected, epsilon = 1e-12);
}

#[test]
fn matrix_log_likelihood_is_rowwise() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![0, 2, 5, 10, 12, 15]).unwrap();
    let mu = RowMajorMatrix::from_row_major(2, 3, vec![1.0, 2.5, 4.0, 9.0, 12.0, 16.0]).unwrap();
    let dispersions = [0.3, 0.1];

    let actual = nbinom_log_likelihood_matrix(&counts, &mu, &dispersions, None).unwrap();
    assert_relative_eq!(
        actual[0],
        nbinom_log_likelihood(&[0, 2, 5], &[1.0, 2.5, 4.0], 0.3).unwrap(),
        epsilon = 1e-12
    );
    assert_relative_eq!(
        actual[1],
        nbinom_log_likelihood(&[10, 12, 15], &[9.0, 12.0, 16.0], 0.1).unwrap(),
        epsilon = 1e-12
    );
}

#[test]
fn negative_twice_log_likelihood_matches_deseq2_deviance_convention() {
    let counts = [1_u32, 3];
    let mu = [1.5, 2.5];
    let dispersion = 0.4;
    let log_like = nbinom_log_likelihood(&counts, &mu, dispersion).unwrap();
    let deviance = nbinom_negative_twice_log_likelihood(&counts, &mu, dispersion).unwrap();
    assert_relative_eq!(deviance, -2.0 * log_like, epsilon = 1e-12);
}

#[test]
fn math_distribution_helpers_wrap_nb_likelihood_primitives() {
    let counts = [1_u32, 3];
    let mu = [1.5, 2.5];
    let weights = [0.25, 1.5];
    let dispersion = 0.4;
    let helpers = negative_binomial_helpers();

    assert_relative_eq!(
        negative_binomial_log_pmf(counts[0], mu[0], dispersion).unwrap(),
        nbinom_log_pmf(counts[0], mu[0], dispersion).unwrap(),
        epsilon = 1e-12
    );
    assert_relative_eq!(
        helpers.log_likelihood(&counts, &mu, dispersion).unwrap(),
        nbinom_log_likelihood(&counts, &mu, dispersion).unwrap(),
        epsilon = 1e-12
    );
    assert_relative_eq!(
        helpers
            .log_likelihood_weighted(&counts, &mu, dispersion, Some(&weights))
            .unwrap(),
        nbinom_log_likelihood_weighted(&counts, &mu, dispersion, Some(&weights)).unwrap(),
        epsilon = 1e-12
    );
    assert_relative_eq!(
        negative_binomial_negative_twice_log_likelihood(&counts, &mu, dispersion).unwrap(),
        nbinom_negative_twice_log_likelihood(&counts, &mu, dispersion).unwrap(),
        epsilon = 1e-12
    );
}

#[test]
fn nb_log_pmf_approaches_poisson_for_small_dispersion() {
    let count = 4_u32;
    let mu = 3.0_f64;
    let nb = nbinom_log_pmf(count, mu, 1e-8).unwrap();
    let poisson = -mu + f64::from(count) * mu.ln() - ln_gamma(f64::from(count) + 1.0);
    assert_relative_eq!(nb, poisson, epsilon = 1e-5);
}

#[test]
fn nb_log_pmf_zero_count_stays_finite_for_tiny_mu_dispersion() {
    let actual = nbinom_log_pmf(0, 1.0e-200, 1.0e-200).unwrap();

    assert!(actual.is_finite());
    assert_relative_eq!(actual, 0.0, epsilon = 1e-12);
}

#[test]
fn nb_log_pmf_stays_finite_when_mu_dispersion_product_overflows() {
    let actual = nbinom_log_pmf(3, 1.0e200, 1.0e200).unwrap();
    let size = 1.0e-200_f64;
    let expected = ln_gamma(3.0 + size)
        - ln_gamma(size)
        - ln_gamma(4.0)
        - size * (200.0_f64 * 10.0_f64.ln() + 200.0_f64 * 10.0_f64.ln());

    assert!(actual.is_finite());
    assert_relative_eq!(actual, expected, epsilon = 1e-12);
}

#[test]
fn nb_log_pmf_keeps_extreme_zero_count_product_finite() {
    let actual = nbinom_log_pmf(0, f64::MAX, f64::MAX).unwrap();

    assert!(actual.is_finite());
    assert_relative_eq!(actual, 0.0, epsilon = 1e-12);
}

#[test]
fn weighted_log_likelihood_rejects_nonfinite_weighted_term() {
    let err = nbinom_log_likelihood_weighted(&[1], &[1.0], 0.3, Some(&[f64::MAX])).unwrap_err();

    assert!(err
        .to_string()
        .contains("negative-binomial weighted log-likelihood term"));
}

#[test]
fn nb_likelihood_validates_inputs() {
    assert!(nbinom_log_pmf(1, 0.0, 0.1).is_err());
    assert!(nbinom_log_pmf(1, 1.0, 0.0).is_err());
    assert!(nbinom_log_likelihood(&[1, 2], &[1.0], 0.1).is_err());
    assert!(nbinom_log_likelihood_weighted(&[1], &[1.0], 0.1, Some(&[-1.0])).is_err());
}
