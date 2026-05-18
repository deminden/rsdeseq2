use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

#[test]
fn lrt_test_matches_chisq_upper_tail() {
    let full = toy_fit(vec![0.0, 0.0], vec![1.0, 1.0], vec![10.0], 1, 2);
    let reduced = toy_fit(vec![0.0], vec![1.0], vec![8.0], 1, 1);

    let lrt = lrt_test(&full, &reduced).unwrap();

    assert_eq!(lrt.degrees_of_freedom, 1);
    assert_eq!(lrt.reduced_converged, vec![true]);
    assert_relative_eq!(lrt.deviance[0].unwrap(), 4.0, epsilon = 1e-12);
    assert_relative_eq!(lrt.pvalue[0].unwrap(), 0.04550026389635853, epsilon = 1e-14);
}

#[test]
fn lrt_test_handles_missing_log_likelihoods() {
    let full = toy_fit(vec![0.0, 0.0], vec![1.0, 1.0], vec![f64::NAN], 1, 2);
    let reduced = toy_fit(vec![0.0], vec![1.0], vec![8.0], 1, 1);

    let lrt = lrt_test(&full, &reduced).unwrap();

    assert_eq!(lrt.deviance[0], None);
    assert_eq!(lrt.pvalue[0], None);
}

#[test]
fn lrt_test_validates_inputs() {
    let full = toy_fit(vec![0.0], vec![1.0], vec![10.0], 1, 1);
    let reduced = toy_fit(vec![0.0], vec![1.0], vec![8.0], 1, 1);
    assert!(lrt_test(&full, &reduced).is_err());

    let bad_reduced = toy_fit(vec![0.0, 0.0], vec![1.0, 1.0], vec![8.0, 9.0], 2, 1);
    assert!(lrt_test(&full, &bad_reduced).is_err());
}

fn toy_fit(
    beta: Vec<f64>,
    beta_se: Vec<f64>,
    log_like: Vec<f64>,
    n_genes: usize,
    n_coef: usize,
) -> NbinomGlmFit {
    let n_samples = 2;
    NbinomGlmFit {
        log_like,
        beta_converged: vec![true; n_genes],
        beta: RowMajorMatrix::from_row_major(n_genes, n_coef, beta).unwrap(),
        beta_se: RowMajorMatrix::from_row_major(n_genes, n_coef, beta_se).unwrap(),
        beta_covariance: None,
        mu: RowMajorMatrix::from_row_major(n_genes, n_samples, vec![1.0; n_genes * n_samples])
            .unwrap(),
        beta_iter: vec![1; n_genes],
        model_matrix: DesignMatrix::from_row_major(
            n_samples,
            n_coef,
            vec![1.0; n_samples * n_coef],
            None,
        )
        .unwrap(),
        n_terms: n_coef,
        hat_diagonal: RowMajorMatrix::from_row_major(
            n_genes,
            n_samples,
            vec![0.5; n_genes * n_samples],
        )
        .unwrap(),
    }
}
