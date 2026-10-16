//! Integration coverage grouped by workflow; shared fixtures remain in this test crate.

use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

fn toy_fit(beta: Vec<f64>, beta_se: Vec<f64>, beta_converged: Vec<bool>) -> NbinomGlmFit {
    let n_genes = beta_converged.len();
    let n_samples = 2;
    NbinomGlmFit {
        log_like: vec![0.0; n_genes],
        beta_converged,
        beta: RowMajorMatrix::from_row_major(n_genes, 1, beta).unwrap(),
        beta_se: RowMajorMatrix::from_row_major(n_genes, 1, beta_se).unwrap(),
        beta_optim_start: RowMajorMatrix::from_elem(n_genes, 1, f64::NAN).unwrap(),
        beta_covariance: None,
        mu: RowMajorMatrix::from_row_major(n_genes, n_samples, vec![1.0; n_genes * n_samples])
            .unwrap(),
        beta_iter: vec![1; n_genes],
        beta_optim_iter: vec![f64::NAN; n_genes],
        beta_optim_start_objective: vec![f64::NAN; n_genes],
        beta_optim_objective: vec![f64::NAN; n_genes],
        beta_optim_gradient_norm: vec![f64::NAN; n_genes],
        model_matrix: DesignMatrix::from_row_major(
            n_samples,
            1,
            vec![1.0, 1.0],
            Some(vec!["Intercept".to_string()]),
        )
        .unwrap(),
        n_terms: 1,
        hat_diagonal: RowMajorMatrix::from_row_major(
            n_genes,
            n_samples,
            vec![0.5; n_genes * n_samples],
        )
        .unwrap(),
    }
}

#[path = "results/tables.rs"]
mod tables;

#[path = "results/expanded_results.rs"]
mod expanded_results;

#[path = "results/replacement.rs"]
mod replacement;

#[path = "results/factor_models.rs"]
mod factor_models;

#[path = "results/additive_models.rs"]
mod additive_models;

#[path = "results/interactions.rs"]
mod interactions;

#[path = "results/formula_models.rs"]
mod formula_models;

#[path = "results/weighted_formulas.rs"]
mod weighted_formulas;

#[path = "results/model_frame.rs"]
mod model_frame;

#[path = "results/lrt.rs"]
mod lrt;

#[path = "results/cooks.rs"]
mod cooks;
