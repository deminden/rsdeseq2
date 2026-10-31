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

fn normal_75_quantile() -> f64 {
    0.674_489_750_196_081_7
}

fn normal_matched_variance(abs_quantile: f64) -> f64 {
    (abs_quantile / normal_75_quantile()).powi(2)
}

#[path = "glm_irls/basic.rs"]
mod basic;
#[path = "glm_irls/expanded.rs"]
mod expanded;
#[path = "glm_irls/irls.rs"]
mod irls;
#[path = "glm_irls/prior.rs"]
mod prior;
#[path = "glm_irls/prior_refits.rs"]
mod prior_refits;
