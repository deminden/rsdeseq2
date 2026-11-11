use rsdeseq2::{
    matrix::{MatrixOrder, MatrixView},
    shrinkage::{MixtureFitOptions, fit_mixture_proportions},
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;
#[derive(Deserialize)]
pub struct Input {
    mixture_fit_cases: Vec<Case>,
}
#[derive(Deserialize)]
struct Case {
    name: String,
    observations: usize,
    components: usize,
    log_likelihoods: Vec<String>,
    observation_weights: Vec<String>,
    dirichlet: Vec<String>,
    initial_weights: Option<Vec<String>>,
    options: Controls,
}
#[derive(Deserialize)]
struct Controls {
    max_iterations: usize,
    gradient_tolerance: String,
    em_iterations: usize,
    active_iterations: usize,
    zero_threshold: String,
    likelihood_offset: String,
}
use super::{parse_numeric_input as number, parse_numeric_inputs as numbers};
pub fn evaluate(input: Input) -> Result<Value, Box<dyn Error>> {
    let mut cases = vec![];
    for c in input.mixture_fit_cases {
        let data = numbers(&c.log_likelihoods)?;
        let a = numbers(&c.observation_weights)?;
        let prior = numbers(&c.dirichlet)?;
        let initial = c.initial_weights.map(|v| numbers(&v)).transpose()?;
        let fit = fit_mixture_proportions(
            MatrixView::new(c.observations, c.components, &data, MatrixOrder::RowMajor)?,
            Some(&a),
            &prior,
            initial.as_deref(),
            MixtureFitOptions {
                positive_coordinate_convergence: false,
                optimizer: Default::default(),
                max_iterations: c.options.max_iterations,
                gradient_tolerance: number(&c.options.gradient_tolerance)?,
                em_iterations: c.options.em_iterations,
                active_iterations: c.options.active_iterations,
                zero_threshold: number(&c.options.zero_threshold)?,
                likelihood_offset: number(&c.options.likelihood_offset)?,
                ..Default::default()
            },
        )?;
        cases.push(json!({"name":c.name,"weights":fit.weights.into_iter().map(super::bits).collect::<Vec<_>>(),
            "log_likelihood":super::bits(fit.log_likelihood),"penalized_log_likelihood":super::bits(fit.penalized_log_likelihood),
            "dual_residual":super::bits(fit.dual_residual),"optimizer_weight_sum":super::bits(fit.optimizer_weight_sum),
            "function_evaluations":super::bits(fit.function_evaluations as f64),"converged":fit.converged,"message":fit.message}));
    }
    Ok(json!({"cases":cases}))
}
