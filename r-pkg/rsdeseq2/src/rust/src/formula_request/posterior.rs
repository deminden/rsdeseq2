//! Primitive transport for standalone posterior kernels; no R callbacks.
use super::{Data, Value, get, integer, named};
use rsdeseq2::math::{normal_log_cdf, normal_log_upper_tail};
use rsdeseq2::matrix::{MatrixOrder, MatrixView};
use rsdeseq2::shrinkage::{MixtureFitOptions, NormalMixture, fit_mixture_proportions, svalues};

fn numbers(value: &Value) -> Result<&[f64], String> {
    match &value.data {
        Data::Number { values, .. } => Ok(values),
        _ => Err("posterior request requires numeric vectors".into()),
    }
}
fn scalar(value: &Value) -> Result<f64, String> {
    match numbers(value)? {
        [value] => Ok(*value),
        _ => Err("posterior control requires one numeric value".into()),
    }
}

pub(super) fn normal_log(request: &Value) -> Result<Value, String> {
    let arguments = numbers(get(request, "arguments")?)?;
    Ok(named(
        &["log_cdf", "log_upper"],
        vec![
            Value::numbers(arguments.iter().map(|v| normal_log_cdf(*v)).collect()),
            Value::numbers(
                arguments
                    .iter()
                    .map(|v| normal_log_upper_tail(*v))
                    .collect(),
            ),
        ],
    ))
}

pub(super) fn fit_mixture(request: &Value) -> Result<Value, String> {
    let controls = get(request, "options")?;
    let options = MixtureFitOptions {
        max_iterations: integer(get(controls, "max_iterations")?)?,
        gradient_tolerance: scalar(get(controls, "gradient_tolerance")?)?,
        em_iterations: integer(get(controls, "em_iterations")?)?,
        active_iterations: integer(get(controls, "active_iterations")?)?,
        zero_threshold: scalar(get(controls, "zero_threshold")?)?,
        likelihood_offset: scalar(get(controls, "likelihood_offset")?)?,
    };
    let initial = get(request, "initial_weights")?;
    let initial = if matches!(initial.data, Data::Null) {
        None
    } else {
        Some(numbers(initial)?)
    };
    let observations = integer(get(request, "observations")?)?;
    let components = integer(get(request, "components")?)?;
    let likelihoods = MatrixView::new(
        observations,
        components,
        numbers(get(request, "log_likelihoods")?)?,
        MatrixOrder::RowMajor,
    )
    .map_err(|e| e.to_string())?;
    let fit = fit_mixture_proportions(
        likelihoods,
        Some(numbers(get(request, "observation_weights")?)?),
        numbers(get(request, "dirichlet")?)?,
        initial,
        options,
    )
    .map_err(|e| e.to_string())?;
    Ok(named(
        &[
            "weights",
            "log_likelihood",
            "penalized_log_likelihood",
            "dual_residual",
            "optimizer_weight_sum",
            "function_evaluations",
            "converged",
            "message",
        ],
        vec![
            Value::numbers(fit.weights),
            Value::numbers(vec![fit.log_likelihood]),
            Value::numbers(vec![fit.penalized_log_likelihood]),
            Value::numbers(vec![fit.dual_residual]),
            Value::numbers(vec![fit.optimizer_weight_sum]),
            Value::numbers(vec![fit.function_evaluations as f64]),
            Value::new(Data::Logical(vec![Some(fit.converged)])),
            Value::new(Data::Text(vec![Some(fit.message)])),
        ],
    ))
}

pub(super) fn normal_mixture(request: &Value) -> Result<Value, String> {
    let prior = NormalMixture::new(
        numbers(get(request, "weights")?)?.to_vec(),
        numbers(get(request, "means")?)?.to_vec(),
        numbers(get(request, "sds")?)?.to_vec(),
    )
    .map_err(|e| e.to_string())?;
    let estimates = numbers(get(request, "estimates")?)?;
    let errors = numbers(get(request, "standard_errors")?)?;
    if estimates.len() != errors.len() {
        return Err("posterior estimates and standard errors must have matching lengths".into());
    }
    let thresholds = numbers(get(request, "thresholds")?)?;
    if thresholds.iter().any(|v| v.is_nan()) {
        return Err("posterior thresholds cannot be missing".into());
    }
    let mut means = Vec::with_capacity(estimates.len());
    let mut sds = Vec::with_capacity(estimates.len());
    let mut weights = vec![];
    let mut cdf = vec![];
    let mut upper = vec![];
    for (&estimate, &se) in estimates.iter().zip(errors) {
        let posterior = prior.posterior(estimate, se).map_err(|e| e.to_string())?;
        means.push(posterior.mean());
        sds.push(posterior.standard_deviation());
        weights.extend_from_slice(posterior.weights());
        for &threshold in thresholds {
            cdf.push(posterior.cdf(threshold));
            upper.push(posterior.survival(threshold));
        }
    }
    Ok(named(
        &[
            "posterior_mean",
            "posterior_sd",
            "posterior_weights",
            "cdf",
            "upper_tail",
        ],
        vec![
            Value::numbers(means),
            Value::numbers(sds),
            Value::numbers(weights),
            Value::numbers(cdf),
            Value::numbers(upper),
        ],
    ))
}

pub(super) fn ordered_svalues(request: &Value) -> Result<Value, String> {
    let rates = numbers(get(request, "rates")?)?
        .iter()
        .map(|v| (!v.is_nan()).then_some(*v))
        .collect::<Vec<_>>();
    let values = svalues(&rates).map_err(|e| e.to_string())?;
    Ok(Value::numbers(
        values
            .into_iter()
            .map(|v| v.unwrap_or(f64::from_bits(0x7ff0_0000_0000_07a2)))
            .collect(),
    ))
}
