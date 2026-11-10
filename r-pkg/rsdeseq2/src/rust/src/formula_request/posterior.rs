//! Primitive transport for standalone posterior kernels; no R callbacks.
use super::{Data, Value, get, named};
use rsdeseq2::shrinkage::{NormalMixture, svalues};

fn numbers(value: &Value) -> Result<&[f64], String> {
    match &value.data {
        Data::Number { values, .. } => Ok(values),
        _ => Err("posterior request requires numeric vectors".into()),
    }
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
