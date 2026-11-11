//! Whole-batch location shrinkage marshalling; no numerical R callbacks.
use super::posterior::{numbers, scalar};
use super::{Data, Value, get, named, text};
use rsdeseq2::shrinkage::{
    LocationLikelihood, LocationShrinkage, LocationShrinkageOptions, MixtureComponents,
    MixtureGrid, MixtureGridOptions, MixturePriorPenalty, shrink_location_batch,
};
fn flag(value: &Value) -> Result<bool, String> {
    match &value.data {
        Data::Logical(v) => match v.as_slice() {
            [Some(v)] => Ok(*v),
            _ => Err("location control requires one nonmissing logical".into()),
        },
        _ => Err("location control requires logical storage".into()),
    }
}
fn optional(value: &Value) -> Result<Option<Vec<f64>>, String> {
    if matches!(value.data, Data::Null) {
        Ok(None)
    } else {
        Ok(Some(numbers(value)?.to_vec()))
    }
}
fn characters(values: Vec<String>) -> Value {
    Value::new(Data::Text(values.into_iter().map(Some).collect()))
}
pub(super) fn supplied_with_policy(
    value: &Value,
    penalty: &MixturePriorPenalty,
    variational: bool,
) -> Result<Option<MixtureGrid>, String> {
    if matches!(value.data, Data::Null) {
        return Ok(None);
    }
    let weights = numbers(get(value, "weights")?)?.to_vec();
    let values = |name| numbers(get(value, name)?).map(<[f64]>::to_vec);
    let components = match text(get(value, "component_class")?)? {
        "normalmix" => MixtureComponents::Normal {
            means: values("means")?,
            sds: values("sds")?,
        },
        "unimix" => MixtureComponents::Uniform {
            lower: values("lower")?,
            upper: values("upper")?,
        },
        "tnormalmix" => MixtureComponents::TruncatedNormal {
            means: values("means")?,
            sds: values("sds")?,
            lower: values("lower")?,
            upper: values("upper")?,
        },
        _ => return Err("unknown supplied-prior class".into()),
    };
    (if variational {
        MixtureGrid::from_supplied_variational(components, weights, penalty)
    } else {
        MixtureGrid::from_supplied(components, weights, penalty)
    })
    .map(Some)
    .map_err(|e| e.to_string())
}
pub(super) fn mixture(components: MixtureComponents, weights: Vec<f64>) -> Value {
    let (class, means, sds, lower, upper) = match components {
        MixtureComponents::Normal { means, sds } => {
            ("normalmix", Some(means), Some(sds), None, None)
        }
        MixtureComponents::Uniform { lower, upper } => {
            ("unimix", None, None, Some(lower), Some(upper))
        }
        MixtureComponents::TruncatedNormal {
            means,
            sds,
            lower,
            upper,
        } => (
            "tnormalmix",
            Some(means),
            Some(sds),
            Some(lower),
            Some(upper),
        ),
    };
    let optional = |v: Option<Vec<f64>>| {
        v.map(Value::numbers)
            .unwrap_or_else(|| Value::new(Data::Null))
    };
    named(
        &[
            "component_class",
            "weights",
            "means",
            "sds",
            "lower",
            "upper",
        ],
        vec![
            characters(vec![class.into()]),
            Value::numbers(weights),
            optional(means),
            optional(sds),
            optional(lower),
            optional(upper),
        ],
    )
}
fn encode(r: LocationShrinkage) -> Value {
    let fit = r
        .fit
        .map(|f| {
            named(
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
                    Value::numbers(f.weights),
                    Value::numbers(vec![f.log_likelihood]),
                    Value::numbers(vec![f.penalized_log_likelihood]),
                    Value::numbers(vec![f.dual_residual]),
                    Value::numbers(vec![f.optimizer_weight_sum]),
                    Value::numbers(vec![f.function_evaluations as f64]),
                    Value::new(Data::Logical(vec![Some(f.converged)])),
                    characters(vec![f.message]),
                ],
            )
        })
        .unwrap_or_else(|| Value::new(Data::Null));
    named(
        &[
            "fitted_prior",
            "fit",
            "posterior_mean",
            "posterior_sd",
            "negative_probability",
            "positive_probability",
            "local_false_discovery_rate",
            "local_false_sign_rate",
            "svalues",
            "qvalues",
            "cdf",
            "log_likelihood",
            "warnings",
        ],
        vec![
            mixture(r.fitted_grid.components, r.fitted_weights),
            fit,
            Value::numbers(r.posterior_mean),
            Value::numbers(r.posterior_sd),
            Value::numbers(r.negative_probability),
            Value::numbers(r.positive_probability),
            Value::numbers(r.local_false_discovery_rate),
            Value::numbers(r.local_false_sign_rate),
            Value::numbers(r.svalues),
            Value::numbers(r.qvalues),
            Value::numbers(r.cdf),
            Value::numbers(vec![r.log_likelihood]),
            characters(r.warnings),
        ],
    )
}
pub(super) fn evaluate(request: &Value) -> Result<Value, String> {
    let range = numbers(get(request, "range")?)?;
    if range.len() != 2 {
        return Err("location range requires two bounds".into());
    }
    let optimizer = super::ashr::field(request, "optimizer")
        .filter(|v| !matches!(v.data, Data::Null))
        .map(|v| super::ashr_optimizer::policy(text(v)?, super::ashr::field(request, "control")))
        .transpose()?
        .unwrap_or_default();
    let nullweight = get(request, "nullweight")?;
    let penalty = if let Some(p) =
        super::ashr::field(request, "penalty").filter(|v| !matches!(v.data, Data::Null))
    {
        if matches!(p.data, Data::Text(_)) && text(p)? == "unit" {
            MixturePriorPenalty::Unit
        } else {
            MixturePriorPenalty::Parameters(numbers(p)?.to_vec())
        }
    } else if matches!(nullweight.data, Data::Null) {
        MixturePriorPenalty::Uniform
    } else {
        MixturePriorPenalty::NullBiased(scalar(nullweight)?)
    };
    let df = get(request, "df")?;
    let mut options = LocationShrinkageOptions {
        grid: MixtureGridOptions {
            family: text(get(request, "family")?)?
                .parse()
                .map_err(|e: rsdeseq2::DeseqError| e.to_string())?,
            scales: optional(get(request, "scales")?)?,
            point_mass: flag(get(request, "point_mass")?)?,
            mode: scalar(get(request, "mode")?)?,
            range: (range[0], range[1]),
            penalty: penalty.clone(),
            concentrated_initialization: !optimizer.is_sqp(),
            variational_prior: optimizer.is_variational(),
            ..Default::default()
        },
        supplied_grid: supplied_with_policy(
            get(request, "supplied_prior")?,
            &penalty,
            optimizer.is_variational(),
        )?,
        fit_options: rsdeseq2::shrinkage::MixtureFitOptions {
            optimizer,
            ..LocationShrinkageOptions::default().fit_options
        },
        fixed_prior: flag(get(request, "fixed_prior")?)?,
        likelihood: if request
            .names
            .as_ref()
            .is_some_and(|n| n.iter().any(|n| n == "custom_likelihood"))
            && !matches!(get(request, "custom_likelihood")?.data, Data::Null)
        {
            LocationLikelihood::Restricted(std::sync::Arc::new(restricted_likelihood(get(
                request,
                "custom_likelihood",
            )?)?))
        } else if matches!(df.data, Data::Null) {
            LocationLikelihood::Normal
        } else {
            LocationLikelihood::StudentT(scalar(df)?)
        },
        alpha: scalar(get(request, "alpha")?)?,
        prune_threshold: scalar(get(request, "prune_threshold")?)?,
        observation_weights: optional(get(request, "observation_weights")?)?,
        thresholds: numbers(get(request, "thresholds")?)?.to_vec(),
    };
    if optimizer.is_sqp() {
        super::ashr_sqp::controls(
            super::ashr::field(request, "control"),
            &mut options.fit_options,
        )?;
    }
    shrink_location_batch(
        numbers(get(request, "estimates")?)?,
        numbers(get(request, "standard_errors")?)?,
        &options,
    )
    .map(encode)
    .map_err(|e| e.to_string())
}

fn restricted_function(value: &Value) -> Result<rsdeseq2::shrinkage::RestrictedFunction, String> {
    let (arguments, expression, constants) = function_parts(value)?;
    rsdeseq2::shrinkage::RestrictedFunction::new(arguments, expression, constants)
        .map_err(|e| e.to_string())
}
type FunctionParts = (
    Vec<String>,
    rsdeseq2::design::expression::Expression,
    Vec<(String, Value)>,
);
pub(super) fn function_parts(value: &Value) -> Result<FunctionParts, String> {
    let arguments = match &get(value, "arguments")?.data {
        Data::Text(v) => v
            .iter()
            .cloned()
            .collect::<Option<Vec<_>>>()
            .ok_or("likelihood argument names cannot be missing")?,
        _ => return Err("likelihood argument names require text storage".into()),
    };
    let constants = get(value, "constants")?;
    let constants = if matches!(&constants.data, Data::List(v) if v.is_empty()) {
        vec![]
    } else {
        let (names, values) = super::fields(constants)?;
        names.iter().cloned().zip(values.iter().cloned()).collect()
    };
    Ok((
        arguments,
        super::ast(get(value, "expression")?, 0)?,
        constants,
    ))
}
pub(super) fn restricted_likelihood(
    value: &Value,
) -> Result<rsdeseq2::shrinkage::RestrictedLikelihood, String> {
    let optional = |name| {
        let v = get(value, name)?;
        if matches!(v.data, Data::Null) {
            Ok(None)
        } else {
            restricted_function(v).map(Some)
        }
    };
    rsdeseq2::shrinkage::RestrictedLikelihood::new(
        text(get(value, "name")?)?.into(),
        restricted_function(get(value, "log_cdf")?)?,
        restricted_function(get(value, "log_density")?)?,
        optional("first")?,
        optional("second")?,
    )
    .map_err(|e| e.to_string())
}
