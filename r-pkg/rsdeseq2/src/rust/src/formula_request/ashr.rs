//! Source-method argument translation. All statistical operations stay in core.
use super::location_shrinkage::{mixture, restricted_likelihood, supplied_with_policy};
use super::posterior::{numbers, scalar};
use super::{Data, Value, fields, get, named, text};
use rsdeseq2::shrinkage::{
    AdaptiveMethod, AdaptiveOutputSelection, AdaptiveShrinkageQuery, LocationLikelihood,
    LocationShrinkageOptions, MixtureGridOptions, MixturePriorPenalty,
    adaptive_shrinkage_with_output,
};
pub(super) fn field<'a>(v: &'a Value, name: &str) -> Option<&'a Value> {
    let Data::List(values) = &v.data else {
        return None;
    };
    v.names
        .as_ref()?
        .iter()
        .position(|n| n == name)
        .and_then(|i| values.get(i))
}
pub(super) fn number(v: &Value, name: &str, default: f64) -> Result<f64, String> {
    field(v, name).map(scalar).unwrap_or(Ok(default))
}
pub(super) fn flag(v: &Value) -> Result<bool, String> {
    match &v.data {
        Data::Logical(v) => match v.as_slice() {
            [Some(v)] => Ok(*v),
            _ => Err("logical option must be one nonmissing value".into()),
        },
        Data::Number { values, .. } => match values.as_slice() {
            [v] if !v.is_nan() => Ok(*v != 0.0),
            _ => Err("logical option must be one nonmissing value".into()),
        },
        _ => Err("logical option requires numeric or logical storage".into()),
    }
}
fn optional_numbers(v: Option<&Value>) -> Result<Option<Vec<f64>>, String> {
    match v {
        None => Ok(None),
        Some(v) if matches!(v.data, Data::Null) => Ok(None),
        Some(v) => numbers(v).map(|v| Some(v.to_vec())),
    }
}
pub(super) fn choice(v: Option<&Value>, choices: &[&str]) -> Result<String, String> {
    let Some(v) = v else {
        return Ok(choices[0].into());
    };
    if matches!(v.data, Data::Null) {
        return Ok(choices[0].into());
    }
    let Data::Text(values) = &v.data else {
        return Err("method option requires text".into());
    };
    if values
        .iter()
        .map(|v| v.as_deref())
        .eq(choices.iter().map(|v| Some(*v)))
    {
        return Ok(choices[0].into());
    }
    let [Some(value)] = values.as_slice() else {
        return Err("method option requires one nonmissing name".into());
    };
    if choices.contains(&value.as_str()) {
        return Ok(value.clone());
    }
    let matches = choices
        .iter()
        .filter(|name| name.starts_with(value))
        .collect::<Vec<_>>();
    if matches.len() == 1 {
        Ok((*matches[0]).into())
    } else {
        Err(format!("unknown or ambiguous method option '{value}'"))
    }
}
fn controls(value: Option<&Value>, options: &mut LocationShrinkageOptions) -> Result<(), String> {
    let Some(value) = value else {
        return Ok(());
    };
    if matches!(&value.data,Data::List(v) if v.is_empty()) {
        return Ok(());
    }
    let (names, values) = fields(value)?;
    for (name, value) in names.iter().zip(values) {
        let count = || -> Result<usize, String> {
            let x = scalar(value)?;
            if !x.is_finite() || x < 0.0 || x.fract() != 0.0 || x > i32::MAX as f64 {
                return Err(format!("control '{name}' requires a nonnegative integer"));
            }
            Ok(x as usize)
        };
        match name.as_str() {
            "verbose" => {
                flag(value)?;
            }
            "eps" => options.fit_options.likelihood_offset = scalar(value)?,
            "numiter.em" => options.fit_options.em_iterations = count()?,
            "maxiter.sqp" => options.fit_options.max_iterations = count()?,
            "maxiter.activeset" => options.fit_options.active_iterations = count()?,
            "convtol.sqp" => options.fit_options.gradient_tolerance = scalar(value)?,
            "zero.threshold.solution" => options.fit_options.zero_threshold = scalar(value)?,
            _ => {
                return Err(format!(
                    "native ashr control '{name}' is not implemented yet"
                ));
            }
        }
    }
    Ok(())
}
pub(super) fn evaluate(request: &Value) -> Result<Value, String> {
    let args = get(request, "arguments")?;
    let (names, _) = fields(args)?;
    let allowed = [
        "betahat",
        "sebetahat",
        "method",
        "mixcompdist",
        "optmethod",
        "df",
        "nullweight",
        "pointmass",
        "prior",
        "mixsd",
        "gridmult",
        "outputlevel",
        "g",
        "fixg",
        "mode",
        "alpha",
        "grange",
        "control",
        "lik",
        "weights",
        "pi_thresh",
    ];
    for name in names {
        if !allowed.contains(&name.as_str()) {
            return Err(format!("unused ashr argument '{name}'"));
        }
    }
    let family = choice(
        field(args, "mixcompdist"),
        &[
            "uniform",
            "halfuniform",
            "normal",
            "+uniform",
            "-uniform",
            "halfnormal",
        ],
    )?;
    let optmethod = choice(
        field(args, "optmethod"),
        &[
            "mixSQP",
            "mixIP",
            "cxxMixSquarem",
            "mixEM",
            "mixVBEM",
            "w_mixEM",
        ],
    )?;
    let optimizer = super::ashr_optimizer::policy(&optmethod, field(args, "control"))?;
    let method = if field(args, "method").is_some() {
        Some(
            if choice(field(args, "method"), &["fdr", "shrink"])? == "shrink" {
                AdaptiveMethod::Shrink
            } else {
                AdaptiveMethod::Fdr
            },
        )
    } else {
        None
    };
    let nullweight = number(args, "nullweight", 10.0)?;
    let penalty = if let Some(v) = field(args, "prior") {
        Some(match &v.data {
            Data::Number { values, .. } => MixturePriorPenalty::Parameters(values.clone()),
            _ => match choice(Some(v), &["nullbiased", "uniform", "unit"])?.as_str() {
                "nullbiased" => MixturePriorPenalty::NullBiased(nullweight),
                "uniform" => MixturePriorPenalty::Uniform,
                _ if optimizer.is_variational() => MixturePriorPenalty::Unit,
                _ => return Err("unit prior requires a variational optimizer".into()),
            },
        })
    } else {
        None
    };
    let point_mass = field(args, "pointmass").map(flag).transpose()?;
    let range =
        optional_numbers(field(args, "grange"))?.unwrap_or(vec![f64::NEG_INFINITY, f64::INFINITY]);
    if range.is_empty() || range.iter().any(|v| v.is_nan()) {
        return Err("grange requires nonmissing numeric bounds".into());
    }
    let lower = range.iter().copied().fold(f64::INFINITY, f64::min);
    let upper = range.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if field(args, "mode").is_some() && field(args, "g").is_some() {
        return Err("specify either mode or g, not both".into());
    }
    let mut options = LocationShrinkageOptions {
        grid: MixtureGridOptions {
            family: family
                .parse()
                .map_err(|e: rsdeseq2::DeseqError| e.to_string())?,
            scales: optional_numbers(field(args, "mixsd"))?,
            grid_multiplier: number(args, "gridmult", std::f64::consts::SQRT_2)?,
            mode: number(args, "mode", 0.0)?,
            range: (lower, upper),
            concentrated_initialization: !optimizer.is_sqp(),
            variational_prior: optimizer.is_variational(),
            ..Default::default()
        }
        .with_ashr_method(method, point_mass, penalty, nullweight)
        .map_err(|e| e.to_string())?,
        fixed_prior: field(args, "fixg").map(flag).transpose()?.unwrap_or(false),
        alpha: number(args, "alpha", 0.0)?,
        prune_threshold: number(args, "pi_thresh", 1e-10)?,
        observation_weights: optional_numbers(field(args, "weights"))?,
        ..Default::default()
    };
    if let Some(g) = field(args, "g") {
        options.supplied_grid =
            supplied_with_policy(g, &options.grid.penalty, optimizer.is_variational())?;
    }
    let df = optional_numbers(field(args, "df"))?;
    let df = match df.as_deref() {
        None | Some([]) => None,
        Some([v]) if v.is_infinite() => None,
        Some([v]) => Some(*v),
        _ => return Err("only one df value can be specified".into()),
    };
    if family == "normal" && df.is_some() {
        return Err("normal mixture for Student-t likelihood is not implemented".into());
    }
    let mut input_warnings = vec![];
    options.likelihood =
        if let Some(lik) = field(args, "lik").filter(|v| !matches!(v.data, Data::Null)) {
            match text(get(lik, "type")?)? {
                "normal" => LocationLikelihood::Normal,
                "t" => LocationLikelihood::StudentT(scalar(get(lik, "df")?)?),
                "restricted" => {
                    let f = restricted_likelihood(get(lik, "model")?)?;
                    if family == "normal" && !matches!(f.name(), "normal" | "normalmix") {
                        return Err("normal mixtures require a normal likelihood".into());
                    }
                    if df.is_some() && !matches!(f.name(), "normal" | "t") {
                        input_warnings
                            .push("Input df is ignored for this likelihood function".to_string());
                    }
                    LocationLikelihood::Restricted(std::sync::Arc::new(f))
                }
                kind => {
                    return Err(format!(
                        "native ashr likelihood '{kind}' is not implemented yet"
                    ));
                }
            }
        } else {
            df.map(LocationLikelihood::StudentT)
                .unwrap_or(LocationLikelihood::Normal)
        };
    options.fit_options.optimizer = optimizer;
    if optimizer.is_sqp() {
        controls(field(args, "control"), &mut options)?;
    }
    let output = match field(args, "outputlevel").map(|v| &v.data) {
        None => AdaptiveOutputSelection::level(2.0),
        Some(Data::Number { values, .. }) => match values.as_slice() {
            [level] => AdaptiveOutputSelection::level(*level),
            _ => return Err("outputlevel requires one numeric level or a text subset".into()),
        },
        Some(Data::Text(names)) => Ok(AdaptiveOutputSelection::names(
            names.iter().flatten().map(String::as_str),
        )),
        Some(Data::Null | Data::Logical(_)) => Ok(AdaptiveOutputSelection::default()),
        _ => return Err("outputlevel requires one numeric level or a text subset".into()),
    }
    .map_err(|e| e.to_string())?;
    let result = adaptive_shrinkage_with_output(
        numbers(get(args, "betahat")?)?,
        numbers(get(args, "sebetahat")?)?,
        &options,
        AdaptiveShrinkageQuery {
            estimates: numbers(get(request, "query_estimates")?)?,
            standard_errors: numbers(get(request, "query_standard_errors")?)?,
            lfc_threshold: scalar(get(request, "lfc_threshold")?)?,
            svalues: flag(get(request, "svalue")?)?,
        },
        output,
    )
    .map_err(|e| e.to_string())?;
    input_warnings.extend(result.fit.warnings);
    let fitted = if result.fitted_prior {
        mixture(result.fit.fitted_grid.components, result.fit.fitted_weights)
    } else {
        Value::new(Data::Null)
    };
    let optional = |values: Option<Vec<f64>>| {
        values
            .map(Value::numbers)
            .unwrap_or_else(|| Value::new(Data::Null))
    };
    Ok(named(
        &[
            "posterior_mean",
            "posterior_sd",
            "svalue",
            "fitted_prior",
            "warnings",
        ],
        vec![
            optional(result.posterior_mean),
            optional(result.posterior_sd),
            optional(result.svalues),
            fitted,
            Value::new(Data::Text(input_warnings.into_iter().map(Some).collect())),
        ],
    ))
}
