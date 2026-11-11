//! Source-compatible primitive optimizer controls, shared by all frontends.
use super::{
    InteriorPointOptions, LocationShrinkageOptions, MixtureFitOptions, MixtureOptimizer,
    PolynomialExtrapolation, SquaremOptions,
};
use crate::{
    DeseqError,
    design::expression::{Data, Value},
};
mod sqp;

/// Translate ashr optimizer controls without executing R. Fixed-prior calls
/// still validate the optimizer name but do not inspect unused control values.
/// Controls use the standalone expression engine's primitive named-list format.
pub fn ashr_mixture_fit_options(
    method: &str,
    control: Option<&Value>,
    fixed_prior: bool,
) -> Result<MixtureFitOptions, DeseqError> {
    let control = if fixed_prior { None } else { control };
    let parse = || -> Result<MixtureFitOptions, String> {
        let optimizer = policy(method, control)?;
        let mut options = MixtureFitOptions {
            optimizer,
            ..LocationShrinkageOptions::default().fit_options
        };
        if optimizer.is_sqp() {
            sqp::controls(control, &mut options)?;
        }
        Ok(options)
    };
    parse().map_err(|reason| DeseqError::InvalidOptions { reason })
}
fn fields(v: &Value) -> Result<(&[String], &[Value]), String> {
    let Data::List(values) = &v.data else {
        return Err("optimizer control requires a named list".into());
    };
    let names = v
        .names
        .as_deref()
        .ok_or("optimizer control requires names")?;
    if names.len() != values.len()
        || names
            .iter()
            .enumerate()
            .any(|(i, n)| names[..i].contains(n))
    {
        return Err("duplicate or invalid optimizer control fields".into());
    }
    Ok((names, values))
}
fn field<'a>(v: &'a Value, name: &str) -> Option<&'a Value> {
    let Data::List(values) = &v.data else {
        return None;
    };
    v.names
        .as_ref()?
        .iter()
        .position(|n| n == name)
        .and_then(|i| values.get(i))
}
fn scalar(v: &Value) -> Result<f64, String> {
    match &v.data {
        Data::Number { values, .. } => match values.as_slice() {
            [v] => Ok(*v),
            _ => Err("posterior control requires one numeric value".into()),
        },
        _ => Err("posterior request requires numeric vectors".into()),
    }
}
fn number(v: &Value, name: &str, default: f64) -> Result<f64, String> {
    field(v, name).map(scalar).unwrap_or(Ok(default))
}
fn flag(v: &Value) -> Result<bool, String> {
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
fn count(value: f64, name: &str) -> Result<usize, String> {
    if !value.is_finite() || value < 0.0 || value.fract() != 0.0 || value > i32::MAX as f64 {
        return Err(format!("control '{name}' requires a nonnegative integer"));
    }
    Ok(value as usize)
}
fn policy(method: &str, control: Option<&Value>) -> Result<MixtureOptimizer, String> {
    if method == "mixSQP" {
        return Ok(MixtureOptimizer::SequentialQuadratic);
    }
    if method == "mixIP" {
        let mut options = InteriorPointOptions::default();
        if let Some(c) = control.filter(|v| !matches!(&v.data,Data::List(v) if v.is_empty())) {
            let (names, _) = fields(c)?;
            for name in names {
                if name != "rtol" {
                    return Err(format!("unknown interior-point control '{name}'"));
                }
            }
            options.relative_tolerance = number(c, "rtol", 1e-6)?;
        }
        return Ok(MixtureOptimizer::InteriorPoint(options));
    }
    let mut options = SquaremOptions::default();
    if let Some(control) =
        control.filter(|v| !matches!(&v.data,Data::List(values) if values.is_empty()))
    {
        let (names, _) = fields(control)?;
        for name in names {
            match name.as_str() {
                "tol" => options.tolerance = Some(number(control, name, 1e-7)?),
                "maxiter" => options.max_evaluations = count(number(control, name, 5000.0)?, name)?,
                "method" => match &field(control, name).unwrap().data {
                    Data::Text(v) => {
                        if method == "cxxMixSquarem" {
                            return Err(
                                "compiled-style EM requires a numeric method control".into()
                            );
                        }
                        options.polynomial_method = if v.as_slice() == [Some("mpe".into())] {
                            PolynomialExtrapolation::MinimalPolynomial
                        } else {
                            PolynomialExtrapolation::ReducedRank
                        };
                    }
                    _ => {
                        let v = number(control, name, 3.0)?;
                        options.method = if [1.0, 2.0, 3.0].contains(&v) {
                            v as u8
                        } else {
                            3
                        };
                    }
                },
                "step.min0" => options.minimum_step = number(control, name, 1.0)?,
                "step.max0" => options.maximum_step = number(control, name, 1.0)?,
                "mstep" => options.step_multiplier = number(control, name, 4.0)?,
                "objfn.inc" => options.objective_increase = number(control, name, 1.0)?,
                "square" => options.square = flag(field(control, name).unwrap())?,
                "trace" => {
                    flag(field(control, name).unwrap())?;
                }
                "kr" => {
                    number(control, name, 1.0)?;
                }
                "K" => options.polynomial_order = count(number(control, name, 1.0)?, name)?,
                _ => return Err(format!("unknown squared-EM control '{name}'")),
            }
        }
    }
    Ok(match method {
        "mixEM" => MixtureOptimizer::Em(options),
        "w_mixEM" => MixtureOptimizer::WeightedEm(options),
        "mixVBEM" => MixtureOptimizer::VariationalEm(options),
        "cxxMixSquarem" => MixtureOptimizer::CxxEm(options),
        _ => return Err("unknown native ashr optimizer".into()),
    })
}
