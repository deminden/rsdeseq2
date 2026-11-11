//! Optimizer selection and primitive controls; no R numerical execution.
use super::ashr::{field, flag, number};
use super::{Data, Value, fields};
use rsdeseq2::shrinkage::{
    InteriorPointOptions, MixtureOptimizer, PolynomialExtrapolation, SquaremOptions,
};
fn count(value: f64, name: &str) -> Result<usize, String> {
    if !value.is_finite() || value < 0.0 || value.fract() != 0.0 || value > i32::MAX as f64 {
        return Err(format!("control '{name}' requires a nonnegative integer"));
    }
    Ok(value as usize)
}
pub(super) fn policy(method: &str, control: Option<&Value>) -> Result<MixtureOptimizer, String> {
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

/// Direct primitive API used to verify optimizer diagnostics independently of
/// component/grid preparation and public result-table presentation.
pub(super) fn fit(request: &Value) -> Result<Value, String> {
    use super::posterior::numbers;
    use super::{get, named, text};
    use rsdeseq2::{
        matrix::{MatrixOrder, MatrixView},
        shrinkage::{MixtureFitOptions, fit_mixture_proportions},
    };
    let count = |name| super::integer(get(request, name)?);
    let matrix = MatrixView::new(
        count("observations")?,
        count("components")?,
        numbers(get(request, "log_likelihoods")?)?,
        MatrixOrder::RowMajor,
    )
    .map_err(|e| e.to_string())?;
    let optional = |name| -> Result<Option<&[f64]>, String> {
        let v = get(request, name)?;
        if matches!(v.data, Data::Null) {
            Ok(None)
        } else {
            numbers(v).map(Some)
        }
    };
    let f = fit_mixture_proportions(
        matrix,
        optional("observation_weights")?,
        numbers(get(request, "dirichlet")?)?,
        optional("initial_weights")?,
        MixtureFitOptions {
            optimizer: policy(text(get(request, "optimizer")?)?, field(request, "control"))?,
            ..Default::default()
        },
    )
    .map_err(|e| e.to_string())?;
    let opt = |v: Option<Vec<f64>>| {
        v.map(Value::numbers)
            .unwrap_or_else(|| Value::new(Data::Null))
    };
    Ok(named(
        &[
            "weights",
            "log_likelihood",
            "penalized_log_likelihood",
            "converged",
            "iterations",
            "evaluations",
            "residual",
            "weight_sum",
            "posterior_concentrations",
            "optimizer_objective",
            "message",
        ],
        vec![
            Value::numbers(f.weights),
            Value::numbers(vec![f.log_likelihood]),
            Value::numbers(vec![f.penalized_log_likelihood]),
            Value::new(Data::Logical(vec![Some(f.converged)])),
            opt(f.iterations.map(|v| vec![v as f64])),
            Value::numbers(vec![f.function_evaluations as f64]),
            Value::numbers(vec![f.dual_residual]),
            Value::numbers(vec![f.optimizer_weight_sum]),
            opt(f.posterior_concentrations),
            opt(f.optimizer_objective.map(|v| vec![v])),
            Value::new(Data::Text(vec![Some(f.message)])),
        ],
    ))
}
