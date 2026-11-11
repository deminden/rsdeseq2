//! Direct primitive mixture transport; control policy lives in standalone core.
use super::{Data, Value};
/// Direct primitive API used to verify optimizer diagnostics independently of
/// component/grid preparation and public result-table presentation.
pub(super) fn fit(request: &Value) -> Result<Value, String> {
    use super::posterior::numbers;
    use super::{get, named, text};
    use rsdeseq2::{
        matrix::{MatrixOrder, MatrixView},
        shrinkage::fit_mixture_proportions,
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
    let options = rsdeseq2::shrinkage::ashr_mixture_fit_options(
        text(get(request, "optimizer")?)?,
        super::ashr::field(request, "control"),
        false,
    )
    .map_err(|e| e.to_string())?;
    let f = fit_mixture_proportions(
        matrix,
        optional("observation_weights")?,
        numbers(get(request, "dirichlet")?)?,
        optional("initial_weights")?,
        options,
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
