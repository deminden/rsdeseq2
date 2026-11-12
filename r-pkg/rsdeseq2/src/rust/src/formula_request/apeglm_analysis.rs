//! Optional posterior argument/output transport. No frontend arithmetic.
use super::ashr::{choice, field, number};
use super::posterior::numbers;
use super::{Data, Value, get, named, text};
use rsdeseq2::{
    matrix::RowMajorMatrix,
    shrinkage::{
        PosteriorAnalysis, PosteriorAnalysisOptions, PosteriorIntervalKind,
        PosteriorIntervalOptions, RegressionLikelihood, RestrictedRegressionLikelihood,
    },
};

pub(super) fn likelihood(model: &Value) -> Result<RegressionLikelihood, String> {
    match text(get(model, "type")?)? {
        "negative_binomial" => Ok(RegressionLikelihood::NegativeBinomial),
        "restricted" => {
            let (a, e, c) = super::location_shrinkage::function_parts(get(model, "model")?)?;
            Ok(RegressionLikelihood::Restricted(std::sync::Arc::new(
                RestrictedRegressionLikelihood::new(a, e, c).map_err(|e| e.to_string())?,
            )))
        }
        _ => Err("unsupported regression likelihood descriptor".into()),
    }
}
fn present(v: Option<&Value>) -> Option<&Value> {
    v.filter(|v| !matches!(v.data, Data::Null))
}
fn count(args: &Value, name: &str, default: usize) -> Result<usize, String> {
    let n = number(args, name, default as f64)?;
    if !n.is_finite() || n < 0.0 || n.fract() != 0.0 || n >= usize::MAX as f64 {
        return Err(format!("{name} requires a nonnegative integer"));
    }
    Ok(n as usize)
}
pub(super) fn options(args: &Value) -> Result<PosteriorAnalysisOptions, String> {
    let kind = match choice(
        field(args, "interval.type"),
        &["laplace", "HPD", "credible"],
    )?
    .as_str()
    {
        "HPD" => PosteriorIntervalKind::HighestDensity,
        "credible" => PosteriorIntervalKind::EqualTail,
        _ => PosteriorIntervalKind::Laplace,
    };
    let explicit = [
        "interval.type",
        "interval.level",
        "ngrid",
        "nsd",
        "ngrid.nuis",
        "nsd.nuis",
    ]
    .iter()
    .any(|n| present(field(args, n)).is_some());
    let interval = if explicit {
        Some(PosteriorIntervalOptions {
            kind,
            level: number(args, "interval.level", 0.95)?,
            grid_intervals: count(args, "ngrid", 50)?,
            grid_standard_deviations: number(args, "nsd", 5.0)?,
            nuisance_points: count(args, "ngrid.nuis", 5)?,
            nuisance_standard_deviations: number(args, "nsd.nuis", 2.0)?,
        })
    } else {
        None
    };
    let contrasts = present(field(args, "contrasts"))
        .map(|v| {
            let [rows, columns] = v.dimensions.ok_or("contrasts requires a matrix")?;
            RowMajorMatrix::from_row_major(rows, columns, numbers(v)?.to_vec())
                .map_err(|e| e.to_string())
        })
        .transpose()?;
    let parameter_standard_errors = present(field(args, "param.sd"))
        .map(|v| numbers(v).map(<[f64]>::to_vec))
        .transpose()?;
    let grid_likelihood = if kind != PosteriorIntervalKind::Laplace {
        present(field(args, "log.lik"))
            .map(likelihood)
            .transpose()?
    } else {
        None
    };
    let mut work = rsdeseq2::shrinkage::PosteriorWorkControl::default()
        .with_interrupt_check(super::poll_interrupt);
    if let Some(value) = field(args, "max.grid.evaluations") {
        work.max_grid_evaluations = if matches!(value.data, Data::Null) {
            None
        } else {
            Some(count(args, "max.grid.evaluations", 10_000_000)?)
        };
    }
    Ok(PosteriorAnalysisOptions {
        work,
        interval,
        contrasts,
        parameter_standard_errors,
        grid_likelihood,
    })
}
fn matrix(v: Option<RowMajorMatrix<f64>>) -> Value {
    v.map(|m| {
        let dimensions = Some([m.n_rows(), m.n_cols()]);
        let mut v = Value::numbers(m.into_values());
        v.dimensions = dimensions;
        v
    })
    .unwrap_or_else(|| Value::new(Data::Null))
}
pub(super) fn encode(analysis: PosteriorAnalysis) -> Value {
    named(
        &["intervals", "contrast_mode", "contrast_sd", "grid"],
        vec![
            matrix(analysis.intervals),
            matrix(analysis.contrast_mode),
            matrix(analysis.contrast_standard_deviation),
            Value::new(Data::List(
                analysis
                    .grid_diagnostics
                    .into_iter()
                    .map(|d| {
                        d.map(|d| {
                            named(
                                &["monotone", "left", "right", "evaluations"],
                                vec![
                                    Value::new(Data::Logical(vec![Some(d.monotone)])),
                                    Value::numbers(vec![
                                        d.outside_left
                                            .unwrap_or(f64::from_bits(0x7ff00000000007a2)),
                                    ]),
                                    Value::numbers(vec![
                                        d.outside_right
                                            .unwrap_or(f64::from_bits(0x7ff00000000007a2)),
                                    ]),
                                    Value::numbers(vec![d.evaluations as f64]),
                                ],
                            )
                        })
                        .unwrap_or_else(|| Value::new(Data::Null))
                    })
                    .collect(),
            )),
        ],
    )
}
