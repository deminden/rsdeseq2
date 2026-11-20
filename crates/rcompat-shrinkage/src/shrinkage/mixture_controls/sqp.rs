//! SQP control names translate to the shared typed fitting options.
use super::{Data, MixtureFitOptions, Value, fields, flag, scalar};
use crate::{Error, errors::invalid_options};
pub(super) fn controls(
    value: Option<&Value>,
    options: &mut MixtureFitOptions,
) -> Result<(), Error> {
    let Some(value) = value else {
        return Ok(());
    };
    if matches!(&value.data,Data::List(v) if v.is_empty()) {
        return Ok(());
    }
    let (names, values) = fields(value)?;
    for (name, value) in names.iter().zip(values) {
        let count = || super::count(scalar(value)?, name);
        match name.as_str() {
            // The core API takes log likelihoods and normalizes them before
            // optimization, as original ashr does before calling mixsqp.
            "verbose" | "normalize.rows" => {
                flag(value)?;
            }
            "eps" => options.likelihood_offset = scalar(value)?,
            "numiter.em" => options.em_iterations = count()?,
            "maxiter.sqp" => options.max_iterations = count()?,
            "maxiter.activeset" => {
                if matches!(value.data, Data::Null) {
                    options.sqp.automatic_active_iterations = true;
                } else {
                    options.active_iterations = count()?;
                    options.sqp.automatic_active_iterations = false;
                }
            }
            "convtol.sqp" => options.gradient_tolerance = scalar(value)?,
            "zero.threshold.solution" => options.zero_threshold = scalar(value)?,
            "convtol.activeset" => options.sqp.active_tolerance = scalar(value)?,
            "zero.threshold.searchdir" => {
                options.sqp.search_direction_threshold = Some(scalar(value)?)
            }
            "suffdecr.linesearch" => options.sqp.sufficient_decrease = scalar(value)?,
            "stepsizereduce" => options.sqp.step_reduction = scalar(value)?,
            "minstepsize" => options.sqp.minimum_step = scalar(value)?,
            "identity.contrib.increase" => options.sqp.ridge_multiplier = scalar(value)?,
            "tol.svd" => options.sqp.svd_tolerance = Some(scalar(value)?),
            _ => return Err(invalid_options(format!("unknown SQP control '{name}'"))),
        }
    }
    Ok(())
}
