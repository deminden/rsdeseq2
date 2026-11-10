//! Descriptor input only; never a foreign function callback.
use rsdeseq2::shrinkage::{RestrictedFunction, RestrictedLikelihood};
use serde_json::Value;
use std::error::Error;
fn function(value: &Value) -> Result<RestrictedFunction, Box<dyn Error>> {
    let arguments = value["arguments"]
        .as_array()
        .ok_or("likelihood argument names must be an array")?
        .iter()
        .map(|s| {
            s.as_str()
                .map(str::to_owned)
                .ok_or("invalid likelihood argument name")
        })
        .collect::<Result<Vec<_>, _>>()?;
    let constants = value["constants"]
        .as_object()
        .map(|v| {
            v.iter()
                .map(|(name, value)| (name.clone(), super::formula::value(value)))
                .collect()
        })
        .unwrap_or_default();
    Ok(RestrictedFunction::new(
        arguments,
        super::formula::expression(&value["expression"]),
        constants,
    )?)
}
pub fn restricted(value: &Value) -> Result<RestrictedLikelihood, Box<dyn Error>> {
    let optional = |name| {
        if value[name].is_null() {
            Ok(None)
        } else {
            function(&value[name]).map(Some)
        }
    };
    Ok(RestrictedLikelihood::new(
        value["name"]
            .as_str()
            .ok_or("missing likelihood name")?
            .into(),
        function(&value["log_cdf"])?,
        function(&value["log_density"])?,
        optional("first")?,
        optional("second")?,
    )?)
}
