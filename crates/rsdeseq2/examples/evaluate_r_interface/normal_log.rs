use rsdeseq2::math::{normal_log_cdf, normal_log_upper_tail};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;
#[derive(Deserialize)]
pub struct Input {
    normal_log_arguments: Vec<String>,
}
pub fn evaluate(input: Input) -> Result<Value, Box<dyn Error>> {
    let values = super::parse_numeric_inputs(&input.normal_log_arguments)?;
    Ok(
        json!({"log_cdf":values.iter().map(|v|super::bits(normal_log_cdf(*v))).collect::<Vec<_>>(),
        "log_upper":values.iter().map(|v|super::bits(normal_log_upper_tail(*v))).collect::<Vec<_>>()}),
    )
}
