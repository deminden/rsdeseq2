//! Exact test-only expanded-factor primitive panel.
use rsdeseq2::design::expression::{EvaluationOptions, frame};
use serde::Deserialize;
use serde_json::{Value, json};

#[derive(Deserialize)]
pub struct Case {
    name: String,
    input: Value,
    null_level: String,
}
pub fn evaluate(cases: &[Case]) -> Value {
    json!(cases.iter().map(|case| {
        let value = super::formula::value(&case.input);
        match frame::expanded_factors(&[value], &case.null_level, &EvaluationOptions { collation: "C".into() }) {
            Ok(values) => json!({"name":case.name,"ok":true,"value":super::formula::descriptor(&values[0])}),
            Err(_) => json!({"name":case.name,"ok":false}),
        }
    }).collect::<Vec<_>>())
}
