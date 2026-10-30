//! Prepared standalone batches for direct interface conformance checks.
use rsdeseq2::prelude::*;
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;

#[derive(Deserialize)]
pub(super) struct BatchInput {
    genes: usize,
    coefficients: usize,
    beta: Vec<String>,
    se: Option<Vec<String>>,
    covariance: Option<Vec<String>>,
    samples: usize,
    counts: Vec<u32>,
    design: Option<Vec<String>>,
    levels: Vec<String>,
    requests: Vec<Request>,
}

#[derive(Deserialize)]
struct Request {
    contrast: Vec<String>,
    threshold: Option<String>,
    alternative: String,
    df: Option<Vec<String>>,
    zero: String,
    numerator: String,
    denominator: String,
}

fn parse(values: Vec<String>) -> Result<Vec<f64>, std::num::ParseFloatError> {
    values.into_iter().map(|v| v.parse()).collect()
}

pub(super) fn evaluate(input: BatchInput) -> Result<Value, Box<dyn Error>> {
    let beta = parse(input.beta)?;
    let se = input.se.map(parse).transpose()?;
    let covariance = input.covariance.map(parse).transpose()?;
    let matrix = |values| {
        MatrixView::new(
            input.genes,
            input.coefficients,
            values,
            MatrixOrder::RowMajor,
        )
    };
    let covariance = covariance
        .as_deref()
        .map(|values| {
            CovarianceView::new(
                input.coefficients,
                MatrixView::new(
                    input.genes,
                    input
                        .coefficients
                        .checked_mul(input.coefficients)
                        .ok_or("covariance shape overflow")?,
                    values,
                    MatrixOrder::RowMajor,
                )?,
            )
            .map_err(|e| -> Box<dyn Error> { Box::new(e) })
        })
        .transpose()?;
    let estimates = InferenceEstimates::new(
        matrix(&beta)?,
        se.as_deref().map(matrix).transpose()?,
        covariance,
    )?;
    let counts = CountMatrix::from_row_major_u32(input.genes, input.samples, input.counts)?;
    let design = input.design.map(parse).transpose()?;
    let design = design
        .as_deref()
        .map(|values| {
            MatrixView::new(
                input.samples,
                input.coefficients,
                values,
                MatrixOrder::RowMajor,
            )
        })
        .transpose()?;
    let view = PreparedContrastView::new(estimates, Some(ContrastCounts::from(&counts)), design)?;
    let levels: Vec<_> = input.levels.iter().map(String::as_str).collect();
    let mut results = Vec::with_capacity(input.requests.len());
    for request in input.requests {
        let alternative = match request.alternative.as_str() {
            "greaterAbs" => WaldAlternative::GreaterAbs,
            "greaterAbsUPSHOT" => WaldAlternative::GreaterAbsUpshot,
            "greaterAbs2014" => WaldAlternative::GreaterAbs2014,
            "lessAbs" => WaldAlternative::LessAbs,
            "greater" => WaldAlternative::Greater,
            "less" => WaldAlternative::Less,
            _ => return Err("unknown Wald alternative".into()),
        };
        let options = if let Some(threshold) = request.threshold {
            Some(
                match request.df.map(parse).transpose()? {
                    None => WaldTestOptions::normal(),
                    Some(df) if df.len() == 1 => WaldTestOptions::t_degrees_of_freedom(df[0]),
                    Some(df) => WaldTestOptions::t_per_gene_degrees_of_freedom(df),
                }
                .with_lfc_threshold(threshold.parse()?, alternative),
            )
        } else {
            None
        };
        let zero = match request.zero.as_str() {
            "none" => ContrastZeroRule::None,
            "numeric" => ContrastZeroRule::Numeric,
            "character" => ContrastZeroRule::Factor {
                levels: &levels,
                numerator: &request.numerator,
                denominator: &request.denominator,
            },
            _ => return Err("unknown contrast zero rule".into()),
        };
        let result = view.evaluate(&parse(request.contrast)?, options.as_ref(), zero)?;
        let hex = |values: Vec<Option<f64>>| {
            values
                .into_iter()
                .map(|v| v.map(super::bits))
                .collect::<Vec<_>>()
        };
        let mut result_json = json!({
            "log2FoldChange": hex(result.projection.log2_fold_change),
            "lfcSE": hex(result.projection.lfc_se),
        });
        if let Some(wald) = result.wald {
            result_json["stat"] = json!(hex(wald.stat));
            result_json["pvalue"] = json!(hex(wald.pvalue));
        }
        results.push(result_json);
    }
    Ok(json!({"results": results}))
}
