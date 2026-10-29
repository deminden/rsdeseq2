//! Standalone numerical conformance driver; it neither links nor launches R.
use rsdeseq2::independent_filtering::{
    IndependentFilteringOptions, evaluate_independent_filtering,
};
use serde::Deserialize;
use serde_json::json;
use std::error::Error;

#[derive(Deserialize)]
struct FilteringInput {
    filter: Vec<String>,
    pvalues: Vec<Option<String>>,
    theta: Option<Vec<String>>,
    alpha: String,
    enabled: bool,
}

fn bits(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

fn filtering(input: FilteringInput) -> Result<serde_json::Value, Box<dyn Error>> {
    let parse = |values: Vec<String>| -> Result<Vec<f64>, std::num::ParseFloatError> {
        values.into_iter().map(|v| v.parse()).collect()
    };
    let filter = parse(input.filter)?;
    let pvalues = input
        .pvalues
        .into_iter()
        .map(|v| v.map(|v| v.parse()).transpose())
        .collect::<Result<Vec<_>, _>>()?;
    let options = IndependentFilteringOptions {
        enabled: input.enabled,
        alpha: input.alpha.parse()?,
        theta: input.theta.map(parse).transpose()?,
    };
    let evaluated = evaluate_independent_filtering(&filter, &pvalues, &options)?;
    let metadata = evaluated.metadata;
    Ok(json!({
            "padj": evaluated.padj.into_iter().map(|v| v.map(bits)).collect::<Vec<_>>(),
            "filtered": evaluated.filtered,
            "enabled": metadata.enabled,
            "theta": metadata.theta.iter().copied().map(bits).collect::<Vec<_>>(),
            "numRej": metadata.num_rejections,
            "selected": metadata.selected_index.map(|v| v + 1),
            "filterTheta": metadata.filter_theta.map(bits),
            "filterThreshold": metadata.filter_threshold.map(bits),
            "lowessX": metadata.lowess_fit.as_ref().map(|_| metadata.theta.iter().copied().map(bits).collect::<Vec<_>>()),
            "lowessY": metadata.lowess_fit.map(|v| v.into_iter().map(bits).collect::<Vec<_>>()),
            "alpha": bits(metadata.alpha),
    }))
}

#[derive(Deserialize)]
struct NormalizationInput {
    genes: usize,
    samples: usize,
    counts: Vec<u32>,
    method: String,
    geo_means: Option<Vec<String>>,
    control_genes: Option<Vec<usize>>,
    size_factors: Option<Vec<String>>,
    normalization_factors: Option<Vec<String>>,
    weights: Option<Vec<String>>,
}

fn normalization(input: NormalizationInput) -> Result<serde_json::Value, Box<dyn Error>> {
    use rsdeseq2::normalization::{
        base_mean, base_mean_with_weights, base_variance, base_variance_with_weights,
        estimate_size_factors_with_options, normalized_counts, normalized_counts_with_factors,
    };
    use rsdeseq2::prelude::{CountMatrix, RowMajorMatrix, SizeFactorMethod};
    let parse = |values: Vec<String>| {
        values
            .into_iter()
            .map(|v| v.parse::<f64>())
            .collect::<Result<Vec<_>, _>>()
    };
    let counts = CountMatrix::from_row_major_u32(input.genes, input.samples, input.counts)?;
    let geo = input.geo_means.map(parse).transpose()?;
    let method = match input.method.as_str() {
        "ratio" => SizeFactorMethod::Ratio,
        "poscounts" => SizeFactorMethod::PosCounts,
        _ => return Err("unknown size-factor method".into()),
    };
    let factors = estimate_size_factors_with_options(
        &counts,
        method,
        geo.as_deref(),
        input.control_genes.as_deref(),
    )?;
    let matrix = |values| -> Result<RowMajorMatrix<f64>, Box<dyn Error>> {
        Ok(RowMajorMatrix::from_row_major(
            input.genes,
            input.samples,
            parse(values)?,
        )?)
    };
    let normalization_factors = input.normalization_factors.map(matrix).transpose()?;
    let weights = input.weights.map(matrix).transpose()?;
    let supplied = input.size_factors.map(parse).transpose()?;
    let normalized = match normalization_factors.as_ref() {
        Some(matrix) => normalized_counts_with_factors(&counts, matrix)?,
        None => normalized_counts(&counts, supplied.as_deref().unwrap_or(&factors))?,
    };
    let (means, variances) = match weights.as_ref() {
        Some(weights) => (
            base_mean_with_weights(&normalized, weights)?,
            base_variance_with_weights(&normalized, weights)?,
        ),
        None => (base_mean(&normalized)?, base_variance(&normalized)?),
    };
    Ok(json!({
        "sizeFactors": factors.into_iter().map(bits).collect::<Vec<_>>(),
        "normalized": normalized.as_slice().iter().copied().map(bits).collect::<Vec<_>>(),
        "baseMean": means.into_iter().map(bits).collect::<Vec<_>>(),
        "baseVar": variances.into_iter().map(bits).collect::<Vec<_>>(),
        "allZero": counts.all_zero_flags(),
    }))
}

#[derive(Deserialize)]
#[serde(untagged)]
enum Input {
    Filtering(FilteringInput),
    Normalization(NormalizationInput),
}

fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("expected an input JSON path")?;
    let input: Input = serde_json::from_slice(&std::fs::read(path)?)?;
    let result = match input {
        Input::Filtering(input) => filtering(input)?,
        Input::Normalization(input) => normalization(input)?,
    };
    println!("{result}");
    Ok(())
}
