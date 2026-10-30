//! Standalone offset products and normal-shrink row policy; no R runtime.
use rsdeseq2::matrix::{CountRows, MatrixOrder, MatrixView};
use rsdeseq2::prelude::{NormalShrinkageRows, RowFlags, normalization_factors_with_offsets};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;
#[derive(Deserialize)]
pub(super) struct Input {
    preparation: Preparation,
}
#[derive(Deserialize)]
#[serde(tag = "kind")]
enum Preparation {
    Offset {
        genes: usize,
        samples: usize,
        size_factors: Option<Vec<String>>,
        factors: Option<Vec<String>>,
        offsets: Vec<String>,
    },
    Rows {
        genes: usize,
        samples: usize,
        counts: Vec<u32>,
        dispersions: Vec<String>,
        fitted: Vec<String>,
        excluded: Option<Vec<bool>>,
    },
}
fn numbers(values: Vec<String>) -> Result<Vec<f64>, Box<dyn Error>> {
    values
        .into_iter()
        .map(|v| Ok(f64::from_bits(u64::from_str_radix(&v, 16)?)))
        .collect()
}
pub(super) fn evaluate(input: Input) -> Result<Value, Box<dyn Error>> {
    let result = (|| -> Result<Value, Box<dyn Error>> {
        match input.preparation {
            Preparation::Offset {
                genes,
                samples,
                size_factors,
                factors,
                offsets,
            } => {
                let sf = size_factors.map(numbers).transpose()?;
                let nf = factors.map(numbers).transpose()?;
                let nf = nf
                    .as_ref()
                    .map(|v| MatrixView::new(genes, samples, v, MatrixOrder::RowMajor))
                    .transpose()?;
                let output = normalization_factors_with_offsets(
                    genes,
                    samples,
                    sf.as_deref(),
                    nf,
                    &numbers(offsets)?,
                )?;
                Ok(
                    json!({"factors":output.as_slice().iter().map(|&v|super::bits(v)).collect::<Vec<_>>()}),
                )
            }
            Preparation::Rows {
                genes,
                samples,
                counts,
                dispersions,
                fitted,
                excluded,
            } => {
                let rows = CountRows::from_unsigned(MatrixView::new(
                    genes,
                    samples,
                    &counts,
                    MatrixOrder::RowMajor,
                )?);
                let output = NormalShrinkageRows::new(
                    &rows,
                    &numbers(dispersions)?,
                    &numbers(fitted)?,
                    excluded.as_deref().map(RowFlags::from_bools),
                )?;
                Ok(json!({"indices":output.indices().iter().map(|i|i+1).collect::<Vec<_>>()}))
            }
        }
    })();
    Ok(result.unwrap_or_else(|_| json!({"error":true})))
}
