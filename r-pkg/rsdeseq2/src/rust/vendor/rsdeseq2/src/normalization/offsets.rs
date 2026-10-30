//! One base-factor/offset composition for native formula and foreign interfaces.
use super::validate_size_factors;
use crate::errors::{DeseqError, invalid_dimensions};
use crate::matrix::{MatrixView, RowMajorMatrix};

/// Compose effective normalization factors as `base * exp(offset)`.
///
/// Gene/sample factors take precedence over sample size factors. Every supplied
/// input is validated; arithmetic uses that one multiplication in both layouts.
/// Offsets must already be aligned to samples. This function does not estimate
/// size factors, apply offsets twice, or normalize by the base factors first.
pub fn normalization_factors_with_offsets(
    genes: usize,
    samples: usize,
    size_factors: Option<&[f64]>,
    factors: Option<MatrixView<'_, f64>>,
    offsets: &[f64],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    if offsets.len() != samples {
        return Err(invalid_dimensions(
            "formula offsets",
            samples,
            offsets.len(),
        ));
    }
    let length = genes
        .checked_mul(samples)
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: "normalization offset shape overflow".into(),
        })?;
    if let Some(sf) = size_factors {
        validate_size_factors(sf, samples)?;
    }
    if let Some(nf) = factors {
        if nf.n_rows() != genes || nf.n_cols() != samples {
            return Err(invalid_dimensions(
                "normalization factors",
                length,
                nf.as_slice().len(),
            ));
        }
        for (index, &value) in nf.as_slice().iter().enumerate() {
            if !value.is_finite() || value <= 0.0 {
                return Err(DeseqError::InvalidSizeFactors {
                    reason: format!(
                        "normalization factor at index {index} must be finite and positive"
                    ),
                });
            }
        }
    } else if size_factors.is_none() {
        return Err(DeseqError::InvalidSizeFactors {
            reason: "offset composition requires base size or normalization factors".into(),
        });
    }
    let multipliers=offsets.iter().enumerate().map(|(sample,&offset)| {
        let multiplier=offset.exp();
        if offset.is_finite() && multiplier.is_finite() && multiplier>0.0 { Ok(multiplier) }
        else { Err(DeseqError::InvalidOptions {reason:format!("formula offset exponentiated to non-finite or non-positive factor at sample {sample}")}) }
    }).collect::<Result<Vec<_>,_>>()?;
    let mut values = Vec::with_capacity(length);
    for gene in 0..genes {
        for (sample, &multiplier) in multipliers.iter().enumerate() {
            let base = if let Some(nf) = factors {
                *nf.at(gene, sample)
            } else {
                size_factors.expect("validated base factors")[sample]
            };
            let value = base * multiplier;
            if !value.is_finite() || value <= 0.0 {
                return Err(DeseqError::InvalidSizeFactors {
                    reason: format!(
                        "effective normalization factor at gene {gene}, sample {sample} must be finite and positive"
                    ),
                });
            }
            values.push(value);
        }
    }
    RowMajorMatrix::from_row_major(genes, samples, values)
}
