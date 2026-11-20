//! One base-factor/offset composition for native formula and foreign interfaces.
use super::{validate_normalization_factor_value, validate_size_factors};
use crate::errors::missing_prerequisite;
use crate::errors::{DeseqError, invalid_dimensions, validate_matrix_shape};
use crate::errors::{ErrorKind, Requirement, invalid_formula, invalid_parameter};
use crate::matrix::{MatrixView, RowMajorMatrix, try_alloc_matrix};
use crate::numeric::math::ProfileMath;

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
    let (base, multipliers) = prepare_offsets(genes, samples, size_factors, factors, offsets)?;
    let mut values = try_alloc_matrix(genes, samples, "effective normalization factors")?;
    match base {
        BaseFactors::Sample(size_factors) => {
            let row = sample_products(size_factors, multipliers)?;
            for _ in 0..genes {
                values.extend_from_slice(&row);
            }
        }
        BaseFactors::Cell(factors) => {
            for gene in 0..genes {
                for (sample, &multiplier) in multipliers.iter().enumerate() {
                    let value = *factors.at(gene, sample) * multiplier;
                    if !valid_factor(value) {
                        return Err(invalid_product(gene, sample));
                    }
                    values.push(value);
                }
            }
        }
    }
    RowMajorMatrix::from_row_major(genes, samples, values)
}

/// Stream effective factors in column-major order. The caller discards output
/// on error. Validation and arithmetic match the owned matrix path.
#[doc(hidden)]
pub fn visit_offset_factors_column_major(
    genes: usize,
    samples: usize,
    size_factors: Option<&[f64]>,
    factors: Option<MatrixView<'_, f64>>,
    offsets: &[f64],
    mut write: impl FnMut(usize, f64),
) -> Result<(), DeseqError> {
    let (base, multipliers) = prepare_offsets(genes, samples, size_factors, factors, offsets)?;
    match base {
        BaseFactors::Sample(size_factors) => {
            let products = sample_products(size_factors, multipliers)?;
            for (sample, value) in products.into_iter().enumerate() {
                for gene in 0..genes {
                    write(sample * genes + gene, value);
                }
            }
        }
        BaseFactors::Cell(factors) => {
            let mut first_error = None;
            let mut genes_to_check = genes;
            for (sample, multiplier) in multipliers.into_iter().enumerate() {
                for gene in 0..genes_to_check {
                    let value = *factors.at(gene, sample) * multiplier;
                    if !valid_factor(value) {
                        first_error = Some((gene, sample));
                        // In later columns, only an earlier row can take precedence.
                        genes_to_check = gene;
                        break;
                    }
                    write(sample * genes + gene, value);
                }
            }
            if let Some((gene, sample)) = first_error {
                return Err(invalid_product(gene, sample));
            }
        }
    }
    Ok(())
}

enum BaseFactors<'a> {
    Sample(&'a [f64]),
    Cell(MatrixView<'a, f64>),
}

fn prepare_offsets<'a>(
    genes: usize,
    samples: usize,
    size_factors: Option<&'a [f64]>,
    factors: Option<MatrixView<'a, f64>>,
    offsets: &[f64],
) -> Result<(BaseFactors<'a>, Vec<f64>), DeseqError> {
    validate_offset_count(offsets, samples)?;
    genes
        .checked_mul(samples)
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: "normalization offset shape overflow".into(),
        })?;
    if let Some(sf) = size_factors {
        validate_size_factors(sf, samples)?;
    }
    let base = if let Some(nf) = factors {
        validate_matrix_shape(
            "normalization factors",
            (genes, samples),
            (nf.n_rows(), nf.n_cols()),
        )?;
        for (index, &value) in nf.as_slice().iter().enumerate() {
            validate_normalization_factor_value(value, index)?;
        }
        BaseFactors::Cell(nf)
    } else if let Some(size_factors) = size_factors {
        BaseFactors::Sample(size_factors)
    } else {
        return Err(missing_prerequisite(
            ErrorKind::Normalization,
            "offset composition",
            "base size or normalization factors",
        ));
    };
    let multipliers = offset_multipliers(offsets)?;
    if genes == 0 {
        return Err(invalid_dimensions("row-major matrix rows", 1, 0));
    }
    if samples == 0 {
        return Err(invalid_dimensions("row-major matrix columns", 1, 0));
    }
    Ok((base, multipliers))
}

pub(crate) fn sample_factors_with_offsets(
    size_factors: &[f64],
    offsets: &[f64],
    samples: usize,
) -> Result<Vec<f64>, DeseqError> {
    validate_offset_count(offsets, samples)?;
    validate_size_factors(size_factors, samples)?;
    sample_products(size_factors, offset_multipliers(offsets)?)
}

fn validate_offset_count(offsets: &[f64], samples: usize) -> Result<(), DeseqError> {
    if offsets.len() != samples {
        return Err(invalid_dimensions(
            "formula offsets",
            samples,
            offsets.len(),
        ));
    }
    Ok(())
}

fn offset_multipliers(offsets: &[f64]) -> Result<Vec<f64>, DeseqError> {
    let mut multipliers = try_alloc_matrix(offsets.len(), 1, "formula offset multipliers")?;
    for (sample, &offset) in offsets.iter().enumerate() {
        let multiplier = offset.exp_p();
        if !offset.is_finite() || !valid_factor(multiplier) {
            return Err(invalid_formula(format!(
                "formula offset exponentiated to non-finite or non-positive factor at sample {sample}"
            )));
        }
        multipliers.push(multiplier);
    }
    Ok(multipliers)
}

fn sample_products(
    size_factors: &[f64],
    mut multipliers: Vec<f64>,
) -> Result<Vec<f64>, DeseqError> {
    for (sample, (multiplier, &size_factor)) in multipliers.iter_mut().zip(size_factors).enumerate()
    {
        let value = size_factor * *multiplier;
        if !valid_factor(value) {
            return Err(invalid_product(0, sample));
        }
        *multiplier = value;
    }
    Ok(multipliers)
}

fn valid_factor(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn invalid_product(gene: usize, sample: usize) -> DeseqError {
    invalid_parameter(
        ErrorKind::Normalization,
        format!("effective normalization factor at gene {gene}, sample {sample}"),
        Requirement::FiniteAndPositive,
    )
}
