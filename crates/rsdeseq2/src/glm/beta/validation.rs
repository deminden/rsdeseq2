use super::*;
use crate::errors::validate_matrix_shape;
// Fit and normalization-factor validation.

pub(super) fn validate_intercept_fit_inputs(
    input: GlmInput<'_>,
    genes: Option<&[usize]>,
) -> Result<(), DeseqError> {
    let GlmInput {
        counts,
        offsets,
        dispersions,
        weights,
    } = input;
    match offsets {
        Offsets::PerSample(factors) => {
            validate_size_factor_offsets(counts, factors)?;
        }
        Offsets::PerCell(factors) => {
            validate_matrix_shape(
                "normalization factors",
                (counts.n_genes(), counts.n_samples()),
                (factors.n_rows(), factors.n_cols()),
            )?;
        }
    }
    if dispersions.len() != counts.n_genes() {
        return Err(invalid_dimensions(
            "dispersions",
            counts.n_genes(),
            dispersions.len(),
        ));
    }
    let n_genes = genes.map_or(counts.n_genes(), <[usize]>::len);
    let source_genes = || (0..n_genes).map(|gene| genes.map_or(gene, |rows| rows[gene]));
    for (gene, source_gene) in source_genes().enumerate() {
        let dispersion = *dispersions.get(source_gene).ok_or_else(|| {
            invalid_dimensions(
                "intercept gene index",
                counts.n_genes(),
                source_gene.saturating_add(1),
            )
        })?;
        validate_positive_finite(dispersion, "dispersion", gene)?;
    }
    if let Some(weights) = weights {
        validate_matrix_shape(
            "weights",
            (counts.n_genes(), counts.n_samples()),
            (weights.n_rows(), weights.n_cols()),
        )?;
        validate_intercept_finite_rows(weights, "weights", genes)?;
    }
    if let Offsets::PerCell(factors) = offsets {
        validate_intercept_finite_rows(factors, "normalization factors", genes)?;
        for source_gene in source_genes() {
            for (sample, &factor) in factors.row(source_gene)?.iter().enumerate() {
                validate_positive_finite(factor, "normalization factor", sample)?;
            }
        }
    }
    Ok(())
}

fn validate_intercept_finite_rows(
    matrix: &RowMajorMatrix<f64>,
    context: &str,
    genes: Option<&[usize]>,
) -> Result<(), DeseqError> {
    let Some(genes) = genes else {
        return matrix.validate_finite(context);
    };
    for (gene, &source_gene) in genes.iter().enumerate() {
        for (sample, &value) in matrix.row(source_gene)?.iter().enumerate() {
            if !value.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: context.into(),
                    index: Some(gene * matrix.n_cols() + sample),
                    value,
                });
            }
        }
    }
    Ok(())
}
