use super::*;
// Fit and normalization-factor validation.

pub(super) fn validate_intercept_fit_inputs(input: GlmInput<'_>) -> Result<(), DeseqError> {
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
            if factors.n_rows() != counts.n_genes() || factors.n_cols() != counts.n_samples() {
                return Err(DeseqError::InvalidDimensions {
                    context: "normalization factors".to_string(),
                    expected: counts.n_genes() * counts.n_samples(),
                    actual: factors.len(),
                });
            }
        }
    }
    if dispersions.len() != counts.n_genes() {
        return Err(invalid_dimensions(
            "dispersions",
            counts.n_genes(),
            dispersions.len(),
        ));
    }
    for (idx, dispersion) in dispersions.iter().copied().enumerate() {
        validate_positive_finite(dispersion, "dispersion", idx)?;
    }
    if let Some(weights) = weights {
        if weights.n_rows() != counts.n_genes() || weights.n_cols() != counts.n_samples() {
            return Err(DeseqError::InvalidDimensions {
                context: "weights".to_string(),
                expected: counts.n_genes() * counts.n_samples(),
                actual: weights.len(),
            });
        }
        weights.validate_finite("weights")?;
    }
    if let Offsets::PerCell(factors) = offsets {
        factors.validate_finite("normalization factors")?;
        for (index, &factor) in factors.as_slice().iter().enumerate() {
            validate_positive_finite(factor, "normalization factor", index % counts.n_samples())?;
        }
    }
    Ok(())
}
