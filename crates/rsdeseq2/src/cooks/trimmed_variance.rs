use super::*;

#[derive(Clone, Copy)]
struct TrimmedVariance {
    trim: f64,
    scale: f64,
    mean: &'static str,
    residual_square: &'static str,
    variance_mean: &'static str,
    variance: &'static str,
}

impl TrimmedVariance {
    const ALL_SAMPLES: Self = Self {
        trim: 1.0 / 8.0,
        scale: 1.51,
        mean: "Cook's trimmed mean",
        residual_square: "Cook's trimmed residual square",
        variance_mean: "Cook's trimmed variance mean",
        variance: "Cook's trimmed variance",
    };

    fn for_cell(samples: usize) -> Self {
        let (trim, scale) = if samples <= 3 {
            (1.0 / 3.0, 2.04)
        } else if samples <= 23 {
            (1.0 / 4.0, 1.86)
        } else {
            (1.0 / 8.0, 1.51)
        };
        Self {
            trim,
            scale,
            mean: "Cook's trimmed cell mean",
            residual_square: "Cook's trimmed cell residual square",
            variance_mean: "Cook's trimmed cell variance mean",
            variance: "Cook's trimmed cell variance",
        }
    }

    fn evaluate(
        self,
        values: impl Iterator<Item = f64> + Clone,
        sorted: &mut Vec<f64>,
        gene: usize,
    ) -> Result<f64, DeseqError> {
        sorted.clear();
        sorted.extend(values.clone());
        let mean = trimmed_scaled_mean(sorted, self.trim, self.mean)?;
        sorted.clear();
        // Keep the original sample traversal for residual errors before sorting.
        for (sample, value) in values.enumerate() {
            let residual = value - mean;
            sorted.push(
                checked::try_mul(residual, residual)
                    .nonfinite_at(self.residual_square, Some(sample))?,
            );
        }
        checked::try_mul(
            self.scale,
            trimmed_scaled_mean(sorted, self.trim, self.variance_mean)?,
        )
        .nonfinite_at(self.variance, Some(gene))
    }
}

pub(super) fn trimmed_variances(
    normalized_counts: &RowMajorMatrix<f64>,
    groups: &[Vec<usize>],
) -> Result<Vec<f64>, DeseqError> {
    let mut replicated = try_alloc_matrix(groups.len(), 1, "Cook's replicated cells")?;
    replicated.extend(
        groups
            .iter()
            .filter(|group| group.len() >= 3)
            .map(|group| (group.as_slice(), TrimmedVariance::for_cell(group.len()))),
    );
    let samples = normalized_counts.n_cols();
    let mut variances = try_filled_vec(normalized_counts.n_rows(), f64::NAN, "Cook's variances")?;
    let first_error = variances
        .par_iter_mut()
        .zip(normalized_counts.as_slice().par_chunks(samples))
        .enumerate()
        .map_init(
            || None,
            |scratch, (gene, (variance, row))| -> Result<(), DeseqError> {
                let sorted = cooks_row_buffer(scratch, samples)?;
                *variance = if replicated.is_empty() {
                    TrimmedVariance::ALL_SAMPLES.evaluate(row.iter().copied(), sorted, gene)?
                } else {
                    let mut maximum = f64::NEG_INFINITY;
                    for &(group, rule) in &replicated {
                        let values = group.iter().map(|&sample| row[sample]);
                        maximum = maximum.max(rule.evaluate(values, sorted, gene)?);
                    }
                    maximum
                };
                Ok(())
            },
        )
        .find_map_first(Result::err);
    if let Some(error) = first_error {
        return Err(error);
    }
    Ok(variances)
}

/// Trim the sorted central values without applying R's full-trim median rule.
fn trimmed_scaled_mean(values: &mut [f64], trim: f64, context: &str) -> Result<f64, DeseqError> {
    if values.is_empty() || values.iter().any(|value| value.is_nan()) {
        return Ok(f64::NAN);
    }
    values.sort_unstable_by(f64::total_cmp);
    let trim_count = (values.len() as f64 * trim).floor() as usize;
    let end = values.len().saturating_sub(trim_count);
    if trim_count >= end {
        return Ok(f64::NAN);
    }
    checked::try_mean(&values[trim_count..end]).nonfinite_at(context, None)
}
