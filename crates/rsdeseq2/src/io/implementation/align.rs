/// Align sample-level factor values to count-matrix sample order.
pub fn align_sample_levels_to_samples(
    levels: &[SampleLevel],
    sample_names: &[String],
) -> Result<Vec<String>, DeseqError> {
    if levels.len() != sample_names.len() {
        return Err(DeseqError::InvalidDimensions {
            context: "sample-level rows".to_string(),
            expected: sample_names.len(),
            actual: levels.len(),
        });
    }
    let mut by_sample = HashMap::with_capacity(levels.len());
    for level in levels {
        if by_sample
            .insert(level.sample.as_str(), level.level.as_str())
            .is_some()
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!("duplicate sample level for {}", level.sample),
            });
        }
    }
    let mut aligned = Vec::with_capacity(sample_names.len());
    for sample in sample_names {
        let level = by_sample
            .get(sample.as_str())
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!("missing sample level for {sample}"),
            })?;
        aligned.push((*level).to_string());
    }
    Ok(aligned)
}

/// Align sample-level numeric values to count-matrix sample order.
pub fn align_sample_numeric_values_to_samples(
    values: &[SampleNumericValue],
    sample_names: &[String],
    context: &str,
) -> Result<Vec<f64>, DeseqError> {
    if values.len() != sample_names.len() {
        return Err(DeseqError::InvalidDimensions {
            context: format!("{context} rows"),
            expected: sample_names.len(),
            actual: values.len(),
        });
    }
    let mut by_sample = HashMap::with_capacity(values.len());
    for value in values {
        if by_sample
            .insert(value.sample.as_str(), value.value)
            .is_some()
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!("duplicate {context} value for sample {}", value.sample),
            });
        }
    }
    let mut aligned = Vec::with_capacity(sample_names.len());
    for sample in sample_names {
        let value = by_sample
            .get(sample.as_str())
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!("missing {context} value for sample {sample}"),
            })?;
        aligned.push(*value);
    }
    Ok(aligned)
}

/// Align gene-level numeric values to count-matrix gene order.
pub fn align_gene_numeric_values_to_genes(
    values: &[GeneNumericValue],
    gene_names: &[String],
    context: &str,
) -> Result<Vec<f64>, DeseqError> {
    if values.len() != gene_names.len() {
        return Err(DeseqError::InvalidDimensions {
            context: format!("{context} rows"),
            expected: gene_names.len(),
            actual: values.len(),
        });
    }
    let mut by_gene = HashMap::with_capacity(values.len());
    for value in values {
        if by_gene.insert(value.gene.as_str(), value.value).is_some() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("duplicate {context} value for gene {}", value.gene),
            });
        }
    }
    let mut aligned = Vec::with_capacity(gene_names.len());
    for gene in gene_names {
        let value = by_gene
            .get(gene.as_str())
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!("missing {context} value for gene {gene}"),
            })?;
        aligned.push(*value);
    }
    Ok(aligned)
}

/// Align a labeled gene x sample matrix to count-matrix gene and sample order.
pub fn align_labeled_assay_matrix_to_counts(
    labeled: LabeledAssayMatrix,
    counts: &CountMatrix,
    context: &str,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let gene_names = counts
        .gene_names()
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!("count gene names are required to align {context}"),
        })?;
    let sample_names = counts
        .sample_names()
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!("count sample names are required to align {context}"),
        })?;
    if labeled.gene_names.len() != gene_names.len() {
        return Err(DeseqError::InvalidDimensions {
            context: format!("{context} gene rows"),
            expected: gene_names.len(),
            actual: labeled.gene_names.len(),
        });
    }
    if labeled.sample_names.len() != sample_names.len() {
        return Err(DeseqError::InvalidDimensions {
            context: format!("{context} sample columns"),
            expected: sample_names.len(),
            actual: labeled.sample_names.len(),
        });
    }
    let mut gene_rows = HashMap::with_capacity(labeled.gene_names.len());
    for (idx, gene) in labeled.gene_names.iter().enumerate() {
        if gene_rows.insert(gene.as_str(), idx).is_some() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("duplicate {context} row for gene {gene}"),
            });
        }
    }
    let mut sample_cols = HashMap::with_capacity(labeled.sample_names.len());
    for (idx, sample) in labeled.sample_names.iter().enumerate() {
        if sample_cols.insert(sample.as_str(), idx).is_some() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("duplicate {context} column for sample {sample}"),
            });
        }
    }
    let mut values = Vec::with_capacity(gene_names.len() * sample_names.len());
    for gene in gene_names {
        let row = gene_rows
            .get(gene.as_str())
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!("missing {context} row for gene {gene}"),
            })?;
        for sample in sample_names {
            let col =
                sample_cols
                    .get(sample.as_str())
                    .ok_or_else(|| DeseqError::InvalidOptions {
                        reason: format!("missing {context} column for sample {sample}"),
                    })?;
            values.push(*labeled.matrix.get(*row, *col).ok_or_else(|| {
                DeseqError::InvalidDimensions {
                    context: format!("{context} aligned value"),
                    expected: labeled.matrix.len(),
                    actual: *row * labeled.matrix.n_cols() + *col,
                }
            })?);
        }
    }
    RowMajorMatrix::from_row_major(gene_names.len(), sample_names.len(), values)
}

/// Align a labeled design matrix to count-matrix sample order.
pub fn align_design_matrix_to_samples(
    labeled: LabeledDesignMatrix,
    sample_names: &[String],
) -> Result<DesignMatrix, DeseqError> {
    if labeled.sample_names.len() != sample_names.len() {
        return Err(DeseqError::InvalidDimensions {
            context: "design sample rows".to_string(),
            expected: sample_names.len(),
            actual: labeled.sample_names.len(),
        });
    }
    let mut by_sample = HashMap::with_capacity(labeled.sample_names.len());
    for (row, sample) in labeled.sample_names.iter().enumerate() {
        if by_sample.insert(sample.as_str(), row).is_some() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("duplicate design row for sample {sample}"),
            });
        }
    }
    let n_coefficients = labeled.design.n_coefficients();
    let mut values = Vec::with_capacity(sample_names.len() * n_coefficients);
    for sample in sample_names {
        let row = by_sample
            .get(sample.as_str())
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!("missing design row for sample {sample}"),
            })?;
        values.extend_from_slice(labeled.design.matrix().row(*row)?);
    }
    DesignMatrix::from_row_major(
        sample_names.len(),
        n_coefficients,
        values,
        labeled.design.coefficient_names().map(<[String]>::to_vec),
    )
}
