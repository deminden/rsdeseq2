use std::collections::{HashMap, HashSet};
use std::path::Path;

use csv::{ReaderBuilder, WriterBuilder};

use crate::cooks::{CooksOutput, CooksRefitPlan};
use crate::core::CountMatrix;
use crate::design::DesignMatrix;
use crate::diagnostics::{
    Deseq2McolsDiagnosticColumn, Deseq2McolsDiagnosticValues, Deseq2McolsDiagnostics,
};
use crate::errors::{invalid_dimensions, DeseqError};
use crate::independent_filtering::IndependentFilteringOutput;
use crate::matrix::RowMajorMatrix;
use crate::results::{DeseqResultColumnValues, DeseqResults};

/// One sample-level factor value loaded from a two-column TSV.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SampleLevel {
    /// Sample label.
    pub sample: String,
    /// Factor level for the sample.
    pub level: String,
}

/// One finite numeric value associated with a sample label.
#[derive(Clone, Debug, PartialEq)]
pub struct SampleNumericValue {
    /// Sample label.
    pub sample: String,
    /// Numeric value for the sample.
    pub value: f64,
}

/// One finite numeric value associated with a gene label.
#[derive(Clone, Debug, PartialEq)]
pub struct GeneNumericValue {
    /// Gene label.
    pub gene: String,
    /// Numeric value for the gene.
    pub value: f64,
}

/// Numeric design matrix plus the row sample labels read from TSV.
#[derive(Clone, Debug, PartialEq)]
pub struct LabeledDesignMatrix {
    /// Design matrix in file row order.
    pub design: DesignMatrix,
    /// One sample label per design row.
    pub sample_names: Vec<String>,
}

/// Gene x sample numeric matrix plus row and column labels read from TSV.
#[derive(Clone, Debug, PartialEq)]
pub struct LabeledAssayMatrix {
    /// Matrix in file row and column order.
    pub matrix: RowMajorMatrix<f64>,
    /// One gene label per matrix row.
    pub gene_names: Vec<String>,
    /// One sample label per matrix column.
    pub sample_names: Vec<String>,
}

/// Read a tab-delimited count matrix with gene IDs in the first column.
pub fn read_count_matrix_tsv(path: impl AsRef<Path>) -> Result<CountMatrix, DeseqError> {
    let path_ref = path.as_ref();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path_ref)?;
    let headers = reader.headers()?.clone();
    if headers.len() < 2 {
        return Err(DeseqError::InvalidCounts {
            reason: "count table must have a gene column and at least one sample column"
                .to_string(),
        });
    }
    let sample_names = headers
        .iter()
        .skip(1)
        .map(ToOwned::to_owned)
        .collect::<Vec<_>>();
    let mut gene_names = Vec::new();
    let mut values = Vec::new();
    for record in reader.records() {
        let record = record?;
        if record.len() != headers.len() {
            return Err(DeseqError::InvalidCounts {
                reason: format!(
                    "record has {} columns but header has {} columns",
                    record.len(),
                    headers.len()
                ),
            });
        }
        gene_names.push(record.get(0).unwrap_or_default().to_string());
        for field in record.iter().skip(1) {
            values.push(parse_count_field(field, path_ref)?);
        }
    }
    CountMatrix::from_row_major_u32_with_names(
        gene_names.len(),
        sample_names.len(),
        values,
        Some(gene_names),
        Some(sample_names),
    )
}

/// Read a tab-delimited numeric design matrix with sample IDs in the first column.
///
/// The first column is treated as row labels and ignored by the primitive Rust
/// core; remaining header fields become coefficient names. Rows must already be
/// in the same sample order as the count matrix used by the caller.
pub fn read_design_matrix_tsv(path: impl AsRef<Path>) -> Result<DesignMatrix, DeseqError> {
    Ok(read_labeled_design_matrix_tsv(path)?.design)
}

/// Read a tab-delimited numeric design matrix and preserve sample row labels.
pub fn read_labeled_design_matrix_tsv(
    path: impl AsRef<Path>,
) -> Result<LabeledDesignMatrix, DeseqError> {
    let path_ref = path.as_ref();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path_ref)?;
    let headers = reader.headers()?.clone();
    if headers.len() < 2 {
        return Err(DeseqError::InvalidDimensions {
            context: "design matrix columns".to_string(),
            expected: 2,
            actual: headers.len(),
        });
    }
    let coefficient_names = headers
        .iter()
        .skip(1)
        .map(ToOwned::to_owned)
        .collect::<Vec<_>>();
    let n_coefficients = coefficient_names.len();
    let mut sample_names = Vec::new();
    let mut values = Vec::new();
    let mut n_samples = 0_usize;
    for record in reader.records() {
        let record = record?;
        if record.len() != headers.len() {
            return Err(DeseqError::InvalidDimensions {
                context: "design matrix record columns".to_string(),
                expected: headers.len(),
                actual: record.len(),
            });
        }
        let sample = record.get(0).unwrap_or_default();
        if sample.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: "design sample names must not be empty".to_string(),
            });
        }
        sample_names.push(sample.to_string());
        for (column, field) in record.iter().skip(1).enumerate() {
            values.push(parse_finite_float_field(
                field,
                path_ref,
                n_samples * n_coefficients + column,
                "design matrix",
            )?);
        }
        n_samples += 1;
    }
    if n_samples == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "design matrix rows".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    let design =
        DesignMatrix::from_row_major(n_samples, n_coefficients, values, Some(coefficient_names))?;
    Ok(LabeledDesignMatrix {
        design,
        sample_names,
    })
}

/// Read a tab-delimited gene/sample normalization-factor matrix.
///
/// The file uses the same shape as `normalizationFactors(dds)`: a leading
/// `gene` column followed by sample columns. Row and column labels are accepted
/// for alignment by the caller; this primitive reader returns only the numeric
/// matrix.
pub fn read_normalization_factors_tsv(
    path: impl AsRef<Path>,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    Ok(read_labeled_normalization_factors_tsv(path)?.matrix)
}

/// Read a tab-delimited gene/sample normalization-factor matrix with labels preserved.
pub fn read_labeled_normalization_factors_tsv(
    path: impl AsRef<Path>,
) -> Result<LabeledAssayMatrix, DeseqError> {
    let path_ref = path.as_ref();
    let factors = read_labeled_assay_matrix_tsv(
        path_ref,
        "normalization factor columns",
        "normalization factor record columns",
        "normalization factor rows",
        "normalization factor matrix",
        "normalization factor gene names",
    )?;
    validate_positive_finite_matrix("normalization factor matrix", &factors.matrix)?;
    Ok(factors)
}

/// Read a tab-delimited gene/sample observation-weight matrix.
///
/// The file uses a leading `gene` column followed by sample columns, matching
/// DESeq2's assay-shaped weight matrix. Row and column labels are accepted for
/// alignment by the caller; this primitive reader returns only numeric weights
/// in row-major order.
pub fn read_observation_weights_tsv(
    path: impl AsRef<Path>,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    Ok(read_labeled_observation_weights_tsv(path)?.matrix)
}

/// Read a tab-delimited gene/sample observation-weight matrix with labels preserved.
pub fn read_labeled_observation_weights_tsv(
    path: impl AsRef<Path>,
) -> Result<LabeledAssayMatrix, DeseqError> {
    let path_ref = path.as_ref();
    let weights = read_labeled_assay_matrix_tsv(
        path_ref,
        "observation weight columns",
        "observation weight record columns",
        "observation weight rows",
        "observation weights",
        "observation weight gene names",
    )?;
    validate_nonnegative_finite_values("observation weight input", weights.matrix.as_slice())?;
    Ok(weights)
}

fn read_labeled_assay_matrix_tsv(
    path_ref: &Path,
    column_context: &str,
    record_context: &str,
    row_context: &str,
    value_context: &str,
    gene_context: &str,
) -> Result<LabeledAssayMatrix, DeseqError> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path_ref)?;
    let headers = reader.headers()?.clone();
    if headers.len() < 2 {
        return Err(DeseqError::InvalidDimensions {
            context: column_context.to_string(),
            expected: 2,
            actual: headers.len(),
        });
    }
    let sample_names = headers
        .iter()
        .skip(1)
        .map(ToOwned::to_owned)
        .collect::<Vec<_>>();
    let n_samples = headers.len() - 1;
    let mut gene_names = Vec::new();
    let mut values = Vec::new();
    let mut n_genes = 0_usize;
    for record in reader.records() {
        let record = record?;
        if record.len() != headers.len() {
            return Err(DeseqError::InvalidDimensions {
                context: record_context.to_string(),
                expected: headers.len(),
                actual: record.len(),
            });
        }
        let gene = record.get(0).unwrap_or_default();
        if gene.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("{gene_context} must not be empty"),
            });
        }
        gene_names.push(gene.to_string());
        for (sample, field) in record.iter().skip(1).enumerate() {
            values.push(parse_finite_float_field(
                field,
                path_ref,
                n_genes * n_samples + sample,
                value_context,
            )?);
        }
        n_genes += 1;
    }
    if n_genes == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: row_context.to_string(),
            expected: 1,
            actual: 0,
        });
    }
    let matrix = RowMajorMatrix::from_row_major(n_genes, n_samples, values)?;
    Ok(LabeledAssayMatrix {
        matrix,
        gene_names,
        sample_names,
    })
}

/// Read a tab-delimited sample-level size-factor table.
///
/// The file shape matches `write_size_factors_tsv`: a leading sample label
/// column and one numeric `size_factor` column. Sample labels are accepted for
/// alignment by the caller; this primitive reader returns only the numeric
/// vector in file order.
pub fn read_size_factors_tsv(path: impl AsRef<Path>) -> Result<Vec<f64>, DeseqError> {
    Ok(read_labeled_size_factors_tsv(path)?
        .into_iter()
        .map(|entry| entry.value)
        .collect())
}

/// Read a tab-delimited sample-level size-factor table with labels preserved.
pub fn read_labeled_size_factors_tsv(
    path: impl AsRef<Path>,
) -> Result<Vec<SampleNumericValue>, DeseqError> {
    let path_ref = path.as_ref();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path_ref)?;
    let headers = reader.headers()?.clone();
    if headers.len() != 2 {
        return Err(DeseqError::InvalidDimensions {
            context: "size-factor columns".to_string(),
            expected: 2,
            actual: headers.len(),
        });
    }
    let mut values = Vec::new();
    for (idx, record) in reader.records().enumerate() {
        let record = record?;
        if record.len() != 2 {
            return Err(DeseqError::InvalidDimensions {
                context: "size-factor record columns".to_string(),
                expected: 2,
                actual: record.len(),
            });
        }
        let sample = record.get(0).unwrap_or_default();
        if sample.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: "size-factor sample names must not be empty".to_string(),
            });
        }
        let value = parse_finite_float_field(
            record.get(1).unwrap_or_default(),
            path_ref,
            idx,
            "size factors",
        )?;
        values.push(SampleNumericValue {
            sample: sample.to_string(),
            value,
        });
    }
    if values.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "size-factor rows".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    validate_positive_finite_values(
        "size-factor input",
        &values.iter().map(|entry| entry.value).collect::<Vec<_>>(),
    )?;
    Ok(values)
}

/// Read one supplied geometric mean per gene.
///
/// The file has a leading gene label column and one numeric value column.
/// Labels are accepted for caller-side alignment; this primitive reader returns
/// the values in file order. Values must be non-negative and finite, matching
/// the size-factor estimator's supported geometric-mean domain.
pub fn read_geometric_means_tsv(path: impl AsRef<Path>) -> Result<Vec<f64>, DeseqError> {
    Ok(read_labeled_geometric_means_tsv(path)?
        .into_iter()
        .map(|entry| entry.value)
        .collect())
}

/// Read one supplied geometric mean per gene with labels preserved.
pub fn read_labeled_geometric_means_tsv(
    path: impl AsRef<Path>,
) -> Result<Vec<GeneNumericValue>, DeseqError> {
    let path_ref = path.as_ref();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path_ref)?;
    let headers = reader.headers()?.clone();
    if headers.len() != 2 {
        return Err(DeseqError::InvalidDimensions {
            context: "geometric-mean columns".to_string(),
            expected: 2,
            actual: headers.len(),
        });
    }
    let mut values = Vec::new();
    for (idx, record) in reader.records().enumerate() {
        let record = record?;
        if record.len() != 2 {
            return Err(DeseqError::InvalidDimensions {
                context: "geometric-mean record columns".to_string(),
                expected: 2,
                actual: record.len(),
            });
        }
        let gene = record.get(0).unwrap_or_default();
        if gene.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: "geometric-mean gene names must not be empty".to_string(),
            });
        }
        let value = parse_finite_float_field(
            record.get(1).unwrap_or_default(),
            path_ref,
            idx,
            "geometric means",
        )?;
        values.push(GeneNumericValue {
            gene: gene.to_string(),
            value,
        });
    }
    if values.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "geometric-mean rows".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    for (idx, value) in values.iter().map(|entry| entry.value).enumerate() {
        if value < 0.0 {
            return Err(DeseqError::InvalidSizeFactors {
                reason: format!("geometric mean at index {idx} must be non-negative"),
            });
        }
    }
    Ok(values)
}

/// Read one finite numeric value per gene with labels preserved.
pub fn read_labeled_gene_numeric_tsv(
    path: impl AsRef<Path>,
    context: &str,
) -> Result<Vec<GeneNumericValue>, DeseqError> {
    let path_ref = path.as_ref();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path_ref)?;
    let headers = reader.headers()?.clone();
    if headers.len() != 2 {
        return Err(DeseqError::InvalidDimensions {
            context: format!("{context} columns"),
            expected: 2,
            actual: headers.len(),
        });
    }
    let mut values = Vec::new();
    for (idx, record) in reader.records().enumerate() {
        let record = record?;
        if record.len() != 2 {
            return Err(DeseqError::InvalidDimensions {
                context: format!("{context} record columns"),
                expected: 2,
                actual: record.len(),
            });
        }
        let gene = record.get(0).unwrap_or_default();
        if gene.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("{context} gene names must not be empty"),
            });
        }
        let value =
            parse_finite_float_field(record.get(1).unwrap_or_default(), path_ref, idx, context)?;
        values.push(GeneNumericValue {
            gene: gene.to_string(),
            value,
        });
    }
    if values.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: format!("{context} rows"),
            expected: 1,
            actual: 0,
        });
    }
    Ok(values)
}

/// Read one finite numeric Wald t degrees-of-freedom value per gene.
///
/// The file has a leading gene label column and one numeric value column.
/// Labels are accepted for caller-side alignment; this primitive reader returns
/// the values in file order.
pub fn read_wald_t_degrees_of_freedom_tsv(path: impl AsRef<Path>) -> Result<Vec<f64>, DeseqError> {
    Ok(read_labeled_wald_t_degrees_of_freedom_tsv(path)?
        .into_iter()
        .map(|entry| entry.value)
        .collect())
}

/// Read one finite numeric Wald t degrees-of-freedom value per gene with labels preserved.
pub fn read_labeled_wald_t_degrees_of_freedom_tsv(
    path: impl AsRef<Path>,
) -> Result<Vec<GeneNumericValue>, DeseqError> {
    let path_ref = path.as_ref();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path_ref)?;
    let headers = reader.headers()?.clone();
    if headers.len() != 2 {
        return Err(DeseqError::InvalidDimensions {
            context: "Wald t degrees-of-freedom columns".to_string(),
            expected: 2,
            actual: headers.len(),
        });
    }
    let mut values = Vec::new();
    for (idx, record) in reader.records().enumerate() {
        let record = record?;
        if record.len() != 2 {
            return Err(DeseqError::InvalidDimensions {
                context: "Wald t degrees-of-freedom record columns".to_string(),
                expected: 2,
                actual: record.len(),
            });
        }
        let gene = record.get(0).unwrap_or_default();
        if gene.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: "Wald t degrees-of-freedom gene names must not be empty".to_string(),
            });
        }
        let value = parse_finite_float_field(
            record.get(1).unwrap_or_default(),
            path_ref,
            idx,
            "Wald t degrees of freedom",
        )?;
        values.push(GeneNumericValue {
            gene: gene.to_string(),
            value,
        });
    }
    if values.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "Wald t degrees-of-freedom rows".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    Ok(values)
}

/// Read one factor level per sample for CLI factor-level contrasts.
///
/// The file has a leading sample label column and one string level column.
/// Labels are preserved so callers can align levels against count-matrix sample
/// columns before applying DESeq2-style factor-level contrast handling.
pub fn read_sample_levels_tsv(path: impl AsRef<Path>) -> Result<Vec<SampleLevel>, DeseqError> {
    let path_ref = path.as_ref();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path_ref)?;
    let headers = reader.headers()?.clone();
    if headers.len() != 2 {
        return Err(DeseqError::InvalidDimensions {
            context: "sample-level columns".to_string(),
            expected: 2,
            actual: headers.len(),
        });
    }
    let mut values = Vec::new();
    for record in reader.records() {
        let record = record?;
        if record.len() != 2 {
            return Err(DeseqError::InvalidDimensions {
                context: "sample-level record columns".to_string(),
                expected: 2,
                actual: record.len(),
            });
        }
        let sample = record.get(0).unwrap_or_default();
        if sample.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: "sample names must not be empty".to_string(),
            });
        }
        let level = record.get(1).unwrap_or_default();
        if level.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: "sample levels must not be empty".to_string(),
            });
        }
        values.push(SampleLevel {
            sample: sample.to_string(),
            level: level.to_string(),
        });
    }
    if values.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "sample-level rows".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    Ok(values)
}

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

fn parse_count_field(field: &str, path: &Path) -> Result<u32, DeseqError> {
    if let Ok(value) = field.parse::<u32>() {
        return Ok(value);
    }

    let value = field.parse::<f64>().map_err(|_| DeseqError::ParseInt {
        context: path.display().to_string(),
        value: field.to_string(),
    })?;
    if !value.is_finite() {
        return Err(DeseqError::NonFiniteValue {
            context: path.display().to_string(),
            index: None,
            value,
        });
    }
    if value < 0.0 || value.fract() != 0.0 || value > f64::from(u32::MAX) {
        return Err(DeseqError::InvalidCounts {
            reason: format!("count value '{field}' must be a non-negative integer"),
        });
    }
    Ok(value as u32)
}

fn parse_finite_float_field(
    field: &str,
    path: &Path,
    index: usize,
    context: &str,
) -> Result<f64, DeseqError> {
    let value = field.parse::<f64>().map_err(|_| DeseqError::ParseFloat {
        context: path.display().to_string(),
        value: field.to_string(),
    })?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value,
        })
    }
}

/// Write a raw count matrix to a tab-delimited file.
///
/// The first column is `gene`, followed by sample columns. Stored row and
/// sample names are used when present; otherwise `gene1`, `gene2`, ... and
/// `sample1`, `sample2`, ... fallback labels are emitted.
pub fn write_count_matrix_tsv(
    path: impl AsRef<Path>,
    counts: &CountMatrix,
) -> Result<(), DeseqError> {
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(matrix_header(
        counts.sample_names(),
        counts.n_samples(),
        "sample",
    ))?;
    for gene in 0..counts.n_genes() {
        let mut record = Vec::with_capacity(counts.n_samples() + 1);
        record.push(matrix_row_name(counts.gene_names(), gene, "gene"));
        record.extend(counts.row_values(gene).iter().map(ToString::to_string));
        writer.write_record(&record)?;
    }
    writer.flush()?;
    Ok(())
}

/// Write DESeq2-style normalized counts to a tab-delimited file.
///
/// This matches the matrix shape returned by DESeq2
/// `counts(dds, normalized = TRUE)`: genes are rows, samples are columns, and
/// gene identifiers are written as the first `gene` column.
pub fn write_normalized_counts_tsv(
    path: impl AsRef<Path>,
    gene_names: Option<&[String]>,
    sample_names: Option<&[String]>,
    normalized_counts: &RowMajorMatrix<f64>,
) -> Result<(), DeseqError> {
    validate_optional_names(
        "normalized count gene names",
        gene_names,
        normalized_counts.n_rows(),
    )?;
    validate_optional_names(
        "normalized count sample names",
        sample_names,
        normalized_counts.n_cols(),
    )?;
    normalized_counts.validate_finite("normalized count export")?;

    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(matrix_header(
        sample_names,
        normalized_counts.n_cols(),
        "sample",
    ))?;
    for gene in 0..normalized_counts.n_rows() {
        let mut record = Vec::with_capacity(normalized_counts.n_cols() + 1);
        record.push(matrix_row_name(gene_names, gene, "gene"));
        record.extend(
            normalized_counts
                .row(gene)?
                .iter()
                .copied()
                .map(|value| value.to_string()),
        );
        writer.write_record(&record)?;
    }
    writer.flush()?;
    Ok(())
}

/// Write gene/sample normalization factors to a tab-delimited file.
///
/// DESeq2 uses `normalizationFactors(dds)` as a genes x samples count-scale
/// factor matrix that preempts size factors. This writer preserves that shape
/// with a leading `gene` column and sample-name columns.
pub fn write_normalization_factors_tsv(
    path: impl AsRef<Path>,
    gene_names: Option<&[String]>,
    sample_names: Option<&[String]>,
    normalization_factors: &RowMajorMatrix<f64>,
) -> Result<(), DeseqError> {
    validate_optional_names(
        "normalization factor gene names",
        gene_names,
        normalization_factors.n_rows(),
    )?;
    validate_optional_names(
        "normalization factor sample names",
        sample_names,
        normalization_factors.n_cols(),
    )?;
    validate_positive_finite_matrix("normalization factor export", normalization_factors)?;

    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(matrix_header(
        sample_names,
        normalization_factors.n_cols(),
        "sample",
    ))?;
    for gene in 0..normalization_factors.n_rows() {
        let mut record = Vec::with_capacity(normalization_factors.n_cols() + 1);
        record.push(matrix_row_name(gene_names, gene, "gene"));
        record.extend(
            normalization_factors
                .row(gene)?
                .iter()
                .copied()
                .map(|value| value.to_string()),
        );
        writer.write_record(&record)?;
    }
    writer.flush()?;
    Ok(())
}

/// Write sample-level size factors to a tab-delimited file.
pub fn write_size_factors_tsv(
    path: impl AsRef<Path>,
    sample_names: Option<&[String]>,
    size_factors: &[f64],
) -> Result<(), DeseqError> {
    validate_optional_names("size-factor sample names", sample_names, size_factors.len())?;
    validate_positive_finite_values("size-factor export", size_factors)?;

    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["sample", "size_factor"])?;
    for (idx, size_factor) in size_factors.iter().copied().enumerate() {
        let sample = sample_names
            .and_then(|names| names.get(idx))
            .cloned()
            .unwrap_or_else(|| format!("sample{}", idx + 1));
        writer.write_record([sample, size_factor.to_string()])?;
    }
    writer.flush()?;
    Ok(())
}

/// Write gene-level base means to a tab-delimited file.
pub fn write_base_mean_tsv(
    path: impl AsRef<Path>,
    gene_names: Option<&[String]>,
    base_mean: &[f64],
) -> Result<(), DeseqError> {
    validate_optional_names("base-mean gene names", gene_names, base_mean.len())?;
    validate_nonnegative_finite_values("base-mean export", base_mean)?;

    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["gene", "base_mean"])?;
    for (idx, value) in base_mean.iter().copied().enumerate() {
        let gene = gene_names
            .and_then(|names| names.get(idx))
            .cloned()
            .unwrap_or_else(|| format!("gene{}", idx + 1));
        writer.write_record([gene, value.to_string()])?;
    }
    writer.flush()?;
    Ok(())
}

/// Write DESeq2-style early row metadata to a tab-delimited file.
///
/// The output mirrors the primitive `mcols(dds)` row metadata shape for
/// `baseMean`, `baseVar`, and `allZero`. Non-finite numeric entries, such as
/// one-sample `baseVar`, are written as `NA`.
pub fn write_base_metadata_tsv(
    path: impl AsRef<Path>,
    gene_names: Option<&[String]>,
    base_mean: &[f64],
    base_variance: &[f64],
    all_zero: &[bool],
) -> Result<(), DeseqError> {
    validate_optional_names("base metadata gene names", gene_names, base_mean.len())?;
    if base_variance.len() != base_mean.len() {
        return Err(invalid_dimensions(
            "base metadata baseVar rows",
            base_mean.len(),
            base_variance.len(),
        ));
    }
    if all_zero.len() != base_mean.len() {
        return Err(invalid_dimensions(
            "base metadata allZero rows",
            base_mean.len(),
            all_zero.len(),
        ));
    }

    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["gene", "baseMean", "baseVar", "allZero"])?;
    for gene in 0..base_mean.len() {
        writer.write_record([
            matrix_row_name(gene_names, gene, "gene"),
            format_finite_or_na(base_mean[gene]),
            format_finite_or_na(base_variance[gene]),
            format_r_logical(all_zero[gene]),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

/// Write the Cook's distance assay to a tab-delimited numeric matrix.
///
/// Missing/non-finite Cook's distances are written as `NA`, matching R-style
/// diagnostic table exports.
pub fn write_cooks_distance_matrix_tsv(
    path: impl AsRef<Path>,
    gene_names: Option<&[String]>,
    sample_names: Option<&[String]>,
    cooks: &RowMajorMatrix<f64>,
) -> Result<(), DeseqError> {
    write_optional_numeric_matrix_tsv(path, gene_names, sample_names, cooks)
}

/// Write the Cook's distance assay from a full Cook's diagnostic output.
///
/// Missing/non-finite Cook's distances are written as `NA`, matching R-style
/// diagnostic table exports.
pub fn write_cooks_distance_tsv(
    path: impl AsRef<Path>,
    gene_names: Option<&[String]>,
    sample_names: Option<&[String]>,
    cooks: &CooksOutput,
) -> Result<(), DeseqError> {
    write_optional_numeric_matrix_tsv(path, gene_names, sample_names, &cooks.cooks)
}

/// Write row-level Cook's diagnostic metadata to a tab-delimited table.
pub fn write_cooks_row_metadata_tsv(
    path: impl AsRef<Path>,
    gene_names: Option<&[String]>,
    cooks: &CooksOutput,
) -> Result<(), DeseqError> {
    validate_optional_names(
        "Cook's row metadata gene names",
        gene_names,
        cooks.cooks.n_rows(),
    )?;
    if cooks.max_cooks.len() != cooks.cooks.n_rows() {
        return Err(invalid_dimensions(
            "Cook's row metadata maxCooks rows",
            cooks.cooks.n_rows(),
            cooks.max_cooks.len(),
        ));
    }
    if cooks.robust_dispersion.len() != cooks.cooks.n_rows() {
        return Err(invalid_dimensions(
            "Cook's row metadata robust dispersion rows",
            cooks.cooks.n_rows(),
            cooks.robust_dispersion.len(),
        ));
    }

    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["gene", "maxCooks", "cooksRobustDispersion"])?;
    for gene in 0..cooks.cooks.n_rows() {
        writer.write_record([
            matrix_row_name(gene_names, gene, "gene"),
            format_optional_finite_or_na(cooks.max_cooks[gene]),
            format_finite_or_na(cooks.robust_dispersion[gene]),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

/// Write sample-level Cook's diagnostic metadata to a tab-delimited table.
pub fn write_cooks_sample_metadata_tsv(
    path: impl AsRef<Path>,
    sample_names: Option<&[String]>,
    cooks: &CooksOutput,
) -> Result<(), DeseqError> {
    validate_optional_names(
        "Cook's sample metadata sample names",
        sample_names,
        cooks.cooks.n_cols(),
    )?;
    if cooks.samples_for_cooks.len() != cooks.cooks.n_cols() {
        return Err(invalid_dimensions(
            "Cook's sample metadata samplesForCooks",
            cooks.cooks.n_cols(),
            cooks.samples_for_cooks.len(),
        ));
    }

    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["sample", "samplesForCooks"])?;
    for sample in 0..cooks.cooks.n_cols() {
        writer.write_record([
            matrix_sample_name(sample_names, sample, "sample"),
            format_r_logical(cooks.samples_for_cooks[sample]),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

/// Write Cook's replacement-count assay to a tab-delimited count matrix.
pub fn write_cooks_replaced_counts_tsv(
    path: impl AsRef<Path>,
    refit_plan: &CooksRefitPlan,
) -> Result<(), DeseqError> {
    write_count_matrix_tsv(path, &refit_plan.replacement.replaced_counts)
}

/// Write Cook's candidate replacement-count assay to a tab-delimited count matrix.
pub fn write_cooks_candidate_replacement_counts_tsv(
    path: impl AsRef<Path>,
    refit_plan: &CooksRefitPlan,
) -> Result<(), DeseqError> {
    write_count_matrix_tsv(path, &refit_plan.replacement.candidate_replacement_counts)
}

/// Write per-cell Cook's outlier flags to a tab-delimited logical matrix.
pub fn write_cooks_outlier_cells_tsv(
    path: impl AsRef<Path>,
    refit_plan: &CooksRefitPlan,
) -> Result<(), DeseqError> {
    write_logical_matrix_tsv(
        path,
        refit_plan.replacement.replaced_counts.gene_names(),
        refit_plan.replacement.replaced_counts.sample_names(),
        &refit_plan.replacement.outlier_cells,
    )
}

/// Write row-level Cook's replacement/refit metadata to a tab-delimited table.
pub fn write_cooks_replacement_row_metadata_tsv(
    path: impl AsRef<Path>,
    refit_plan: &CooksRefitPlan,
) -> Result<(), DeseqError> {
    let n_rows = refit_plan.replacement.replaced_counts.n_genes();
    validate_cooks_refit_plan_rows(refit_plan, n_rows)?;
    let refit_rows = refit_plan
        .refit_rows
        .iter()
        .copied()
        .collect::<HashSet<_>>();
    let new_all_zero_rows = refit_plan
        .new_all_zero_rows
        .iter()
        .copied()
        .collect::<HashSet<_>>();

    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record([
        "gene",
        "replace",
        "refitReplace",
        "newAllZero",
        "replacedAllZero",
        "replacedBaseMean",
        "replacedBaseVar",
        "postRefitMaxCooks",
    ])?;
    for gene in 0..n_rows {
        writer.write_record([
            matrix_row_name(
                refit_plan.replacement.replaced_counts.gene_names(),
                gene,
                "gene",
            ),
            format_optional_r_logical(refit_plan.replacement.replace[gene]),
            format_r_logical(refit_rows.contains(&gene)),
            format_r_logical(new_all_zero_rows.contains(&gene)),
            format_r_logical(refit_plan.replaced_all_zero[gene]),
            format_finite_or_na(refit_plan.replaced_base_mean[gene]),
            format_finite_or_na(refit_plan.replaced_base_var[gene]),
            format_optional_finite_or_na(refit_plan.post_refit_max_cooks[gene]),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

fn write_logical_matrix_tsv(
    path: impl AsRef<Path>,
    gene_names: Option<&[String]>,
    sample_names: Option<&[String]>,
    matrix: &RowMajorMatrix<bool>,
) -> Result<(), DeseqError> {
    validate_optional_names("logical matrix gene names", gene_names, matrix.n_rows())?;
    validate_optional_names("logical matrix sample names", sample_names, matrix.n_cols())?;

    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(matrix_header(sample_names, matrix.n_cols(), "sample"))?;
    for gene in 0..matrix.n_rows() {
        let mut record = Vec::with_capacity(matrix.n_cols() + 1);
        record.push(matrix_row_name(gene_names, gene, "gene"));
        record.extend(matrix.row(gene)?.iter().copied().map(format_r_logical));
        writer.write_record(&record)?;
    }
    writer.flush()?;
    Ok(())
}

fn write_optional_numeric_matrix_tsv(
    path: impl AsRef<Path>,
    gene_names: Option<&[String]>,
    sample_names: Option<&[String]>,
    matrix: &RowMajorMatrix<f64>,
) -> Result<(), DeseqError> {
    validate_optional_names("numeric matrix gene names", gene_names, matrix.n_rows())?;
    validate_optional_names("numeric matrix sample names", sample_names, matrix.n_cols())?;

    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(matrix_header(sample_names, matrix.n_cols(), "sample"))?;
    for gene in 0..matrix.n_rows() {
        let mut record = Vec::with_capacity(matrix.n_cols() + 1);
        record.push(matrix_row_name(gene_names, gene, "gene"));
        record.extend(matrix.row(gene)?.iter().copied().map(format_finite_or_na));
        writer.write_record(&record)?;
    }
    writer.flush()?;
    Ok(())
}

fn validate_cooks_refit_plan_rows(
    refit_plan: &CooksRefitPlan,
    expected_rows: usize,
) -> Result<(), DeseqError> {
    if refit_plan.replacement.replace.len() != expected_rows {
        return Err(invalid_dimensions(
            "Cook's replacement replace rows",
            expected_rows,
            refit_plan.replacement.replace.len(),
        ));
    }
    if refit_plan.replaced_all_zero.len() != expected_rows {
        return Err(invalid_dimensions(
            "Cook's replacement all-zero rows",
            expected_rows,
            refit_plan.replaced_all_zero.len(),
        ));
    }
    if refit_plan.replaced_base_mean.len() != expected_rows {
        return Err(invalid_dimensions(
            "Cook's replacement baseMean rows",
            expected_rows,
            refit_plan.replaced_base_mean.len(),
        ));
    }
    if refit_plan.replaced_base_var.len() != expected_rows {
        return Err(invalid_dimensions(
            "Cook's replacement baseVar rows",
            expected_rows,
            refit_plan.replaced_base_var.len(),
        ));
    }
    if refit_plan.post_refit_max_cooks.len() != expected_rows {
        return Err(invalid_dimensions(
            "Cook's replacement post-refit maxCooks rows",
            expected_rows,
            refit_plan.post_refit_max_cooks.len(),
        ));
    }
    Ok(())
}

fn matrix_header(
    sample_names: Option<&[String]>,
    n_samples: usize,
    fallback_prefix: &str,
) -> Vec<String> {
    let mut header = Vec::with_capacity(n_samples + 1);
    header.push("gene".to_string());
    header.extend(
        (0..n_samples).map(|sample| matrix_sample_name(sample_names, sample, fallback_prefix)),
    );
    header
}

fn matrix_row_name(gene_names: Option<&[String]>, gene: usize, fallback_prefix: &str) -> String {
    gene_names
        .and_then(|names| names.get(gene))
        .cloned()
        .unwrap_or_else(|| format!("{}{}", fallback_prefix, gene + 1))
}

fn matrix_sample_name(
    sample_names: Option<&[String]>,
    sample: usize,
    fallback_prefix: &str,
) -> String {
    sample_names
        .and_then(|names| names.get(sample))
        .cloned()
        .unwrap_or_else(|| format!("{}{}", fallback_prefix, sample + 1))
}

fn validate_optional_names(
    context: &str,
    names: Option<&[String]>,
    expected: usize,
) -> Result<(), DeseqError> {
    if let Some(names) = names {
        if names.len() != expected {
            return Err(invalid_dimensions(context, expected, names.len()));
        }
    }
    Ok(())
}

fn validate_positive_finite_matrix(
    context: &str,
    matrix: &RowMajorMatrix<f64>,
) -> Result<(), DeseqError> {
    validate_positive_finite_values(context, matrix.as_slice())
}

fn validate_positive_finite_values(context: &str, values: &[f64]) -> Result<(), DeseqError> {
    for (idx, value) in values.iter().copied().enumerate() {
        if !value.is_finite() || value <= 0.0 {
            return Err(DeseqError::InvalidSizeFactors {
                reason: format!("{context} value at index {idx} must be finite and positive"),
            });
        }
    }
    Ok(())
}

fn validate_nonnegative_finite_values(context: &str, values: &[f64]) -> Result<(), DeseqError> {
    for (idx, value) in values.iter().copied().enumerate() {
        if !value.is_finite() || value < 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: context.to_string(),
                index: Some(idx),
                value,
            });
        }
    }
    Ok(())
}

fn format_finite_or_na(value: f64) -> String {
    if value.is_finite() {
        value.to_string()
    } else {
        "NA".to_string()
    }
}

fn format_optional_finite_or_na(value: Option<f64>) -> String {
    value
        .map(format_finite_or_na)
        .unwrap_or_else(|| "NA".to_string())
}

fn format_r_logical(value: bool) -> String {
    if value {
        "TRUE".to_string()
    } else {
        "FALSE".to_string()
    }
}

fn format_optional_r_logical(value: Option<bool>) -> String {
    value
        .map(format_r_logical)
        .unwrap_or_else(|| "NA".to_string())
}

/// Write an assembled DESeq2-shaped result table to a tab-delimited file.
///
/// Gene identifiers are written as the first `gene` column, while statistical
/// columns follow [`DeseqResults::column_names`]. Missing numeric or logical
/// values are written as `NA`, matching R-style result-table exports.
pub fn write_deseq_results_tsv(
    path: impl AsRef<Path>,
    results: &DeseqResults,
) -> Result<(), DeseqError> {
    write_deseq_results_with_row_header_tsv(path, results, "gene", "gene")
}

/// Write a DESeq2 `results(tidy = TRUE)`-style table to a tab-delimited file.
///
/// Row identifiers are written as the first `row` column, matching the tidy
/// result-table shape from DESeq2. Statistical columns follow
/// [`DeseqResults::column_names`], and missing values are written as `NA`.
pub fn write_deseq_results_tidy_tsv(
    path: impl AsRef<Path>,
    results: &DeseqResults,
) -> Result<(), DeseqError> {
    write_deseq_results_with_row_header_tsv(path, results, "row", "row")
}

fn write_deseq_results_with_row_header_tsv(
    path: impl AsRef<Path>,
    results: &DeseqResults,
    row_header: &str,
    fallback_prefix: &str,
) -> Result<(), DeseqError> {
    validate_deseq_results_export(results)?;
    let frame = results.data_frame();
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    let mut header = Vec::with_capacity(frame.columns.len() + 1);
    header.push(row_header.to_string());
    header.extend(
        frame
            .columns
            .iter()
            .map(|column| column.metadata.name.clone()),
    );
    writer.write_record(&header)?;

    for row_idx in 0..frame.row_names.len() {
        let mut record = Vec::with_capacity(frame.columns.len() + 1);
        record.push(
            frame.row_names[row_idx]
                .clone()
                .unwrap_or_else(|| format!("{}{}", fallback_prefix, row_idx + 1)),
        );
        for column in &frame.columns {
            record.push(format_result_column_value(&column.values, row_idx));
        }
        writer.write_record(&record)?;
    }
    writer.flush()?;
    Ok(())
}

fn validate_deseq_results_export(results: &DeseqResults) -> Result<(), DeseqError> {
    for (idx, row) in results.rows.iter().enumerate() {
        if !row.base_mean.is_finite() || row.base_mean < 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: "result export baseMean".to_string(),
                index: Some(idx),
                value: row.base_mean,
            });
        }
        validate_optional_export_finite(row.log2_fold_change, "result export log2FoldChange", idx)?;
        validate_optional_export_finite(row.lfc_se, "result export lfcSE", idx)?;
        validate_optional_export_finite(row.stat, "result export stat", idx)?;
        validate_optional_export_probability(row.pvalue, "result export pvalue", idx)?;
        validate_optional_export_probability(row.padj, "result export padj", idx)?;
        validate_optional_export_positive(row.dispersion, "result export dispersion", idx)?;
        validate_optional_export_nonnegative(row.max_cooks, "result export maxCooks", idx)?;
    }
    Ok(())
}

fn validate_optional_export_finite(
    value: Option<f64>,
    context: &str,
    idx: usize,
) -> Result<(), DeseqError> {
    if let Some(value) = value {
        if !value.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: context.to_string(),
                index: Some(idx),
                value,
            });
        }
    }
    Ok(())
}

fn validate_optional_export_probability(
    value: Option<f64>,
    context: &str,
    idx: usize,
) -> Result<(), DeseqError> {
    if let Some(value) = value {
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(DeseqError::InvalidOptions {
                reason: format!("{context} at index {idx} must be finite and within [0, 1]"),
            });
        }
    }
    Ok(())
}

fn validate_optional_export_positive(
    value: Option<f64>,
    context: &str,
    idx: usize,
) -> Result<(), DeseqError> {
    if let Some(value) = value {
        if !value.is_finite() || value <= 0.0 {
            return Err(DeseqError::InvalidDispersion {
                reason: format!("{context} at index {idx} must be finite and positive"),
            });
        }
    }
    Ok(())
}

fn validate_optional_export_nonnegative(
    value: Option<f64>,
    context: &str,
    idx: usize,
) -> Result<(), DeseqError> {
    if let Some(value) = value {
        if !value.is_finite() || value < 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: context.to_string(),
                index: Some(idx),
                value,
            });
        }
    }
    Ok(())
}

/// Write DESeq2-style result column metadata to a tab-delimited file.
///
/// The output mirrors the `type` and `description` columns available from
/// `mcols(res)`, with an added `name` column for the result-table column name.
pub fn write_deseq_result_column_metadata_tsv(
    path: impl AsRef<Path>,
    results: &DeseqResults,
) -> Result<(), DeseqError> {
    let metadata = results.column_metadata();
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["name", "type", "description"])?;
    for column in metadata {
        writer.write_record([column.name, column.column_type, column.description])?;
    }
    writer.flush()?;
    Ok(())
}

/// Write table-level result metadata to a tab-delimited key/value file.
pub fn write_deseq_result_table_metadata_tsv(
    path: impl AsRef<Path>,
    results: &DeseqResults,
) -> Result<(), DeseqError> {
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["name", "value"])?;
    for entry in results.metadata.scalar_metadata() {
        writer.write_record([entry.name, entry.value])?;
    }
    writer.flush()?;
    Ok(())
}

/// Write Cook's replacement/refit scalar metadata to a tab-delimited key/value file.
pub fn write_cooks_replacement_metadata_tsv(
    path: impl AsRef<Path>,
    refit_plan: &CooksRefitPlan,
) -> Result<(), DeseqError> {
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["name", "value"])?;
    for entry in refit_plan.scalar_metadata() {
        writer.write_record([entry.name, entry.value])?;
    }
    writer.flush()?;
    Ok(())
}

/// Write DESeq2-shaped fit diagnostics to a tab-delimited file.
///
/// The columns follow [`Deseq2McolsDiagnostics::present_column_names`].
/// Missing optional numeric values are written as `NA`; logical values use
/// R-style `TRUE`/`FALSE` strings.
pub fn write_deseq_mcols_diagnostics_tsv(
    path: impl AsRef<Path>,
    gene_names: Option<&[String]>,
    diagnostics: &Deseq2McolsDiagnostics,
) -> Result<(), DeseqError> {
    let frame = diagnostics.data_frame();
    let n_rows = diagnostic_frame_row_count(&frame.columns)?;
    if let Some(names) = gene_names {
        if names.len() != n_rows {
            return Err(invalid_dimensions(
                "diagnostic gene names",
                n_rows,
                names.len(),
            ));
        }
    }
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    let mut header = Vec::with_capacity(frame.columns.len() + 1);
    header.push("gene".to_string());
    header.extend(frame.columns.iter().map(|column| column.name.to_string()));
    writer.write_record(&header)?;

    for row_idx in 0..n_rows {
        let mut record = Vec::with_capacity(frame.columns.len() + 1);
        record.push(
            gene_names
                .and_then(|names| names.get(row_idx))
                .cloned()
                .unwrap_or_else(|| format!("gene{}", row_idx + 1)),
        );
        for column in &frame.columns {
            record.push(format_diagnostic_column_value(&column.values, row_idx));
        }
        writer.write_record(&record)?;
    }
    writer.flush()?;
    Ok(())
}

/// Write DESeq2-style independent-filtering rejection counts.
///
/// The output mirrors `metadata(res)$filterNumRej`, with `theta` and `numRej`
/// columns.
pub fn write_independent_filter_num_rej_tsv(
    path: impl AsRef<Path>,
    filtering: &IndependentFilteringOutput,
) -> Result<(), DeseqError> {
    validate_independent_filtering_export(filtering)?;
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["theta", "numRej"])?;
    for row in filtering.filter_num_rej() {
        writer.write_record([row.theta.to_string(), row.num_rejections.to_string()])?;
    }
    writer.flush()?;
    Ok(())
}

/// Write DESeq2-style independent-filtering lowess metadata.
///
/// The output mirrors `metadata(res)$lo.fit`, with `x` and `y` columns.
pub fn write_independent_filter_lowess_tsv(
    path: impl AsRef<Path>,
    filtering: &IndependentFilteringOutput,
) -> Result<(), DeseqError> {
    validate_independent_filtering_export(filtering)?;
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["x", "y"])?;
    for row in filtering.lowess_fit_table() {
        writer.write_record([row.theta.to_string(), row.fitted_rejections.to_string()])?;
    }
    writer.flush()?;
    Ok(())
}

/// Write scalar DESeq2-style independent-filtering metadata.
///
/// The output contains `name` and `value` columns for scalar entries such as
/// `filterThreshold`, `filterTheta`, and `alpha`.
pub fn write_independent_filter_metadata_tsv(
    path: impl AsRef<Path>,
    filtering: &IndependentFilteringOutput,
) -> Result<(), DeseqError> {
    validate_independent_filtering_export(filtering)?;
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["name", "value"])?;
    for entry in filtering.scalar_metadata() {
        writer.write_record([entry.name, entry.value.to_string()])?;
    }
    writer.flush()?;
    Ok(())
}

fn validate_independent_filtering_export(
    filtering: &IndependentFilteringOutput,
) -> Result<(), DeseqError> {
    if filtering.theta.len() != filtering.num_rejections.len() {
        return Err(invalid_dimensions(
            "independent-filter numRej rows",
            filtering.theta.len(),
            filtering.num_rejections.len(),
        ));
    }
    for (idx, theta) in filtering.theta.iter().copied().enumerate() {
        if !theta.is_finite() || !(0.0..=1.0).contains(&theta) {
            return Err(DeseqError::InvalidOptions {
                reason: format!("independent-filter theta at index {idx} must be within [0, 1]"),
            });
        }
    }
    if let Some(lowess_fit) = &filtering.lowess_fit {
        if lowess_fit.len() != filtering.theta.len() {
            return Err(invalid_dimensions(
                "independent-filter lowess rows",
                filtering.theta.len(),
                lowess_fit.len(),
            ));
        }
        validate_nonnegative_finite_values("independent-filter lowess export", lowess_fit)?;
    }
    if let Some(selected_index) = filtering.selected_index {
        if selected_index >= filtering.theta.len() {
            return Err(invalid_dimensions(
                "independent-filter selected index",
                filtering.theta.len().saturating_sub(1),
                selected_index,
            ));
        }
    }
    validate_optional_export_probability(filtering.filter_theta, "independent-filter theta", 0)?;
    validate_optional_export_nonnegative(
        filtering.filter_threshold,
        "independent-filter threshold",
        0,
    )?;
    validate_optional_export_probability(Some(filtering.alpha), "independent-filter alpha", 0)?;
    Ok(())
}

fn diagnostic_frame_row_count(
    columns: &[Deseq2McolsDiagnosticColumn],
) -> Result<usize, DeseqError> {
    let Some(first) = columns.first() else {
        return Ok(0);
    };
    let expected = diagnostic_column_len(first);
    for column in columns.iter().skip(1) {
        let actual = diagnostic_column_len(column);
        if actual != expected {
            return Err(invalid_dimensions(
                format!("diagnostic column {}", column.name),
                expected,
                actual,
            ));
        }
    }
    Ok(expected)
}

fn diagnostic_column_len(column: &Deseq2McolsDiagnosticColumn) -> usize {
    match &column.values {
        Deseq2McolsDiagnosticValues::Numeric(values) => values.len(),
        Deseq2McolsDiagnosticValues::OptionalNumeric(values) => values.len(),
        Deseq2McolsDiagnosticValues::Integer(values) => values.len(),
        Deseq2McolsDiagnosticValues::Logical(values) => values.len(),
    }
}

fn format_result_column_value(values: &DeseqResultColumnValues, row_idx: usize) -> String {
    match values {
        DeseqResultColumnValues::Numeric(values) => values
            .get(row_idx)
            .copied()
            .flatten()
            .filter(|value| value.is_finite())
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NA".to_string()),
        DeseqResultColumnValues::Logical(values) => values
            .get(row_idx)
            .copied()
            .flatten()
            .map(|value| {
                if value {
                    "TRUE".to_string()
                } else {
                    "FALSE".to_string()
                }
            })
            .unwrap_or_else(|| "NA".to_string()),
    }
}

fn format_diagnostic_column_value(values: &Deseq2McolsDiagnosticValues, row_idx: usize) -> String {
    match values {
        Deseq2McolsDiagnosticValues::Numeric(values) => values
            .get(row_idx)
            .copied()
            .filter(|value| value.is_finite())
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NA".to_string()),
        Deseq2McolsDiagnosticValues::OptionalNumeric(values) => values
            .get(row_idx)
            .copied()
            .flatten()
            .filter(|value| value.is_finite())
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NA".to_string()),
        Deseq2McolsDiagnosticValues::Integer(values) => values
            .get(row_idx)
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NA".to_string()),
        Deseq2McolsDiagnosticValues::Logical(values) => values
            .get(row_idx)
            .map(|value| {
                if *value {
                    "TRUE".to_string()
                } else {
                    "FALSE".to_string()
                }
            })
            .unwrap_or_else(|| "NA".to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagnostics::Deseq2McolsDiagnostics;
    use crate::independent_filtering::IndependentFilteringOutput;
    use crate::results::{DeseqResultRow, DeseqResults};
    use std::fs;
    use std::path::PathBuf;

    #[test]
    fn parse_count_field_accepts_integer_scientific_notation() {
        let path = Path::new("counts.tsv");
        assert_eq!(parse_count_field("1e+05", path).unwrap(), 100_000);
        assert_eq!(parse_count_field("42.0", path).unwrap(), 42);
    }

    #[test]
    fn parse_count_field_rejects_non_integer_numeric_values() {
        let path = Path::new("counts.tsv");
        assert!(matches!(
            parse_count_field("1.5", path),
            Err(DeseqError::InvalidCounts { .. })
        ));
        assert!(matches!(
            parse_count_field("-1", path),
            Err(DeseqError::InvalidCounts { .. })
        ));
    }

    #[test]
    fn read_design_matrix_tsv_reads_numeric_matrix_and_names() {
        let path = unique_test_path("design.tsv");
        fs::write(
            &path,
            concat!(
                "sample\tIntercept\tcondition_B_vs_A\n",
                "s1\t1\t0\n",
                "s2\t1\t0\n",
                "s3\t1\t1\n",
            ),
        )
        .unwrap();

        let design = read_design_matrix_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(design.n_samples(), 3);
        assert_eq!(design.n_coefficients(), 2);
        assert_eq!(
            design.coefficient_names().unwrap(),
            &["Intercept".to_string(), "condition_B_vs_A".to_string()]
        );
        assert_eq!(design.matrix().as_slice(), &[1.0, 0.0, 1.0, 0.0, 1.0, 1.0]);
    }

    #[test]
    fn read_labeled_design_matrix_tsv_preserves_sample_names() {
        let path = unique_test_path("labeled_design.tsv");
        fs::write(
            &path,
            concat!(
                "sample\tIntercept\tcondition_B_vs_A\n",
                "s1\t1\t0\n",
                "s2\t1\t0\n",
                "s3\t1\t1\n",
            ),
        )
        .unwrap();

        let labeled = read_labeled_design_matrix_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(
            labeled.sample_names,
            vec!["s1".to_string(), "s2".to_string(), "s3".to_string()]
        );
        assert_eq!(
            labeled.design.matrix().as_slice(),
            &[1.0, 0.0, 1.0, 0.0, 1.0, 1.0]
        );
    }

    #[test]
    fn align_design_matrix_to_samples_uses_count_sample_order() {
        let design = DesignMatrix::from_row_major(
            3,
            2,
            vec![1.0, 1.0, 1.0, 0.0, 1.0, 2.0],
            Some(vec!["Intercept".to_string(), "condition".to_string()]),
        )
        .unwrap();
        let labeled = LabeledDesignMatrix {
            design,
            sample_names: vec!["s3".to_string(), "s1".to_string(), "s2".to_string()],
        };
        let samples = vec!["s1".to_string(), "s2".to_string(), "s3".to_string()];

        let aligned = align_design_matrix_to_samples(labeled, &samples).unwrap();

        assert_eq!(aligned.matrix().as_slice(), &[1.0, 0.0, 1.0, 2.0, 1.0, 1.0]);
        assert_eq!(
            aligned.coefficient_names().unwrap(),
            &["Intercept".to_string(), "condition".to_string()]
        );
    }

    #[test]
    fn align_design_matrix_to_samples_rejects_missing_and_duplicate_samples() {
        let samples = vec!["s1".to_string(), "s2".to_string()];
        let missing = LabeledDesignMatrix {
            design: DesignMatrix::from_row_major(2, 1, vec![1.0, 1.0], Some(vec!["x".to_string()]))
                .unwrap(),
            sample_names: vec!["s1".to_string(), "s3".to_string()],
        };
        assert!(align_design_matrix_to_samples(missing, &samples).is_err());

        let duplicated = LabeledDesignMatrix {
            design: DesignMatrix::from_row_major(2, 1, vec![1.0, 1.0], Some(vec!["x".to_string()]))
                .unwrap(),
            sample_names: vec!["s1".to_string(), "s1".to_string()],
        };
        assert!(align_design_matrix_to_samples(duplicated, &samples).is_err());
    }

    #[test]
    fn read_design_matrix_tsv_validates_shape_and_values() {
        let bad_value = unique_test_path("bad_design_value.tsv");
        fs::write(
            &bad_value,
            concat!("sample\tIntercept\tcondition_B_vs_A\n", "s1\t1\tNA\n",),
        )
        .unwrap();
        assert!(matches!(
            read_design_matrix_tsv(&bad_value),
            Err(DeseqError::ParseFloat { .. })
        ));
        let _ = fs::remove_file(&bad_value);

        let bad_shape = unique_test_path("bad_design_shape.tsv");
        fs::write(
            &bad_shape,
            concat!("sample\tIntercept\tcondition_B_vs_A\n", "s1\t1\n",),
        )
        .unwrap();
        assert!(read_design_matrix_tsv(&bad_shape).is_err());
        let _ = fs::remove_file(&bad_shape);
    }

    #[test]
    fn read_normalization_factors_tsv_reads_positive_matrix() {
        let path = unique_test_path("read_normalization_factors.tsv");
        fs::write(
            &path,
            concat!(
                "gene\tsample_1\tsample_2\n",
                "gene_a\t1\t2\n",
                "gene_b\t0.5\t4\n",
            ),
        )
        .unwrap();

        let factors = read_normalization_factors_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(factors.n_rows(), 2);
        assert_eq!(factors.n_cols(), 2);
        assert_eq!(factors.as_slice(), &[1.0, 2.0, 0.5, 4.0]);
    }

    #[test]
    fn read_labeled_normalization_factors_tsv_preserves_names() {
        let path = unique_test_path("read_labeled_normalization_factors.tsv");
        fs::write(
            &path,
            concat!(
                "gene\tsample_1\tsample_2\n",
                "gene_a\t1\t2\n",
                "gene_b\t0.5\t4\n",
            ),
        )
        .unwrap();

        let factors = read_labeled_normalization_factors_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(
            factors.gene_names,
            vec!["gene_a".to_string(), "gene_b".to_string()]
        );
        assert_eq!(
            factors.sample_names,
            vec!["sample_1".to_string(), "sample_2".to_string()]
        );
        assert_eq!(factors.matrix.as_slice(), &[1.0, 2.0, 0.5, 4.0]);
    }

    #[test]
    fn align_labeled_assay_matrix_to_counts_uses_gene_and_sample_order() {
        let counts = CountMatrix::from_row_major_u32_with_names(
            2,
            2,
            vec![1, 2, 3, 4],
            Some(vec!["gene_a".to_string(), "gene_b".to_string()]),
            Some(vec!["sample_1".to_string(), "sample_2".to_string()]),
        )
        .unwrap();
        let labeled = LabeledAssayMatrix {
            matrix: RowMajorMatrix::from_row_major(2, 2, vec![4.0, 3.0, 2.0, 1.0]).unwrap(),
            gene_names: vec!["gene_b".to_string(), "gene_a".to_string()],
            sample_names: vec!["sample_2".to_string(), "sample_1".to_string()],
        };

        let aligned =
            align_labeled_assay_matrix_to_counts(labeled, &counts, "test matrix").unwrap();

        assert_eq!(aligned.as_slice(), &[1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn align_labeled_assay_matrix_to_counts_rejects_missing_and_duplicate_names() {
        let counts = CountMatrix::from_row_major_u32_with_names(
            1,
            2,
            vec![1, 2],
            Some(vec!["gene_a".to_string()]),
            Some(vec!["sample_1".to_string(), "sample_2".to_string()]),
        )
        .unwrap();
        let missing = LabeledAssayMatrix {
            matrix: RowMajorMatrix::from_row_major(1, 2, vec![1.0, 2.0]).unwrap(),
            gene_names: vec!["gene_b".to_string()],
            sample_names: vec!["sample_1".to_string(), "sample_2".to_string()],
        };
        assert!(align_labeled_assay_matrix_to_counts(missing, &counts, "test matrix").is_err());

        let duplicated = LabeledAssayMatrix {
            matrix: RowMajorMatrix::from_row_major(1, 2, vec![1.0, 2.0]).unwrap(),
            gene_names: vec!["gene_a".to_string()],
            sample_names: vec!["sample_1".to_string(), "sample_1".to_string()],
        };
        assert!(align_labeled_assay_matrix_to_counts(duplicated, &counts, "test matrix").is_err());
    }

    #[test]
    fn read_normalization_factors_tsv_validates_shape_and_values() {
        let bad_value = unique_test_path("bad_normalization_factor_value.tsv");
        fs::write(
            &bad_value,
            concat!("gene\tsample_1\tsample_2\n", "gene_a\t1\t0\n",),
        )
        .unwrap();
        assert!(matches!(
            read_normalization_factors_tsv(&bad_value),
            Err(DeseqError::InvalidSizeFactors { .. })
        ));
        let _ = fs::remove_file(&bad_value);

        let bad_shape = unique_test_path("bad_normalization_factor_shape.tsv");
        fs::write(
            &bad_shape,
            concat!("gene\tsample_1\tsample_2\n", "gene_a\t1\n",),
        )
        .unwrap();
        assert!(read_normalization_factors_tsv(&bad_shape).is_err());
        let _ = fs::remove_file(&bad_shape);
    }

    #[test]
    fn read_observation_weights_tsv_reads_nonnegative_matrix() {
        let path = unique_test_path("read_observation_weights.tsv");
        fs::write(
            &path,
            concat!(
                "gene\tsample_1\tsample_2\n",
                "gene_a\t1\t0\n",
                "gene_b\t0.5\t4\n",
            ),
        )
        .unwrap();

        let weights = read_observation_weights_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(weights.n_rows(), 2);
        assert_eq!(weights.n_cols(), 2);
        assert_eq!(weights.as_slice(), &[1.0, 0.0, 0.5, 4.0]);
    }

    #[test]
    fn read_labeled_observation_weights_tsv_preserves_names() {
        let path = unique_test_path("read_labeled_observation_weights.tsv");
        fs::write(
            &path,
            concat!(
                "gene\tsample_1\tsample_2\n",
                "gene_a\t1\t0\n",
                "gene_b\t0.5\t4\n",
            ),
        )
        .unwrap();

        let weights = read_labeled_observation_weights_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(
            weights.gene_names,
            vec!["gene_a".to_string(), "gene_b".to_string()]
        );
        assert_eq!(
            weights.sample_names,
            vec!["sample_1".to_string(), "sample_2".to_string()]
        );
        assert_eq!(weights.matrix.as_slice(), &[1.0, 0.0, 0.5, 4.0]);
    }

    #[test]
    fn read_observation_weights_tsv_validates_shape_and_values() {
        let bad_value = unique_test_path("bad_observation_weight_value.tsv");
        fs::write(
            &bad_value,
            concat!("gene\tsample_1\tsample_2\n", "gene_a\t1\t-0.1\n",),
        )
        .unwrap();
        assert!(matches!(
            read_observation_weights_tsv(&bad_value),
            Err(DeseqError::NonFiniteValue { .. })
        ));
        let _ = fs::remove_file(&bad_value);

        let bad_shape = unique_test_path("bad_observation_weight_shape.tsv");
        fs::write(
            &bad_shape,
            concat!("gene\tsample_1\tsample_2\n", "gene_a\t1\n",),
        )
        .unwrap();
        assert!(read_observation_weights_tsv(&bad_shape).is_err());
        let _ = fs::remove_file(&bad_shape);
    }

    #[test]
    fn read_size_factors_tsv_reads_positive_values() {
        let path = unique_test_path("read_size_factors.tsv");
        fs::write(
            &path,
            concat!(
                "sample\tsize_factor\n",
                "sample_1\t1\n",
                "sample_2\t0.5\n",
                "sample_3\t2\n",
            ),
        )
        .unwrap();

        let size_factors = read_size_factors_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(size_factors, vec![1.0, 0.5, 2.0]);
    }

    #[test]
    fn read_labeled_size_factors_tsv_preserves_sample_names() {
        let path = unique_test_path("read_labeled_size_factors.tsv");
        fs::write(
            &path,
            concat!("sample\tsize_factor\n", "sample_1\t1\n", "sample_2\t0.5\n",),
        )
        .unwrap();

        let size_factors = read_labeled_size_factors_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(
            size_factors,
            vec![
                SampleNumericValue {
                    sample: "sample_1".to_string(),
                    value: 1.0
                },
                SampleNumericValue {
                    sample: "sample_2".to_string(),
                    value: 0.5
                },
            ]
        );
    }

    #[test]
    fn align_sample_numeric_values_to_samples_uses_count_sample_order() {
        let values = vec![
            SampleNumericValue {
                sample: "sample_3".to_string(),
                value: 3.0,
            },
            SampleNumericValue {
                sample: "sample_1".to_string(),
                value: 1.0,
            },
            SampleNumericValue {
                sample: "sample_2".to_string(),
                value: 2.0,
            },
        ];
        let samples = vec![
            "sample_1".to_string(),
            "sample_2".to_string(),
            "sample_3".to_string(),
        ];

        assert_eq!(
            align_sample_numeric_values_to_samples(&values, &samples, "size-factor").unwrap(),
            vec![1.0, 2.0, 3.0]
        );
    }

    #[test]
    fn align_sample_numeric_values_to_samples_rejects_missing_and_duplicate_samples() {
        let samples = vec!["sample_1".to_string(), "sample_2".to_string()];
        let missing = vec![
            SampleNumericValue {
                sample: "sample_1".to_string(),
                value: 1.0,
            },
            SampleNumericValue {
                sample: "sample_3".to_string(),
                value: 3.0,
            },
        ];
        assert!(align_sample_numeric_values_to_samples(&missing, &samples, "size-factor").is_err());

        let duplicated = vec![
            SampleNumericValue {
                sample: "sample_1".to_string(),
                value: 1.0,
            },
            SampleNumericValue {
                sample: "sample_1".to_string(),
                value: 2.0,
            },
        ];
        assert!(
            align_sample_numeric_values_to_samples(&duplicated, &samples, "size-factor").is_err()
        );
    }

    #[test]
    fn read_size_factors_tsv_validates_shape_and_values() {
        let bad_value = unique_test_path("bad_size_factor_value.tsv");
        fs::write(
            &bad_value,
            concat!("sample\tsize_factor\n", "sample_1\t0\n",),
        )
        .unwrap();
        assert!(matches!(
            read_size_factors_tsv(&bad_value),
            Err(DeseqError::InvalidSizeFactors { .. })
        ));
        let _ = fs::remove_file(&bad_value);

        let bad_shape = unique_test_path("bad_size_factor_shape.tsv");
        fs::write(
            &bad_shape,
            concat!("sample\tsize_factor\n", "sample_1\t1\textra\n",),
        )
        .unwrap();
        assert!(read_size_factors_tsv(&bad_shape).is_err());
        let _ = fs::remove_file(&bad_shape);
    }

    #[test]
    fn read_geometric_means_tsv_reads_nonnegative_values() {
        let path = unique_test_path("read_geometric_means.tsv");
        fs::write(
            &path,
            concat!(
                "gene\tgeo_mean\n",
                "gene_1\t1\n",
                "gene_2\t0\n",
                "gene_3\t2.5\n",
            ),
        )
        .unwrap();

        let geometric_means = read_geometric_means_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(geometric_means, vec![1.0, 0.0, 2.5]);
    }

    #[test]
    fn read_labeled_geometric_means_tsv_preserves_gene_names() {
        let path = unique_test_path("read_labeled_geometric_means.tsv");
        fs::write(
            &path,
            concat!("gene\tgeo_mean\n", "gene_1\t1\n", "gene_2\t0\n",),
        )
        .unwrap();

        let geometric_means = read_labeled_geometric_means_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(
            geometric_means,
            vec![
                GeneNumericValue {
                    gene: "gene_1".to_string(),
                    value: 1.0
                },
                GeneNumericValue {
                    gene: "gene_2".to_string(),
                    value: 0.0
                },
            ]
        );
    }

    #[test]
    fn align_gene_numeric_values_to_genes_uses_count_gene_order() {
        let values = vec![
            GeneNumericValue {
                gene: "gene_3".to_string(),
                value: 3.0,
            },
            GeneNumericValue {
                gene: "gene_1".to_string(),
                value: 1.0,
            },
            GeneNumericValue {
                gene: "gene_2".to_string(),
                value: 2.0,
            },
        ];
        let genes = vec![
            "gene_1".to_string(),
            "gene_2".to_string(),
            "gene_3".to_string(),
        ];

        assert_eq!(
            align_gene_numeric_values_to_genes(&values, &genes, "geometric-mean").unwrap(),
            vec![1.0, 2.0, 3.0]
        );
    }

    #[test]
    fn align_gene_numeric_values_to_genes_rejects_missing_and_duplicate_genes() {
        let genes = vec!["gene_1".to_string(), "gene_2".to_string()];
        let missing = vec![
            GeneNumericValue {
                gene: "gene_1".to_string(),
                value: 1.0,
            },
            GeneNumericValue {
                gene: "gene_3".to_string(),
                value: 3.0,
            },
        ];
        assert!(align_gene_numeric_values_to_genes(&missing, &genes, "geometric-mean").is_err());

        let duplicated = vec![
            GeneNumericValue {
                gene: "gene_1".to_string(),
                value: 1.0,
            },
            GeneNumericValue {
                gene: "gene_1".to_string(),
                value: 2.0,
            },
        ];
        assert!(align_gene_numeric_values_to_genes(&duplicated, &genes, "geometric-mean").is_err());
    }

    #[test]
    fn read_geometric_means_tsv_validates_shape_and_values() {
        let bad_value = unique_test_path("bad_geometric_mean_value.tsv");
        fs::write(&bad_value, concat!("gene\tgeo_mean\n", "gene_1\t-1\n")).unwrap();
        assert!(matches!(
            read_geometric_means_tsv(&bad_value),
            Err(DeseqError::InvalidSizeFactors { .. })
        ));
        let _ = fs::remove_file(&bad_value);

        let bad_shape = unique_test_path("bad_geometric_mean_shape.tsv");
        fs::write(
            &bad_shape,
            concat!("gene\tgeo_mean\n", "gene_1\t1\textra\n"),
        )
        .unwrap();
        assert!(read_geometric_means_tsv(&bad_shape).is_err());
        let _ = fs::remove_file(&bad_shape);
    }

    #[test]
    fn read_wald_t_degrees_of_freedom_tsv_reads_finite_values() {
        let path = unique_test_path("read_wald_t_df.tsv");
        fs::write(&path, concat!("gene\tdf\n", "gene_1\t4\n", "gene_2\t2.5\n")).unwrap();

        let degrees_of_freedom = read_wald_t_degrees_of_freedom_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(degrees_of_freedom, vec![4.0, 2.5]);
    }

    #[test]
    fn read_labeled_wald_t_degrees_of_freedom_tsv_preserves_gene_names() {
        let path = unique_test_path("read_labeled_wald_t_df.tsv");
        fs::write(&path, concat!("gene\tdf\n", "gene_1\t4\n", "gene_2\t2.5\n")).unwrap();

        let degrees_of_freedom = read_labeled_wald_t_degrees_of_freedom_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(
            degrees_of_freedom,
            vec![
                GeneNumericValue {
                    gene: "gene_1".to_string(),
                    value: 4.0
                },
                GeneNumericValue {
                    gene: "gene_2".to_string(),
                    value: 2.5
                },
            ]
        );
    }

    #[test]
    fn read_wald_t_degrees_of_freedom_tsv_validates_shape_and_values() {
        let bad_value = unique_test_path("bad_wald_t_df_value.tsv");
        fs::write(&bad_value, concat!("gene\tdf\n", "gene_1\tNA\n")).unwrap();
        assert!(read_wald_t_degrees_of_freedom_tsv(&bad_value).is_err());
        let _ = fs::remove_file(&bad_value);

        let bad_shape = unique_test_path("bad_wald_t_df_shape.tsv");
        fs::write(&bad_shape, concat!("gene\tdf\n", "gene_1\t4\textra\n")).unwrap();
        assert!(read_wald_t_degrees_of_freedom_tsv(&bad_shape).is_err());
        let _ = fs::remove_file(&bad_shape);
    }

    #[test]
    fn read_sample_levels_tsv_reads_labeled_string_levels() {
        let path = unique_test_path("read_sample_levels.tsv");
        fs::write(
            &path,
            concat!(
                "sample\tcondition\n",
                "sample_1\tA\n",
                "sample_2\tB\n",
                "sample_3\tA\n",
            ),
        )
        .unwrap();

        let levels = read_sample_levels_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(
            levels,
            vec![
                SampleLevel {
                    sample: "sample_1".to_string(),
                    level: "A".to_string()
                },
                SampleLevel {
                    sample: "sample_2".to_string(),
                    level: "B".to_string()
                },
                SampleLevel {
                    sample: "sample_3".to_string(),
                    level: "A".to_string()
                },
            ]
        );
    }

    #[test]
    fn read_sample_levels_tsv_validates_shape_and_values() {
        let bad_value = unique_test_path("bad_sample_level_value.tsv");
        fs::write(&bad_value, concat!("sample\tcondition\n", "sample_1\t\n")).unwrap();
        assert!(read_sample_levels_tsv(&bad_value).is_err());
        let _ = fs::remove_file(&bad_value);

        let bad_shape = unique_test_path("bad_sample_level_shape.tsv");
        fs::write(
            &bad_shape,
            concat!("sample\tcondition\n", "sample_1\tA\textra\n"),
        )
        .unwrap();
        assert!(read_sample_levels_tsv(&bad_shape).is_err());
        let _ = fs::remove_file(&bad_shape);
    }

    #[test]
    fn align_sample_levels_to_samples_uses_count_sample_order() {
        let levels = vec![
            SampleLevel {
                sample: "sample_3".to_string(),
                level: "C".to_string(),
            },
            SampleLevel {
                sample: "sample_1".to_string(),
                level: "A".to_string(),
            },
            SampleLevel {
                sample: "sample_2".to_string(),
                level: "B".to_string(),
            },
        ];
        let samples = vec![
            "sample_1".to_string(),
            "sample_2".to_string(),
            "sample_3".to_string(),
        ];

        assert_eq!(
            align_sample_levels_to_samples(&levels, &samples).unwrap(),
            vec!["A".to_string(), "B".to_string(), "C".to_string()]
        );
    }

    #[test]
    fn align_sample_levels_to_samples_rejects_missing_and_duplicate_samples() {
        let samples = vec!["sample_1".to_string(), "sample_2".to_string()];
        let missing = vec![
            SampleLevel {
                sample: "sample_1".to_string(),
                level: "A".to_string(),
            },
            SampleLevel {
                sample: "sample_3".to_string(),
                level: "C".to_string(),
            },
        ];
        assert!(align_sample_levels_to_samples(&missing, &samples).is_err());

        let duplicated = vec![
            SampleLevel {
                sample: "sample_1".to_string(),
                level: "A".to_string(),
            },
            SampleLevel {
                sample: "sample_1".to_string(),
                level: "B".to_string(),
            },
        ];
        assert!(align_sample_levels_to_samples(&duplicated, &samples).is_err());
    }

    #[test]
    fn write_count_matrices_tsv_preserves_names_and_fallbacks() {
        let raw_path = unique_test_path("counts.tsv");
        let normalized_path = unique_test_path("normalized_counts.tsv");
        let unnamed_path = unique_test_path("unnamed_normalized_counts.tsv");
        let counts = CountMatrix::from_row_major_u32_with_names(
            2,
            3,
            vec![10, 20, 30, 5, 10, 20],
            Some(vec!["gene_a".to_string(), "gene_b".to_string()]),
            Some(vec![
                "sample_1".to_string(),
                "sample_2".to_string(),
                "sample_3".to_string(),
            ]),
        )
        .unwrap();
        let normalized =
            RowMajorMatrix::from_row_major(2, 3, vec![10.0, 10.0, 6.0, 5.0, 5.0, 4.0]).unwrap();

        write_count_matrix_tsv(&raw_path, &counts).unwrap();
        write_normalized_counts_tsv(
            &normalized_path,
            counts.gene_names(),
            counts.sample_names(),
            &normalized,
        )
        .unwrap();
        write_normalized_counts_tsv(&unnamed_path, None, None, &normalized).unwrap();

        let raw = fs::read_to_string(&raw_path).unwrap();
        let normalized_text = fs::read_to_string(&normalized_path).unwrap();
        let unnamed = fs::read_to_string(&unnamed_path).unwrap();
        let _ = fs::remove_file(&raw_path);
        let _ = fs::remove_file(&normalized_path);
        let _ = fs::remove_file(&unnamed_path);

        assert_eq!(
            raw,
            concat!(
                "gene\tsample_1\tsample_2\tsample_3\n",
                "gene_a\t10\t20\t30\n",
                "gene_b\t5\t10\t20\n",
            )
        );
        assert_eq!(
            normalized_text,
            concat!(
                "gene\tsample_1\tsample_2\tsample_3\n",
                "gene_a\t10\t10\t6\n",
                "gene_b\t5\t5\t4\n",
            )
        );
        assert_eq!(
            unnamed,
            concat!(
                "gene\tsample1\tsample2\tsample3\n",
                "gene1\t10\t10\t6\n",
                "gene2\t5\t5\t4\n",
            )
        );
    }

    #[test]
    fn write_normalized_counts_tsv_validates_names_and_finite_values() {
        let path = unique_test_path("bad_normalized_counts.tsv");
        let normalized = RowMajorMatrix::from_row_major(1, 2, vec![1.0, f64::INFINITY]).unwrap();

        assert!(matches!(
            write_normalized_counts_tsv(&path, Some(&["gene_a".to_string()]), None, &normalized),
            Err(DeseqError::NonFiniteValue { .. })
        ));
        assert!(matches!(
            write_normalized_counts_tsv(
                &path,
                Some(&["gene_a".to_string(), "gene_b".to_string()]),
                None,
                &RowMajorMatrix::from_row_major(1, 2, vec![1.0, 2.0]).unwrap(),
            ),
            Err(DeseqError::InvalidDimensions { .. })
        ));
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn write_normalization_factors_tsv_writes_deseq2_shape_and_validates_values() {
        let path = unique_test_path("normalization_factors.tsv");
        let bad_path = unique_test_path("bad_normalization_factors.tsv");
        let genes = vec!["gene_a".to_string(), "gene_b".to_string()];
        let samples = vec!["sample_1".to_string(), "sample_2".to_string()];
        let factors = RowMajorMatrix::from_row_major(2, 2, vec![1.0, 2.0, 0.5, 4.0]).unwrap();

        write_normalization_factors_tsv(&path, Some(&genes), Some(&samples), &factors).unwrap();

        let text = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(&path);
        assert_eq!(
            text,
            concat!(
                "gene\tsample_1\tsample_2\n",
                "gene_a\t1\t2\n",
                "gene_b\t0.5\t4\n",
            )
        );
        assert!(matches!(
            write_normalization_factors_tsv(
                &bad_path,
                None,
                None,
                &RowMajorMatrix::from_row_major(1, 2, vec![1.0, 0.0]).unwrap(),
            ),
            Err(DeseqError::InvalidSizeFactors { .. })
        ));
        assert!(matches!(
            write_normalization_factors_tsv(
                &bad_path,
                Some(&genes),
                None,
                &RowMajorMatrix::from_row_major(1, 2, vec![1.0, 2.0]).unwrap(),
            ),
            Err(DeseqError::InvalidDimensions { .. })
        ));
        let _ = fs::remove_file(&bad_path);
    }

    #[test]
    fn write_size_factors_tsv_validates_names_and_values() {
        let path = unique_test_path("bad_size_factors.tsv");
        let sample_names = vec!["sample_1".to_string(), "sample_2".to_string()];

        assert!(matches!(
            write_size_factors_tsv(&path, Some(&sample_names), &[1.0]),
            Err(DeseqError::InvalidDimensions { .. })
        ));
        assert!(matches!(
            write_size_factors_tsv(&path, None, &[1.0, 0.0]),
            Err(DeseqError::InvalidSizeFactors { .. })
        ));
        assert!(matches!(
            write_size_factors_tsv(&path, None, &[1.0, f64::NAN]),
            Err(DeseqError::InvalidSizeFactors { .. })
        ));
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn write_base_mean_tsv_validates_names_and_values() {
        let path = unique_test_path("bad_base_mean.tsv");
        let gene_names = vec!["gene_a".to_string(), "gene_b".to_string()];

        assert!(matches!(
            write_base_mean_tsv(&path, Some(&gene_names), &[1.0]),
            Err(DeseqError::InvalidDimensions { .. })
        ));
        assert!(matches!(
            write_base_mean_tsv(&path, None, &[1.0, -1.0]),
            Err(DeseqError::NonFiniteValue { .. })
        ));
        assert!(matches!(
            write_base_mean_tsv(&path, None, &[1.0, f64::INFINITY]),
            Err(DeseqError::NonFiniteValue { .. })
        ));
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn write_base_metadata_tsv_writes_deseq2_shape_and_validates_lengths() {
        let path = unique_test_path("base_metadata.tsv");
        let bad_path = unique_test_path("bad_base_metadata.tsv");
        let genes = vec!["gene_a".to_string(), "gene_b".to_string()];

        write_base_metadata_tsv(
            &path,
            Some(&genes),
            &[10.0, f64::NAN],
            &[2.5, f64::INFINITY],
            &[false, true],
        )
        .unwrap();

        let text = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(&path);
        assert_eq!(
            text,
            concat!(
                "gene\tbaseMean\tbaseVar\tallZero\n",
                "gene_a\t10\t2.5\tFALSE\n",
                "gene_b\tNA\tNA\tTRUE\n",
            )
        );
        assert!(matches!(
            write_base_metadata_tsv(&bad_path, Some(&genes), &[1.0], &[2.0], &[false]),
            Err(DeseqError::InvalidDimensions { .. })
        ));
        assert!(matches!(
            write_base_metadata_tsv(&bad_path, None, &[1.0, 2.0], &[3.0], &[false, true]),
            Err(DeseqError::InvalidDimensions { .. })
        ));
        let _ = fs::remove_file(&bad_path);
    }

    #[test]
    fn write_deseq_results_tsv_writes_columns_missing_values_and_logicals() {
        let path = unique_test_path("deseq_results.tsv");
        let tidy_path = unique_test_path("deseq_results_tidy.tsv");
        let metadata_path = unique_test_path("deseq_result_column_metadata.tsv");
        let table_metadata_path = unique_test_path("deseq_result_table_metadata.tsv");
        let results = DeseqResults {
            rows: vec![
                DeseqResultRow {
                    gene: Some("gene_a".to_string()),
                    base_mean: 10.0,
                    log2_fold_change: Some(1.25),
                    lfc_se: None,
                    stat: Some(-2.0),
                    pvalue: None,
                    padj: Some(0.05),
                    dispersion: Some(0.1),
                    converged: Some(true),
                    max_cooks: Some(4.0),
                    cooks_outlier: Some(false),
                    filtered: None,
                },
                DeseqResultRow {
                    gene: None,
                    base_mean: 20.0,
                    log2_fold_change: None,
                    lfc_se: Some(0.5),
                    stat: None,
                    pvalue: Some(0.8),
                    padj: None,
                    dispersion: None,
                    converged: Some(false),
                    max_cooks: None,
                    cooks_outlier: None,
                    filtered: None,
                },
            ],
            metadata: crate::results::DeseqResultsTableMetadata {
                test_type: Some(crate::options::TestType::Wald),
                result_name: Some("condition_B_vs_A".to_string()),
                comparison: Some("coefficient condition_B_vs_A".to_string()),
                lfc_threshold: 1.5,
                alt_hypothesis: Some("greater".to_string()),
                ..crate::results::DeseqResultsTableMetadata::default()
            },
            ..DeseqResults::default()
        };

        write_deseq_results_tsv(&path, &results).unwrap();
        write_deseq_results_tidy_tsv(&tidy_path, &results).unwrap();
        write_deseq_result_column_metadata_tsv(&metadata_path, &results).unwrap();
        write_deseq_result_table_metadata_tsv(&table_metadata_path, &results).unwrap();

        let text = fs::read_to_string(&path).unwrap();
        let tidy_text = fs::read_to_string(&tidy_path).unwrap();
        let metadata = fs::read_to_string(&metadata_path).unwrap();
        let table_metadata = fs::read_to_string(&table_metadata_path).unwrap();
        let _ = fs::remove_file(&path);
        let _ = fs::remove_file(&tidy_path);
        let _ = fs::remove_file(&metadata_path);
        let _ = fs::remove_file(&table_metadata_path);
        assert_eq!(
            text,
            concat!(
                "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\t",
                "dispersion\tconverged\tmaxCooks\tcooksOutlier\n",
                "gene_a\t10\t1.25\tNA\t-2\tNA\t0.05\t0.1\tTRUE\t4\tFALSE\n",
                "gene2\t20\tNA\t0.5\tNA\t0.8\tNA\tNA\tFALSE\tNA\tNA\n"
            )
        );
        assert_eq!(
            tidy_text,
            concat!(
                "row\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\t",
                "dispersion\tconverged\tmaxCooks\tcooksOutlier\n",
                "gene_a\t10\t1.25\tNA\t-2\tNA\t0.05\t0.1\tTRUE\t4\tFALSE\n",
                "row2\t20\tNA\t0.5\tNA\t0.8\tNA\tNA\tFALSE\tNA\tNA\n"
            )
        );
        assert!(metadata.starts_with("name\ttype\tdescription\n"));
        assert!(metadata.contains("baseMean\tresults\tmean of normalized counts for all samples\n"));
        assert!(metadata.contains("converged\tdiagnostic\twhether beta fitting converged\n"));
        assert_eq!(
            table_metadata,
            concat!(
                "name\tvalue\n",
                "testType\tWald\n",
                "resultName\tcondition_B_vs_A\n",
                "comparison\tcoefficient condition_B_vs_A\n",
                "lfcThreshold\t1.5\n",
                "altHypothesis\tgreater\n",
                "pAdjustMethod\tBH\n",
            )
        );
    }

    #[test]
    fn write_deseq_results_tsv_rejects_invalid_numeric_rows() {
        let path = unique_test_path("bad_deseq_results.tsv");
        let mut results = DeseqResults {
            rows: vec![DeseqResultRow {
                gene: Some("gene_a".to_string()),
                base_mean: 10.0,
                log2_fold_change: Some(1.25),
                lfc_se: Some(0.5),
                stat: Some(-2.0),
                pvalue: Some(0.04),
                padj: Some(0.05),
                dispersion: Some(0.1),
                converged: Some(true),
                max_cooks: Some(4.0),
                cooks_outlier: Some(false),
                filtered: None,
            }],
            ..DeseqResults::default()
        };

        results.rows[0].pvalue = Some(1.2);
        assert!(write_deseq_results_tsv(&path, &results).is_err());
        results.rows[0].pvalue = Some(0.04);
        results.rows[0].max_cooks = Some(f64::NAN);
        assert!(write_deseq_results_tsv(&path, &results).is_err());
        results.rows[0].max_cooks = Some(4.0);
        results.rows[0].base_mean = -1.0;
        assert!(write_deseq_results_tsv(&path, &results).is_err());
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn write_deseq_mcols_diagnostics_tsv_writes_present_columns_and_na_values() {
        let path = unique_test_path("mcols_diagnostics.tsv");
        let gene_names = vec!["gene_a".to_string(), "gene_b".to_string()];
        let diagnostics = Deseq2McolsDiagnostics {
            disp_gene_est: Some(vec![f64::NAN, 0.2]),
            disp_gene_iter: Some(vec![0, 3]),
            disp_outlier: Some(vec![false, true]),
            max_cooks: Some(vec![None, Some(4.5)]),
            ..Deseq2McolsDiagnostics::default()
        };

        write_deseq_mcols_diagnostics_tsv(&path, Some(&gene_names), &diagnostics).unwrap();

        let text = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(&path);
        assert_eq!(
            text,
            concat!(
                "gene\tdispGeneEst\tdispGeneIter\tdispOutlier\tmaxCooks\n",
                "gene_a\tNA\t0\tFALSE\tNA\n",
                "gene_b\t0.2\t3\tTRUE\t4.5\n"
            )
        );
    }

    #[test]
    fn write_cooks_diagnostics_tsv_writes_matrix_and_metadata() {
        let cooks_path = unique_test_path("cooks_distance.tsv");
        let row_path = unique_test_path("cooks_row_metadata.tsv");
        let sample_path = unique_test_path("cooks_sample_metadata.tsv");
        let gene_names = vec!["gene_a".to_string(), "gene_b".to_string()];
        let sample_names = vec!["sample_1".to_string(), "sample_2".to_string()];
        let cooks = CooksOutput {
            cooks: RowMajorMatrix::from_row_major(2, 2, vec![0.5, f64::NAN, 2.0, 4.0]).unwrap(),
            max_cooks: vec![None, Some(4.0)],
            robust_dispersion: vec![0.04, 0.25],
            samples_for_cooks: vec![true, false],
        };

        write_cooks_distance_tsv(&cooks_path, Some(&gene_names), Some(&sample_names), &cooks)
            .unwrap();
        write_cooks_row_metadata_tsv(&row_path, Some(&gene_names), &cooks).unwrap();
        write_cooks_sample_metadata_tsv(&sample_path, Some(&sample_names), &cooks).unwrap();

        let cooks_text = fs::read_to_string(&cooks_path).unwrap();
        let row_text = fs::read_to_string(&row_path).unwrap();
        let sample_text = fs::read_to_string(&sample_path).unwrap();
        let _ = fs::remove_file(&cooks_path);
        let _ = fs::remove_file(&row_path);
        let _ = fs::remove_file(&sample_path);

        assert_eq!(
            cooks_text,
            concat!(
                "gene\tsample_1\tsample_2\n",
                "gene_a\t0.5\tNA\n",
                "gene_b\t2\t4\n",
            )
        );
        assert_eq!(
            row_text,
            concat!(
                "gene\tmaxCooks\tcooksRobustDispersion\n",
                "gene_a\tNA\t0.04\n",
                "gene_b\t4\t0.25\n",
            )
        );
        assert_eq!(
            sample_text,
            concat!(
                "sample\tsamplesForCooks\n",
                "sample_1\tTRUE\n",
                "sample_2\tFALSE\n",
            )
        );
    }

    #[test]
    fn write_cooks_diagnostics_tsv_rejects_misaligned_metadata() {
        let path = unique_test_path("bad_cooks_metadata.tsv");
        let mut cooks = CooksOutput {
            cooks: RowMajorMatrix::from_row_major(2, 2, vec![0.5, 1.0, 2.0, 4.0]).unwrap(),
            max_cooks: vec![Some(1.0)],
            robust_dispersion: vec![0.04, 0.25],
            samples_for_cooks: vec![true, false],
        };

        assert!(write_cooks_row_metadata_tsv(&path, None, &cooks).is_err());
        cooks.max_cooks = vec![Some(1.0), Some(4.0)];
        cooks.robust_dispersion = vec![0.04];
        assert!(write_cooks_row_metadata_tsv(&path, None, &cooks).is_err());
        cooks.robust_dispersion = vec![0.04, 0.25];
        cooks.samples_for_cooks = vec![true];
        assert!(write_cooks_sample_metadata_tsv(&path, None, &cooks).is_err());
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn write_cooks_replacement_metadata_tsv_writes_scalar_summary() {
        let path = unique_test_path("cooks_replacement_metadata.tsv");
        let replaced_path = unique_test_path("cooks_replaced_counts.tsv");
        let candidate_path = unique_test_path("cooks_candidate_counts.tsv");
        let outlier_path = unique_test_path("cooks_outlier_cells.tsv");
        let row_metadata_path = unique_test_path("cooks_replacement_row_metadata.tsv");
        let counts = CountMatrix::from_row_major_u32(
            2,
            4,
            vec![
                10, 30, 20, 50, //
                0, 0, 0, 0,
            ],
        )
        .unwrap();
        let size_factors = vec![1.0, 1.0, 1.0, 1.0];
        let normalized = crate::normalization::normalized_counts(&counts, &size_factors).unwrap();
        let cooks = RowMajorMatrix::from_row_major(
            2,
            4,
            vec![
                0.0, 9.0, 0.0, 0.0, //
                8.0, 0.0, 0.0, 0.0,
            ],
        )
        .unwrap();
        let design = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();
        let plan = crate::cooks::prepare_cooks_replacement_refit(
            &counts,
            &normalized,
            &size_factors,
            None,
            &cooks,
            &design,
            &crate::cooks::CooksReplacementOptions {
                trim: 0.25,
                cooks_cutoff: 5.0,
                min_replicates: 3,
                which_samples: None,
            },
        )
        .unwrap();

        write_cooks_replacement_metadata_tsv(&path, &plan).unwrap();
        write_cooks_replaced_counts_tsv(&replaced_path, &plan).unwrap();
        write_cooks_candidate_replacement_counts_tsv(&candidate_path, &plan).unwrap();
        write_cooks_outlier_cells_tsv(&outlier_path, &plan).unwrap();
        write_cooks_replacement_row_metadata_tsv(&row_metadata_path, &plan).unwrap();

        let text = fs::read_to_string(&path).unwrap();
        let replaced = fs::read_to_string(&replaced_path).unwrap();
        let candidate = fs::read_to_string(&candidate_path).unwrap();
        let outlier = fs::read_to_string(&outlier_path).unwrap();
        let row_metadata = fs::read_to_string(&row_metadata_path).unwrap();
        let _ = fs::remove_file(&path);
        let _ = fs::remove_file(&replaced_path);
        let _ = fs::remove_file(&candidate_path);
        let _ = fs::remove_file(&outlier_path);
        let _ = fs::remove_file(&row_metadata_path);
        assert_eq!(
            text,
            concat!(
                "name\tvalue\n",
                "nRefit\t2\n",
                "nRefitRows\t1\n",
                "nNewAllZero\t1\n",
                "nOutlierCells\t2\n",
                "nReplacedCells\t2\n",
                "nReplaceableSamples\t4\n",
                "shouldRefit\ttrue\n",
            )
        );
        assert_eq!(
            replaced,
            concat!(
                "gene\tsample1\tsample2\tsample3\tsample4\n",
                "gene1\t10\t25\t20\t50\n",
                "gene2\t0\t0\t0\t0\n",
            )
        );
        assert_eq!(
            candidate,
            concat!(
                "gene\tsample1\tsample2\tsample3\tsample4\n",
                "gene1\t25\t25\t25\t25\n",
                "gene2\t0\t0\t0\t0\n",
            )
        );
        assert_eq!(
            outlier,
            concat!(
                "gene\tsample1\tsample2\tsample3\tsample4\n",
                "gene1\tFALSE\tTRUE\tFALSE\tFALSE\n",
                "gene2\tTRUE\tFALSE\tFALSE\tFALSE\n",
            )
        );
        assert_eq!(
            row_metadata,
            concat!(
                "gene\treplace\trefitReplace\tnewAllZero\treplacedAllZero\t",
                "replacedBaseMean\treplacedBaseVar\tpostRefitMaxCooks\n",
                "gene1\tTRUE\tTRUE\tFALSE\tFALSE\t26.25\t289.58333333333337\tNA\n",
                "gene2\tTRUE\tFALSE\tTRUE\tTRUE\t0\t0\tNA\n",
            )
        );
    }

    #[test]
    fn write_deseq_mcols_diagnostics_tsv_rejects_misaligned_columns() {
        let path = unique_test_path("mcols_bad_diagnostics.tsv");
        let diagnostics = Deseq2McolsDiagnostics {
            disp_gene_est: Some(vec![0.1, 0.2]),
            disp_gene_iter: Some(vec![1]),
            ..Deseq2McolsDiagnostics::default()
        };

        assert!(write_deseq_mcols_diagnostics_tsv(&path, None, &diagnostics).is_err());
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn write_deseq_mcols_diagnostics_tsv_rejects_misaligned_gene_names() {
        let path = unique_test_path("mcols_bad_names.tsv");
        let gene_names = vec!["gene_a".to_string()];
        let diagnostics = Deseq2McolsDiagnostics {
            disp_gene_est: Some(vec![0.1, 0.2]),
            ..Deseq2McolsDiagnostics::default()
        };

        assert!(write_deseq_mcols_diagnostics_tsv(&path, Some(&gene_names), &diagnostics).is_err());
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn write_independent_filter_metadata_tsv_writes_deseq2_shapes() {
        let num_rej_path = unique_test_path("filter_num_rej.tsv");
        let lowess_path = unique_test_path("filter_lowess.tsv");
        let scalar_path = unique_test_path("filter_metadata.tsv");
        let filtering = IndependentFilteringOutput {
            enabled: true,
            theta: vec![0.0, 0.5],
            num_rejections: vec![1, 3],
            selected_index: Some(1),
            filter_theta: Some(0.5),
            filter_threshold: Some(10.0),
            lowess_fit: Some(vec![1.25, 2.75]),
            alpha: 0.1,
        };

        write_independent_filter_num_rej_tsv(&num_rej_path, &filtering).unwrap();
        write_independent_filter_lowess_tsv(&lowess_path, &filtering).unwrap();
        write_independent_filter_metadata_tsv(&scalar_path, &filtering).unwrap();

        let num_rej = fs::read_to_string(&num_rej_path).unwrap();
        let lowess = fs::read_to_string(&lowess_path).unwrap();
        let scalar = fs::read_to_string(&scalar_path).unwrap();
        let _ = fs::remove_file(&num_rej_path);
        let _ = fs::remove_file(&lowess_path);
        let _ = fs::remove_file(&scalar_path);
        assert_eq!(num_rej, "theta\tnumRej\n0\t1\n0.5\t3\n");
        assert_eq!(lowess, "x\ty\n0\t1.25\n0.5\t2.75\n");
        assert_eq!(
            scalar,
            "name\tvalue\nfilterThreshold\t10\nfilterTheta\t0.5\nalpha\t0.1\n"
        );
    }

    #[test]
    fn write_independent_filter_metadata_tsv_rejects_invalid_shapes() {
        let path = unique_test_path("bad_filter_metadata.tsv");
        let mut filtering = IndependentFilteringOutput {
            enabled: true,
            theta: vec![0.0, 0.5],
            num_rejections: vec![1],
            selected_index: Some(1),
            filter_theta: Some(0.5),
            filter_threshold: Some(10.0),
            lowess_fit: Some(vec![1.25, 2.75]),
            alpha: 0.1,
        };

        assert!(write_independent_filter_num_rej_tsv(&path, &filtering).is_err());
        filtering.num_rejections = vec![1, 3];
        filtering.lowess_fit = Some(vec![1.25]);
        assert!(write_independent_filter_lowess_tsv(&path, &filtering).is_err());
        filtering.lowess_fit = Some(vec![1.25, f64::NAN]);
        assert!(write_independent_filter_lowess_tsv(&path, &filtering).is_err());
        filtering.lowess_fit = Some(vec![1.25, 2.75]);
        filtering.theta[1] = 1.2;
        assert!(write_independent_filter_metadata_tsv(&path, &filtering).is_err());
        filtering.theta[1] = 0.5;
        filtering.alpha = f64::NAN;
        assert!(write_independent_filter_metadata_tsv(&path, &filtering).is_err());
        let _ = fs::remove_file(&path);
    }

    fn unique_test_path(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "rsdeseq2_{}_{}_{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("test"),
            name
        ))
    }
}
