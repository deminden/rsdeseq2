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

/// Read one finite numeric value per sample with labels preserved.
pub fn read_labeled_sample_numeric_tsv(
    path: impl AsRef<Path>,
    context: &str,
) -> Result<Vec<SampleNumericValue>, DeseqError> {
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
        let sample = record.get(0).unwrap_or_default();
        if sample.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("{context} sample names must not be empty"),
            });
        }
        let value =
            parse_finite_float_field(record.get(1).unwrap_or_default(), path_ref, idx, context)?;
        values.push(SampleNumericValue {
            sample: sample.to_string(),
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
