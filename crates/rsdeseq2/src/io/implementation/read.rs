use super::*;
use crate::errors::{ErrorKind, Requirement, invalid_parameter};
pub(super) use std::collections::{HashMap, HashSet};
pub(super) use std::io::Read;
pub(super) use std::path::Path;

pub(super) use csv::{ReaderBuilder, WriterBuilder};

pub(super) use crate::cooks::{CooksOutput, CooksRefitPlan};
pub(super) use crate::design::DesignMatrix;
pub(super) use crate::diagnostics::{
    Deseq2McolsDiagnosticColumn, Deseq2McolsDiagnosticValues, Deseq2McolsDiagnostics,
};
pub(super) use crate::errors::{DeseqError, invalid_dimensions};
pub(super) use crate::independent_filtering::IndependentFilteringOutput;
pub(super) use crate::matrix::CountMatrix;
pub(super) use crate::matrix::RowMajorMatrix;
pub(super) use crate::results::{DeseqResultColumnValues, DeseqResults};

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
    read_count_matrix_tsv_records(tsv_from_path(path_ref)?, path_ref)
}

/// Like [`read_count_matrix_tsv`], reading from any [`Read`] source such as an in-memory
/// buffer; error messages name the source `<reader>`.
pub fn read_count_matrix_tsv_from_reader(reader: impl Read) -> Result<CountMatrix, DeseqError> {
    read_count_matrix_tsv_records(tsv_from_reader(reader), Path::new("<reader>"))
}

fn read_count_matrix_tsv_records<R: Read>(
    reader: csv::Reader<R>,
    path_ref: &Path,
) -> Result<CountMatrix, DeseqError> {
    let TsvMatrix {
        rows: gene_names,
        columns: sample_names,
        values,
    } = read_labelled_matrix(reader, MatrixTsvSchema::Counts, |_, field| {
        parse_count_field(field, path_ref)
    })?;
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
    read_labeled_design_matrix_tsv_records(tsv_from_path(path_ref)?, path_ref)
}

/// Like [`read_labeled_design_matrix_tsv`], reading from any [`Read`] source such as an in-memory
/// buffer; error messages name the source `<reader>`.
pub fn read_labeled_design_matrix_tsv_from_reader(
    reader: impl Read,
) -> Result<LabeledDesignMatrix, DeseqError> {
    read_labeled_design_matrix_tsv_records(tsv_from_reader(reader), Path::new("<reader>"))
}

fn read_labeled_design_matrix_tsv_records<R: Read>(
    reader: csv::Reader<R>,
    path_ref: &Path,
) -> Result<LabeledDesignMatrix, DeseqError> {
    let TsvMatrix {
        rows: sample_names,
        columns: coefficient_names,
        values,
    } = read_labelled_matrix(
        reader,
        MatrixTsvSchema::Numeric {
            columns: "design matrix columns",
            record: "design matrix record columns",
            rows: "design matrix rows",
            label: "design sample names",
        },
        |index, field| parse_finite_float_field(field, path_ref, index, "design matrix"),
    )?;
    let design = DesignMatrix::from_row_major(
        sample_names.len(),
        coefficient_names.len(),
        values,
        Some(coefficient_names),
    )?;
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
    let TsvMatrix {
        rows: gene_names,
        columns: sample_names,
        values,
    } = read_labelled_matrix(
        tsv_from_path(path_ref)?,
        MatrixTsvSchema::Numeric {
            columns: column_context,
            record: record_context,
            rows: row_context,
            label: gene_context,
        },
        |index, field| parse_finite_float_field(field, path_ref, index, value_context),
    )?;
    let matrix = RowMajorMatrix::from_row_major(gene_names.len(), sample_names.len(), values)?;
    Ok(LabeledAssayMatrix {
        matrix,
        gene_names,
        sample_names,
    })
}

enum MatrixTsvSchema<'a> {
    Counts,
    Numeric {
        columns: &'a str,
        record: &'a str,
        rows: &'a str,
        label: &'a str,
    },
}

struct TsvMatrix<T> {
    rows: Vec<String>,
    columns: Vec<String>,
    values: Vec<T>,
}

fn read_labelled_matrix<R: Read, T>(
    mut reader: csv::Reader<R>,
    schema: MatrixTsvSchema<'_>,
    mut cell: impl FnMut(usize, &str) -> Result<T, DeseqError>,
) -> Result<TsvMatrix<T>, DeseqError> {
    let headers = reader.headers()?;
    let width = headers.len();
    if width < 2 {
        return Err(match schema {
            MatrixTsvSchema::Counts => DeseqError::InvalidCounts {
                reason: "count table must have a gene column and at least one sample column".into(),
            },
            MatrixTsvSchema::Numeric { columns, .. } => invalid_dimensions(columns, 2, width),
        });
    }
    let columns = headers.iter().skip(1).map(ToOwned::to_owned).collect();
    let mut rows = Vec::new();
    let mut values = Vec::new();
    let mut record = csv::StringRecord::new();
    while reader.read_record(&mut record)? {
        if record.len() != width {
            return Err(match schema {
                MatrixTsvSchema::Counts => DeseqError::InvalidCounts {
                    reason: format!(
                        "record has {} columns but header has {width} columns",
                        record.len()
                    ),
                },
                MatrixTsvSchema::Numeric {
                    record: context, ..
                } => invalid_dimensions(context, width, record.len()),
            });
        }
        let label = record.get(0).unwrap_or_default();
        if label.is_empty()
            && let MatrixTsvSchema::Numeric { label: context, .. } = schema
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!("{context} must not be empty"),
            });
        }
        rows.push(label.to_owned());
        for field in record.iter().skip(1) {
            values.push(cell(values.len(), field)?);
        }
    }
    if rows.is_empty()
        && let MatrixTsvSchema::Numeric { rows: context, .. } = schema
    {
        return Err(invalid_dimensions(context, 1, 0));
    }
    Ok(TsvMatrix {
        rows,
        columns,
        values,
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
    read_labeled_size_factors_tsv_records(tsv_from_path(path_ref)?, path_ref)
}

/// Like [`read_labeled_size_factors_tsv`], reading from any [`Read`] source such as an in-memory
/// buffer; error messages name the source `<reader>`.
pub fn read_labeled_size_factors_tsv_from_reader(
    reader: impl Read,
) -> Result<Vec<SampleNumericValue>, DeseqError> {
    read_labeled_size_factors_tsv_records(tsv_from_reader(reader), Path::new("<reader>"))
}

fn read_labeled_size_factors_tsv_records<R: Read>(
    reader: csv::Reader<R>,
    path_ref: &Path,
) -> Result<Vec<SampleNumericValue>, DeseqError> {
    let values = read_labelled(
        reader,
        "size-factor",
        "size-factor sample names must not be empty",
        |idx, sample, value| {
            Ok(SampleNumericValue {
                sample: sample.to_string(),
                value: parse_finite_float_field(value, path_ref, idx, "size factors")?,
            })
        },
    )?;
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
    read_labeled_geometric_means_tsv_records(tsv_from_path(path_ref)?, path_ref)
}

/// Like [`read_labeled_geometric_means_tsv`], reading from any [`Read`] source such as an in-memory
/// buffer; error messages name the source `<reader>`.
pub fn read_labeled_geometric_means_tsv_from_reader(
    reader: impl Read,
) -> Result<Vec<GeneNumericValue>, DeseqError> {
    read_labeled_geometric_means_tsv_records(tsv_from_reader(reader), Path::new("<reader>"))
}

fn read_labeled_geometric_means_tsv_records<R: Read>(
    reader: csv::Reader<R>,
    path_ref: &Path,
) -> Result<Vec<GeneNumericValue>, DeseqError> {
    let values = read_labelled(
        reader,
        "geometric-mean",
        "geometric-mean gene names must not be empty",
        |idx, gene, value| {
            Ok(GeneNumericValue {
                gene: gene.to_string(),
                value: parse_finite_float_field(value, path_ref, idx, "geometric means")?,
            })
        },
    )?;
    for (idx, value) in values.iter().map(|entry| entry.value).enumerate() {
        if value < 0.0 {
            return Err(invalid_parameter(
                ErrorKind::Normalization,
                format!("geometric mean at index {idx}"),
                Requirement::NonNegative,
            ));
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
    read_labeled_sample_numeric_tsv_records(tsv_from_path(path_ref)?, path_ref, context)
}

/// Like [`read_labeled_sample_numeric_tsv`], reading from any [`Read`] source such as an in-memory
/// buffer; error messages name the source `<reader>`.
pub fn read_labeled_sample_numeric_tsv_from_reader(
    reader: impl Read,
    context: &str,
) -> Result<Vec<SampleNumericValue>, DeseqError> {
    read_labeled_sample_numeric_tsv_records(tsv_from_reader(reader), Path::new("<reader>"), context)
}

fn read_labeled_sample_numeric_tsv_records<R: Read>(
    reader: csv::Reader<R>,
    path_ref: &Path,
    context: &str,
) -> Result<Vec<SampleNumericValue>, DeseqError> {
    read_labelled(
        reader,
        context,
        &format!("{context} sample names must not be empty"),
        |idx, sample, value| {
            Ok(SampleNumericValue {
                sample: sample.to_string(),
                value: parse_finite_float_field(value, path_ref, idx, context)?,
            })
        },
    )
}

/// Read one finite numeric value per gene with labels preserved.
pub fn read_labeled_gene_numeric_tsv(
    path: impl AsRef<Path>,
    context: &str,
) -> Result<Vec<GeneNumericValue>, DeseqError> {
    let path_ref = path.as_ref();
    read_labeled_gene_numeric_tsv_records(tsv_from_path(path_ref)?, path_ref, context)
}

/// Like [`read_labeled_gene_numeric_tsv`], reading from any [`Read`] source such as an in-memory
/// buffer; error messages name the source `<reader>`.
pub fn read_labeled_gene_numeric_tsv_from_reader(
    reader: impl Read,
    context: &str,
) -> Result<Vec<GeneNumericValue>, DeseqError> {
    read_labeled_gene_numeric_tsv_records(tsv_from_reader(reader), Path::new("<reader>"), context)
}

fn read_labeled_gene_numeric_tsv_records<R: Read>(
    reader: csv::Reader<R>,
    path_ref: &Path,
    context: &str,
) -> Result<Vec<GeneNumericValue>, DeseqError> {
    read_labelled(
        reader,
        context,
        &format!("{context} gene names must not be empty"),
        |idx, gene, value| {
            Ok(GeneNumericValue {
                gene: gene.to_string(),
                value: parse_finite_float_field(value, path_ref, idx, context)?,
            })
        },
    )
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
    read_labeled_wald_t_degrees_of_freedom_tsv_records(tsv_from_path(path_ref)?, path_ref)
}

/// Like [`read_labeled_wald_t_degrees_of_freedom_tsv`], reading from any [`Read`] source such as an in-memory
/// buffer; error messages name the source `<reader>`.
pub fn read_labeled_wald_t_degrees_of_freedom_tsv_from_reader(
    reader: impl Read,
) -> Result<Vec<GeneNumericValue>, DeseqError> {
    read_labeled_wald_t_degrees_of_freedom_tsv_records(
        tsv_from_reader(reader),
        Path::new("<reader>"),
    )
}

fn read_labeled_wald_t_degrees_of_freedom_tsv_records<R: Read>(
    reader: csv::Reader<R>,
    path_ref: &Path,
) -> Result<Vec<GeneNumericValue>, DeseqError> {
    read_labelled(
        reader,
        "Wald t degrees-of-freedom",
        "Wald t degrees-of-freedom gene names must not be empty",
        |idx, gene, value| {
            Ok(GeneNumericValue {
                gene: gene.to_string(),
                value: parse_finite_float_field(value, path_ref, idx, "Wald t degrees of freedom")?,
            })
        },
    )
}

/// Read one factor level per sample for CLI factor-level contrasts.
///
/// The file has a leading sample label column and one string level column.
/// Labels are preserved so callers can align levels against count-matrix sample
/// columns before applying DESeq2-style factor-level contrast handling.
pub fn read_sample_levels_tsv(path: impl AsRef<Path>) -> Result<Vec<SampleLevel>, DeseqError> {
    let path_ref = path.as_ref();
    read_sample_levels_tsv_records(tsv_from_path(path_ref)?)
}

/// Like [`read_sample_levels_tsv`], reading from any [`Read`] source such as an in-memory
/// buffer; error messages name the source `<reader>`.
pub fn read_sample_levels_tsv_from_reader(
    reader: impl Read,
) -> Result<Vec<SampleLevel>, DeseqError> {
    read_sample_levels_tsv_records(tsv_from_reader(reader))
}

fn read_sample_levels_tsv_records<R: Read>(
    reader: csv::Reader<R>,
) -> Result<Vec<SampleLevel>, DeseqError> {
    read_labelled(
        reader,
        "sample-level",
        "sample names must not be empty",
        |row, sample, level| {
            if level.is_empty() {
                return Err(DeseqError::EmptySampleLevel {
                    row,
                    sample: sample.to_string(),
                });
            }
            Ok(SampleLevel {
                sample: sample.to_string(),
                level: level.to_string(),
            })
        },
    )
}

/// Read a two-column labelled table: a header, then one non-empty label and
/// one value per row, at least one row. `cell` turns the row index, label
/// and value text into an entry. Errors name `{what} columns`,
/// `{what} record columns` and `{what} rows`.
fn read_labelled<R: Read, T>(
    mut reader: csv::Reader<R>,
    what: &str,
    empty_label: &str,
    mut cell: impl FnMut(usize, &str, &str) -> Result<T, DeseqError>,
) -> Result<Vec<T>, DeseqError> {
    let headers = reader.headers()?.clone();
    if headers.len() != 2 {
        return Err(DeseqError::InvalidDimensions {
            context: format!("{what} columns"),
            expected: 2,
            actual: headers.len(),
        });
    }
    let mut values = Vec::new();
    let mut record = csv::StringRecord::new();
    while reader.read_record(&mut record)? {
        let idx = values.len();
        if record.len() != 2 {
            return Err(DeseqError::InvalidDimensions {
                context: format!("{what} record columns"),
                expected: 2,
                actual: record.len(),
            });
        }
        let label = record.get(0).unwrap_or_default();
        if label.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: empty_label.to_string(),
            });
        }
        values.push(cell(idx, label, record.get(1).unwrap_or_default())?);
    }
    if values.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: format!("{what} rows"),
            expected: 1,
            actual: 0,
        });
    }
    Ok(values)
}

fn tsv_from_path(path: &Path) -> Result<csv::Reader<std::fs::File>, csv::Error> {
    ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
}

fn tsv_from_reader<R: Read>(reader: R) -> csv::Reader<R> {
    ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(reader)
}
