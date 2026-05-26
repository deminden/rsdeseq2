use std::path::Path;

use csv::{ReaderBuilder, WriterBuilder};

use crate::core::CountMatrix;
use crate::design::DesignMatrix;
use crate::diagnostics::{
    Deseq2McolsDiagnosticColumn, Deseq2McolsDiagnosticValues, Deseq2McolsDiagnostics,
};
use crate::errors::{invalid_dimensions, DeseqError};
use crate::independent_filtering::IndependentFilteringOutput;
use crate::matrix::RowMajorMatrix;
use crate::results::{DeseqResultColumnValues, DeseqResults};

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
    DesignMatrix::from_row_major(n_samples, n_coefficients, values, Some(coefficient_names))
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
    let path_ref = path.as_ref();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path_ref)?;
    let headers = reader.headers()?.clone();
    if headers.len() < 2 {
        return Err(DeseqError::InvalidDimensions {
            context: "normalization factor columns".to_string(),
            expected: 2,
            actual: headers.len(),
        });
    }
    let n_samples = headers.len() - 1;
    let mut values = Vec::new();
    let mut n_genes = 0_usize;
    for record in reader.records() {
        let record = record?;
        if record.len() != headers.len() {
            return Err(DeseqError::InvalidDimensions {
                context: "normalization factor record columns".to_string(),
                expected: headers.len(),
                actual: record.len(),
            });
        }
        for (sample, field) in record.iter().skip(1).enumerate() {
            values.push(parse_finite_float_field(
                field,
                path_ref,
                n_genes * n_samples + sample,
                "normalization factor matrix",
            )?);
        }
        n_genes += 1;
    }
    if n_genes == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "normalization factor rows".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    let factors = RowMajorMatrix::from_row_major(n_genes, n_samples, values)?;
    validate_positive_finite_matrix("normalization factor matrix", &factors)?;
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
    let path_ref = path.as_ref();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path_ref)?;
    let headers = reader.headers()?.clone();
    if headers.len() < 2 {
        return Err(DeseqError::InvalidDimensions {
            context: "observation weight columns".to_string(),
            expected: 2,
            actual: headers.len(),
        });
    }
    let n_samples = headers.len() - 1;
    let mut values = Vec::new();
    let mut n_genes = 0_usize;
    for record in reader.records() {
        let record = record?;
        if record.len() != headers.len() {
            return Err(DeseqError::InvalidDimensions {
                context: "observation weight record columns".to_string(),
                expected: headers.len(),
                actual: record.len(),
            });
        }
        for (sample, field) in record.iter().skip(1).enumerate() {
            values.push(parse_finite_float_field(
                field,
                path_ref,
                n_genes * n_samples + sample,
                "observation weights",
            )?);
        }
        n_genes += 1;
    }
    if n_genes == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "observation weight rows".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    let weights = RowMajorMatrix::from_row_major(n_genes, n_samples, values)?;
    validate_nonnegative_finite_values("observation weight input", weights.as_slice())?;
    Ok(weights)
}

/// Read a tab-delimited sample-level size-factor table.
///
/// The file shape matches `write_size_factors_tsv`: a leading sample label
/// column and one numeric `size_factor` column. Sample labels are accepted for
/// alignment by the caller; this primitive reader returns only the numeric
/// vector in file order.
pub fn read_size_factors_tsv(path: impl AsRef<Path>) -> Result<Vec<f64>, DeseqError> {
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
        values.push(parse_finite_float_field(
            record.get(1).unwrap_or_default(),
            path_ref,
            idx,
            "size factors",
        )?);
    }
    if values.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "size-factor rows".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    validate_positive_finite_values("size-factor input", &values)?;
    Ok(values)
}

/// Read one supplied geometric mean per gene.
///
/// The file has a leading gene label column and one numeric value column.
/// Labels are accepted for caller-side alignment; this primitive reader returns
/// the values in file order. Values must be non-negative and finite, matching
/// the size-factor estimator's supported geometric-mean domain.
pub fn read_geometric_means_tsv(path: impl AsRef<Path>) -> Result<Vec<f64>, DeseqError> {
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
        values.push(parse_finite_float_field(
            record.get(1).unwrap_or_default(),
            path_ref,
            idx,
            "geometric means",
        )?);
    }
    if values.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "geometric-mean rows".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    for (idx, value) in values.iter().copied().enumerate() {
        if value < 0.0 {
            return Err(DeseqError::InvalidSizeFactors {
                reason: format!("geometric mean at index {idx} must be non-negative"),
            });
        }
    }
    Ok(values)
}

/// Read one finite numeric Wald t degrees-of-freedom value per gene.
///
/// The file has a leading gene label column and one numeric value column.
/// Labels are accepted for caller-side alignment; this primitive reader returns
/// the values in file order.
pub fn read_wald_t_degrees_of_freedom_tsv(path: impl AsRef<Path>) -> Result<Vec<f64>, DeseqError> {
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
        values.push(parse_finite_float_field(
            record.get(1).unwrap_or_default(),
            path_ref,
            idx,
            "Wald t degrees of freedom",
        )?);
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

fn format_r_logical(value: bool) -> String {
    if value {
        "TRUE".to_string()
    } else {
        "FALSE".to_string()
    }
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
