use std::path::Path;

use csv::{ReaderBuilder, WriterBuilder};

use crate::core::CountMatrix;
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
    for (idx, value) in matrix.as_slice().iter().copied().enumerate() {
        if !value.is_finite() || value <= 0.0 {
            return Err(DeseqError::InvalidSizeFactors {
                reason: format!("{context} value at index {idx} must be finite and positive"),
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
    let n_rows = diagnostic_frame_row_count(&frame.columns);
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
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["name", "value"])?;
    for entry in filtering.scalar_metadata() {
        writer.write_record([entry.name, entry.value.to_string()])?;
    }
    writer.flush()?;
    Ok(())
}

fn diagnostic_frame_row_count(columns: &[Deseq2McolsDiagnosticColumn]) -> usize {
    columns
        .first()
        .map(|column| match &column.values {
            Deseq2McolsDiagnosticValues::Numeric(values) => values.len(),
            Deseq2McolsDiagnosticValues::OptionalNumeric(values) => values.len(),
            Deseq2McolsDiagnosticValues::Integer(values) => values.len(),
            Deseq2McolsDiagnosticValues::Logical(values) => values.len(),
        })
        .unwrap_or(0)
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
                    base_mean: f64::NAN,
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
                "gene2\tNA\tNA\t0.5\tNA\t0.8\tNA\tNA\tFALSE\tNA\tNA\n"
            )
        );
        assert_eq!(
            tidy_text,
            concat!(
                "row\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\t",
                "dispersion\tconverged\tmaxCooks\tcooksOutlier\n",
                "gene_a\t10\t1.25\tNA\t-2\tNA\t0.05\t0.1\tTRUE\t4\tFALSE\n",
                "row2\tNA\tNA\t0.5\tNA\t0.8\tNA\tNA\tFALSE\tNA\tNA\n"
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

    fn unique_test_path(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "rsdeseq2_{}_{}_{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("test"),
            name
        ))
    }
}
