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

pub fn write_optional_numeric_matrix_tsv(
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
    if let Some(names) = names
        && names.len() != expected {
            return Err(invalid_dimensions(context, expected, names.len()));
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
