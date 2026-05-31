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
