fn low_count_outlier_heuristic_spares_row(counts: &[u32], cooks: &[f64]) -> bool {
    let Some((max_cook_sample, _)) = cooks
        .iter()
        .copied()
        .enumerate()
        .filter(|(_, value)| value.is_finite())
        .max_by(|(_, left), (_, right)| left.total_cmp(right))
    else {
        return false;
    };
    let outlier_count = counts[max_cook_sample];
    counts
        .iter()
        .filter(|count| **count > outlier_count)
        .count()
        >= 3
}

fn validate_cooks_heuristic_inputs(
    results: &DeseqResults,
    counts: &CountMatrix,
    cooks: &RowMajorMatrix<f64>,
) -> Result<(), DeseqError> {
    if results.rows.len() != counts.n_genes() {
        return Err(invalid_dimensions(
            "Cook's heuristic result rows",
            counts.n_genes(),
            results.rows.len(),
        ));
    }
    if cooks.n_rows() != counts.n_genes() {
        return Err(invalid_dimensions(
            "Cook's heuristic Cook's rows",
            counts.n_genes(),
            cooks.n_rows(),
        ));
    }
    if cooks.n_cols() != counts.n_samples() {
        return Err(invalid_dimensions(
            "Cook's heuristic Cook's columns",
            counts.n_samples(),
            cooks.n_cols(),
        ));
    }
    Ok(())
}

fn numeric_column<F>(rows: &[DeseqResultRow], selector: F) -> DeseqResultColumnValues
where
    F: Fn(&DeseqResultRow) -> Option<f64>,
{
    DeseqResultColumnValues::Numeric(rows.iter().map(selector).collect())
}

fn logical_column<F>(rows: &[DeseqResultRow], selector: F) -> DeseqResultColumnValues
where
    F: Fn(&DeseqResultRow) -> Option<bool>,
{
    DeseqResultColumnValues::Logical(rows.iter().map(selector).collect())
}

fn validate_result_inputs(
    base_mean: &[f64],
    fit: &NbinomGlmFit,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<(), DeseqError> {
    let n_genes = fit.beta.n_rows();
    if base_mean.len() != n_genes {
        return Err(invalid_dimensions(
            "result baseMean",
            n_genes,
            base_mean.len(),
        ));
    }
    for (idx, value) in base_mean.iter().copied().enumerate() {
        if !value.is_finite() || value < 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: "result baseMean".to_string(),
                index: Some(idx),
                value,
            });
        }
    }
    if fit.beta_se.n_rows() != n_genes || fit.beta_se.n_cols() != fit.beta.n_cols() {
        return Err(invalid_dimensions(
            "result betaSE matrix values",
            fit.beta.len(),
            fit.beta_se.len(),
        ));
    }
    if fit.beta_converged.len() != n_genes {
        return Err(invalid_dimensions(
            "result beta convergence flags",
            n_genes,
            fit.beta_converged.len(),
        ));
    }
    if let Some(names) = gene_names
        && names.len() != n_genes {
            return Err(invalid_dimensions(
                "result gene names",
                n_genes,
                names.len(),
            ));
        }
    if let Some(values) = dispersions
        && values.len() != n_genes {
            return Err(invalid_dimensions(
                "result dispersions",
                n_genes,
                values.len(),
            ));
        }
    Ok(())
}

fn validate_wald_output(wald: &WaldOutput, n_genes: usize) -> Result<(), DeseqError> {
    if wald.stat.len() != n_genes {
        return Err(invalid_dimensions(
            "Wald result statistic rows",
            n_genes,
            wald.stat.len(),
        ));
    }
    if wald.pvalue.len() != n_genes {
        return Err(invalid_dimensions(
            "Wald result p-value rows",
            n_genes,
            wald.pvalue.len(),
        ));
    }
    if let Some(df) = &wald.degrees_of_freedom {
        if df.len() != n_genes {
            return Err(invalid_dimensions(
                "Wald result degrees-of-freedom rows",
                n_genes,
                df.len(),
            ));
        }
        validate_optional_positive_finite(df, "Wald result degrees of freedom")?;
    }
    validate_optional_finite(&wald.stat, "Wald result statistic")?;
    validate_optional_probability(&wald.pvalue, "Wald result p-value")?;
    Ok(())
}

fn validate_wald_contrast_output(
    contrast: &WaldContrastOutput,
    n_genes: usize,
) -> Result<(), DeseqError> {
    if contrast.log2_fold_change.len() != n_genes {
        return Err(invalid_dimensions(
            "Wald contrast estimate rows",
            n_genes,
            contrast.log2_fold_change.len(),
        ));
    }
    if contrast.lfc_se.len() != n_genes {
        return Err(invalid_dimensions(
            "Wald contrast SE rows",
            n_genes,
            contrast.lfc_se.len(),
        ));
    }
    validate_optional_finite(&contrast.log2_fold_change, "Wald contrast estimate")?;
    validate_optional_finite(&contrast.lfc_se, "Wald contrast SE")?;
    validate_wald_output(&contrast.wald, n_genes)
}

fn validate_lrt_output(lrt: &LrtOutput, n_genes: usize) -> Result<(), DeseqError> {
    if lrt.deviance.len() != n_genes {
        return Err(invalid_dimensions(
            "LRT statistic rows",
            n_genes,
            lrt.deviance.len(),
        ));
    }
    if lrt.pvalue.len() != n_genes {
        return Err(invalid_dimensions(
            "LRT p-value rows",
            n_genes,
            lrt.pvalue.len(),
        ));
    }
    if lrt.reduced_converged.len() != n_genes {
        return Err(invalid_dimensions(
            "LRT reduced convergence flags",
            n_genes,
            lrt.reduced_converged.len(),
        ));
    }
    if lrt.degrees_of_freedom == 0 {
        return Err(DeseqError::InvalidOptions {
            reason: "LRT degrees of freedom must be positive".to_string(),
        });
    }
    validate_optional_finite(&lrt.deviance, "LRT statistic")?;
    validate_optional_probability(&lrt.pvalue, "LRT p-value")?;
    Ok(())
}

fn validate_optional_finite(values: &[Option<f64>], context: &str) -> Result<(), DeseqError> {
    for (idx, value) in values.iter().copied().enumerate() {
        if let Some(value) = value
            && !value.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: context.to_string(),
                    index: Some(idx),
                    value,
                });
            }
    }
    Ok(())
}

fn validate_optional_positive_finite(
    values: &[Option<f64>],
    context: &str,
) -> Result<(), DeseqError> {
    for (idx, value) in values.iter().copied().enumerate() {
        if let Some(value) = value
            && (!value.is_finite() || value <= 0.0) {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("{context} at index {idx} must be positive and finite"),
                });
            }
    }
    Ok(())
}

fn validate_optional_probability(values: &[Option<f64>], context: &str) -> Result<(), DeseqError> {
    for (idx, value) in values.iter().copied().enumerate() {
        if let Some(value) = value
            && (!value.is_finite() || !(0.0..=1.0).contains(&value)) {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("{context} at index {idx} must be finite and within [0, 1]"),
                });
            }
    }
    Ok(())
}

fn wald_table_metadata(fit: &NbinomGlmFit, coefficient: usize) -> DeseqResultsTableMetadata {
    DeseqResultsTableMetadata {
        test_type: Some(TestType::Wald),
        result_name: Some(result_name_for_coefficient(fit, coefficient)),
        ..DeseqResultsTableMetadata::default()
    }
}

fn lrt_table_metadata(fit: &NbinomGlmFit, coefficient: usize) -> DeseqResultsTableMetadata {
    DeseqResultsTableMetadata {
        test_type: Some(TestType::Lrt),
        result_name: Some(result_name_for_coefficient(fit, coefficient)),
        comparison: Some("full model versus reduced model".to_string()),
        ..DeseqResultsTableMetadata::default()
    }
}

fn result_name_for_coefficient(fit: &NbinomGlmFit, coefficient: usize) -> String {
    fit.model_matrix
        .coefficient_names()
        .and_then(|names| names.get(coefficient))
        .cloned()
        .unwrap_or_else(|| format!("coefficient_{coefficient}"))
}

fn result_column_type(name: &str) -> &'static str {
    match name {
        "dispersion" | "converged" | "maxCooks" | "cooksOutlier" | "filtered" => "diagnostic",
        _ => "results",
    }
}

fn result_column_description(name: &str, metadata: &DeseqResultsTableMetadata) -> String {
    // Column descriptions combine stable DESeq2 column names with the selected
    // contrast/test metadata so sidecar exports stay informative but compact.
    match name {
        "baseMean" => "mean of normalized counts for all samples".to_string(),
        "log2FoldChange" => effect_description(metadata, "log2 fold change (MLE)"),
        "lfcSE" => effect_description(metadata, "standard error"),
        "stat" => statistic_description(metadata),
        "pvalue" => pvalue_description(metadata),
        "padj" => format!("{} adjusted p-values", metadata.p_adjust_method),
        "dispersion" => "final dispersion estimate".to_string(),
        "converged" => "whether beta fitting converged".to_string(),
        "maxCooks" => "maximum Cook's distance over eligible samples".to_string(),
        "cooksOutlier" => "whether Cook's cutoff masked the p-value".to_string(),
        "filtered" => "whether independent filtering removed this row".to_string(),
        _ => "result column".to_string(),
    }
}

fn effect_description(metadata: &DeseqResultsTableMetadata, prefix: &str) -> String {
    match effect_description_label(metadata) {
        Some(label) => format!("{prefix}: {label}"),
        None => prefix.to_string(),
    }
}

fn statistic_description(metadata: &DeseqResultsTableMetadata) -> String {
    match metadata.test_type {
        Some(TestType::Wald) => {
            labelled_description("Wald statistic", test_description_label(metadata))
        }
        Some(TestType::Lrt) => {
            labelled_description("LRT statistic", test_description_label(metadata))
        }
        None => "test statistic".to_string(),
    }
}

fn pvalue_description(metadata: &DeseqResultsTableMetadata) -> String {
    match metadata.test_type {
        Some(TestType::Wald) => {
            labelled_description("Wald test p-value", test_description_label(metadata))
        }
        Some(TestType::Lrt) => {
            labelled_description("LRT p-value", test_description_label(metadata))
        }
        None => "Wald or likelihood-ratio test p-value".to_string(),
    }
}

fn labelled_description(prefix: &str, label: Option<&str>) -> String {
    match label {
        Some(label) => format!("{prefix}: {label}"),
        None => prefix.to_string(),
    }
}

fn effect_description_label(metadata: &DeseqResultsTableMetadata) -> Option<&str> {
    match metadata.test_type {
        Some(TestType::Lrt) => metadata
            .result_name
            .as_deref()
            .or(metadata.comparison.as_deref()),
        _ => metadata
            .comparison
            .as_deref()
            .or(metadata.result_name.as_deref()),
    }
}

fn test_description_label(metadata: &DeseqResultsTableMetadata) -> Option<&str> {
    metadata
        .comparison
        .as_deref()
        .or(metadata.result_name.as_deref())
}

fn test_type_label(test_type: TestType) -> &'static str {
    match test_type {
        TestType::Wald => "Wald",
        TestType::Lrt => "LRT",
    }
}

fn wald_alternative_name(alternative: WaldAlternative) -> &'static str {
    match alternative {
        WaldAlternative::GreaterAbs => "greaterAbs",
        WaldAlternative::GreaterAbsUpshot => "greaterAbsUPSHOT",
        WaldAlternative::GreaterAbs2014 => "greaterAbs2014",
        WaldAlternative::LessAbs => "lessAbs",
        WaldAlternative::Greater => "greater",
        WaldAlternative::Less => "less",
    }
}

fn finite_option(value: f64) -> Option<f64> {
    value.is_finite().then_some(value)
}
