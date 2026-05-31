/// Build DESeq2-shaped Wald result rows from precomputed Wald statistics.
pub fn build_wald_results_from_wald(
    base_mean: &[f64],
    fit: &NbinomGlmFit,
    coefficient: usize,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
    wald: &WaldOutput,
) -> Result<DeseqResults, DeseqError> {
    let n_genes = fit.beta.n_rows();
    validate_result_inputs(base_mean, fit, gene_names, dispersions)?;
    validate_wald_output(wald, n_genes)?;
    if coefficient >= fit.beta.n_cols() {
        return Err(DeseqError::InvalidDimensions {
            context: "Wald result coefficient index".to_string(),
            expected: fit.beta.n_cols().saturating_sub(1),
            actual: coefficient,
        });
    }
    let padj = bh_adjust(&wald.pvalue);

    let mut rows = Vec::with_capacity(n_genes);
    for gene in 0..n_genes {
        let beta = fit.beta.row(gene)?[coefficient];
        let beta_se = fit.beta_se.row(gene)?[coefficient];
        rows.push(DeseqResultRow {
            gene: gene_names.and_then(|names| names.get(gene)).cloned(),
            base_mean: base_mean[gene],
            log2_fold_change: finite_option(beta),
            lfc_se: finite_option(beta_se),
            stat: wald.stat[gene],
            pvalue: wald.pvalue[gene],
            padj: padj[gene],
            dispersion: dispersions
                .and_then(|values| values.get(gene).copied())
                .and_then(finite_option),
            converged: fit.beta_converged.get(gene).copied(),
            max_cooks: None,
            cooks_outlier: None,
            filtered: None,
        });
    }
    Ok(DeseqResults {
        rows,
        metadata: wald_table_metadata(fit, coefficient),
        independent_filtering: None,
    })
}

/// Build DESeq2-shaped Wald result rows from a precomputed numeric contrast.
///
/// This is the result-table companion to the primitive numeric contrast helper.
/// It does not parse R-style contrast specifications.
pub fn build_wald_contrast_results(
    base_mean: &[f64],
    fit: &NbinomGlmFit,
    contrast: &WaldContrastOutput,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    let n_genes = fit.beta.n_rows();
    validate_result_inputs(base_mean, fit, gene_names, dispersions)?;
    validate_wald_contrast_output(contrast, n_genes)?;
    let padj = bh_adjust(&contrast.wald.pvalue);

    let mut rows = Vec::with_capacity(n_genes);
    for gene in 0..n_genes {
        rows.push(DeseqResultRow {
            gene: gene_names.and_then(|names| names.get(gene)).cloned(),
            base_mean: base_mean[gene],
            log2_fold_change: contrast.log2_fold_change[gene],
            lfc_se: contrast.lfc_se[gene],
            stat: contrast.wald.stat[gene],
            pvalue: contrast.wald.pvalue[gene],
            padj: padj[gene],
            dispersion: dispersions
                .and_then(|values| values.get(gene).copied())
                .and_then(finite_option),
            converged: fit.beta_converged.get(gene).copied(),
            max_cooks: None,
            cooks_outlier: None,
            filtered: None,
        });
    }
    Ok(DeseqResults {
        rows,
        metadata: DeseqResultsTableMetadata {
            test_type: Some(TestType::Wald),
            result_name: Some("contrast".to_string()),
            comparison: Some("primitive numeric contrast".to_string()),
            ..DeseqResultsTableMetadata::default()
        },
        independent_filtering: None,
    })
}

/// Build DESeq2-shaped LRT result rows for one reported full-model coefficient.
///
/// DESeq2 reports full-model beta and SE columns alongside the model-level LRT
/// statistic and p-value. This function follows that shape for primitive Rust
/// matrices.
pub fn build_lrt_results(
    base_mean: &[f64],
    full_fit: &NbinomGlmFit,
    lrt: &LrtOutput,
    coefficient: usize,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    validate_result_inputs(base_mean, full_fit, gene_names, dispersions)?;
    if coefficient >= full_fit.beta.n_cols() {
        return Err(DeseqError::InvalidDimensions {
            context: "LRT result coefficient index".to_string(),
            expected: full_fit.beta.n_cols().saturating_sub(1),
            actual: coefficient,
        });
    }
    validate_lrt_output(lrt, full_fit.beta.n_rows())?;
    let padj = bh_adjust(&lrt.pvalue);
    let mut rows = Vec::with_capacity(full_fit.beta.n_rows());
    for gene in 0..full_fit.beta.n_rows() {
        let beta = full_fit.beta.row(gene)?[coefficient];
        let beta_se = full_fit.beta_se.row(gene)?[coefficient];
        rows.push(DeseqResultRow {
            gene: gene_names.and_then(|names| names.get(gene)).cloned(),
            base_mean: base_mean[gene],
            log2_fold_change: finite_option(beta),
            lfc_se: finite_option(beta_se),
            stat: lrt.deviance[gene],
            pvalue: lrt.pvalue[gene],
            padj: padj[gene],
            dispersion: dispersions
                .and_then(|values| values.get(gene).copied())
                .and_then(finite_option),
            converged: full_fit.beta_converged.get(gene).copied(),
            max_cooks: None,
            cooks_outlier: None,
            filtered: None,
        });
    }
    Ok(DeseqResults {
        rows,
        metadata: lrt_table_metadata(full_fit, coefficient),
        independent_filtering: None,
    })
}

/// Build DESeq2-shaped LRT result rows with contrast effect-size columns.
///
/// The LRT statistic and p-value still come from the full-vs-reduced model
/// comparison, while `log2FoldChange` and `lfcSE` report the caller-supplied
/// full-model coefficient contrast. This mirrors the `results()` shape where
/// an LRT result can display a requested effect size without changing the
/// likelihood-ratio test itself.
pub fn build_lrt_contrast_results(
    base_mean: &[f64],
    full_fit: &NbinomGlmFit,
    lrt: &LrtOutput,
    contrast: &WaldContrastOutput,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    let n_genes = full_fit.beta.n_rows();
    validate_result_inputs(base_mean, full_fit, gene_names, dispersions)?;
    validate_lrt_output(lrt, n_genes)?;
    validate_wald_contrast_output(contrast, n_genes)?;
    let padj = bh_adjust(&lrt.pvalue);

    let mut rows = Vec::with_capacity(n_genes);
    for gene in 0..n_genes {
        rows.push(DeseqResultRow {
            gene: gene_names.and_then(|names| names.get(gene)).cloned(),
            base_mean: base_mean[gene],
            log2_fold_change: contrast.log2_fold_change[gene],
            lfc_se: contrast.lfc_se[gene],
            stat: lrt.deviance[gene],
            pvalue: lrt.pvalue[gene],
            padj: padj[gene],
            dispersion: dispersions
                .and_then(|values| values.get(gene).copied())
                .and_then(finite_option),
            converged: full_fit.beta_converged.get(gene).copied(),
            max_cooks: None,
            cooks_outlier: None,
            filtered: None,
        });
    }
    Ok(DeseqResults {
        rows,
        metadata: DeseqResultsTableMetadata {
            test_type: Some(TestType::Lrt),
            result_name: Some("contrast".to_string()),
            comparison: Some("primitive numeric contrast".to_string()),
            ..DeseqResultsTableMetadata::default()
        },
        independent_filtering: None,
    })
}

/// Resolve a Cook's cutoff option for a model with `m` samples and `p` columns.
pub fn resolve_cooks_cutoff(
    cutoff: CooksCutoff,
    n_samples: usize,
    n_coefficients: usize,
) -> Result<Option<f64>, DeseqError> {
    match cutoff {
        CooksCutoff::Disabled => Ok(None),
        CooksCutoff::Threshold(value) => {
            if !value.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: "Cook's cutoff".to_string(),
                    index: None,
                    value,
                });
            }
            Ok(Some(value))
        }
        CooksCutoff::Default => default_cooks_cutoff(n_samples, n_coefficients),
    }
}

/// DESeq2 default Cook's cutoff, `qf(.99, p, m - p)`.
pub fn default_cooks_cutoff(
    n_samples: usize,
    n_coefficients: usize,
) -> Result<Option<f64>, DeseqError> {
    if n_samples <= n_coefficients {
        return Ok(None);
    }
    let df1 = n_coefficients as f64;
    let df2 = (n_samples - n_coefficients) as f64;
    if !df1.is_finite() || df1 <= 0.0 || !df2.is_finite() || df2 <= 0.0 {
        return Err(DeseqError::InvalidDimensions {
            context: "Cook's cutoff degrees of freedom".to_string(),
            expected: n_samples,
            actual: n_coefficients,
        });
    }
    let distribution =
        FisherSnedecor::new(df1, df2).map_err(|error| DeseqError::InvalidDimensions {
            context: format!("Cook's cutoff F distribution: {error}"),
            expected: n_samples,
            actual: n_coefficients,
        })?;
    let cutoff = distribution.inverse_cdf(0.99);
    if !cutoff.is_finite() {
        return Err(DeseqError::NonFiniteValue {
            context: "Cook's default cutoff".to_string(),
            index: None,
            value: cutoff,
        });
    }
    Ok(Some(cutoff))
}

/// Apply Cook's p-value filtering and recompute BH-adjusted p-values.
///
/// This mirrors the default `results()` behavior where rows with `maxCooks`
/// above the selected cutoff have their p-value set to missing before
/// p-value adjustment. Count replacement is handled separately.
pub fn apply_cooks_cutoff(
    results: &mut DeseqResults,
    cutoff: Option<f64>,
) -> Result<(), DeseqError> {
    let Some(cutoff) = cutoff else {
        return Ok(());
    };
    if !cutoff.is_finite() {
        return Err(DeseqError::NonFiniteValue {
            context: "Cook's cutoff".to_string(),
            index: None,
            value: cutoff,
        });
    }

    for (idx, row) in results.rows.iter_mut().enumerate() {
        row.cooks_outlier = match row.max_cooks {
            Some(value) => {
                if !value.is_finite() {
                    return Err(DeseqError::NonFiniteValue {
                        context: "result maxCooks".to_string(),
                        index: Some(idx),
                        value,
                    });
                }
                Some(value > cutoff)
            }
            None => None,
        };
        if row.cooks_outlier == Some(true) {
            row.pvalue = None;
        }
    }
    recompute_padj(results)?;
    Ok(())
}

/// Apply Cook's p-value filtering with DESeq2's two-group low-count heuristic.
///
/// DESeq2's `results()` applies this heuristic only for formula designs with a
/// single two-level factor. The Rust core cannot infer that from primitive
/// matrices, so callers should use this helper only after establishing that
/// condition. For rows above the Cook's cutoff, the helper finds the sample
/// with the maximum Cook's distance. If at least three counts in that row are
/// larger than the count in that sample, the row is not masked.
pub fn apply_cooks_cutoff_with_low_count_heuristic(
    results: &mut DeseqResults,
    cutoff: Option<f64>,
    counts: &CountMatrix,
    cooks: &RowMajorMatrix<f64>,
) -> Result<(), DeseqError> {
    let Some(cutoff) = cutoff else {
        return Ok(());
    };
    if !cutoff.is_finite() {
        return Err(DeseqError::NonFiniteValue {
            context: "Cook's cutoff".to_string(),
            index: None,
            value: cutoff,
        });
    }
    validate_cooks_heuristic_inputs(results, counts, cooks)?;

    for (gene, row) in results.rows.iter_mut().enumerate() {
        let is_outlier = match row.max_cooks {
            Some(value) => {
                if !value.is_finite() {
                    return Err(DeseqError::NonFiniteValue {
                        context: "result maxCooks".to_string(),
                        index: Some(gene),
                        value,
                    });
                }
                Some(value > cutoff)
            }
            None => None,
        };
        row.cooks_outlier = match is_outlier {
            Some(true) => {
                let spare =
                    low_count_outlier_heuristic_spares_row(counts.row(gene)?, cooks.row(gene)?);
                if spare {
                    Some(false)
                } else {
                    row.pvalue = None;
                    Some(true)
                }
            }
            other => other,
        };
    }
    recompute_padj(results)?;
    Ok(())
}

/// Recompute BH-adjusted p-values from the current result p-values.
pub fn recompute_padj(results: &mut DeseqResults) -> Result<(), DeseqError> {
    let pvalues = results
        .rows
        .iter()
        .map(|row| row.pvalue)
        .collect::<Vec<_>>();
    validate_optional_probability(&pvalues, "result p-value")?;
    let padj = bh_adjust(&pvalues);
    for (row, adjusted) in results.rows.iter_mut().zip(padj) {
        row.padj = adjusted;
    }
    Ok(())
}
