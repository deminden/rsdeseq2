/// Build DESeq2-shaped Wald result rows from precomputed Wald statistics.
pub fn build_wald_results_from_wald(
    base_mean: &[f64],
    fit: &NbinomGlmFit,
    coefficient: usize,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
    wald: &WaldOutput,
) -> Result<DeseqResults, DeseqError> {
    assemble_wald_coefficient(
        base_mean,
        fit,
        coefficient,
        wald,
        RowAssembly::adjusted(gene_names, dispersions),
    )
}

pub(crate) fn build_wald_results_from_wald_unadjusted(
    base_mean: &[f64],
    fit: &NbinomGlmFit,
    coefficient: usize,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
    wald: &WaldOutput,
) -> Result<DeseqResults, DeseqError> {
    assemble_wald_coefficient(
        base_mean,
        fit,
        coefficient,
        wald,
        RowAssembly::unadjusted(gene_names, dispersions),
    )
}

fn assemble_wald_coefficient(
    base_mean: &[f64],
    fit: &NbinomGlmFit,
    coefficient: usize,
    wald: &WaldOutput,
    options: RowAssembly<'_>,
) -> Result<DeseqResults, DeseqError> {
    let RowAssembly {
        gene_names,
        dispersions,
        adjust,
    } = options;
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
    let padj = adjust.then(|| bh_adjust(&wald.pvalue));

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
            padj: padj.as_ref().and_then(|values| values[gene]),
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
    build_wald_contrast_results_view(
        base_mean,
        &GlmInferenceView::from(fit),
        contrast,
        gene_names,
        dispersions,
    )
}

pub(crate) fn build_wald_contrast_results_view(
    base_mean: &[f64],
    fit: &GlmInferenceView<'_>,
    contrast: &WaldContrastOutput,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    assemble_wald_contrast(
        base_mean,
        fit,
        contrast,
        RowAssembly::adjusted(gene_names, dispersions),
    )
}

pub(crate) fn build_wald_contrast_results_unadjusted(
    base_mean: &[f64],
    fit: &NbinomGlmFit,
    contrast: &WaldContrastOutput,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    assemble_wald_contrast(
        base_mean,
        &GlmInferenceView::from(fit),
        contrast,
        RowAssembly::unadjusted(gene_names, dispersions),
    )
}

fn assemble_wald_contrast(
    base_mean: &[f64],
    fit: &GlmInferenceView<'_>,
    contrast: &WaldContrastOutput,
    options: RowAssembly<'_>,
) -> Result<DeseqResults, DeseqError> {
    let RowAssembly {
        gene_names,
        dispersions,
        adjust,
    } = options;
    let n_genes = fit.beta.n_rows();
    validate_result_view(base_mean, fit, gene_names, dispersions)?;
    validate_wald_contrast_output(contrast, n_genes)?;
    let padj = adjust.then(|| bh_adjust(&contrast.wald.pvalue));

    let mut rows = Vec::with_capacity(n_genes);
    for gene in 0..n_genes {
        rows.push(DeseqResultRow {
            gene: gene_names.and_then(|names| names.get(gene)).cloned(),
            base_mean: base_mean[gene],
            log2_fold_change: contrast.log2_fold_change[gene],
            lfc_se: contrast.lfc_se[gene],
            stat: contrast.wald.stat[gene],
            pvalue: contrast.wald.pvalue[gene],
            padj: padj.as_ref().and_then(|values| values[gene]),
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
