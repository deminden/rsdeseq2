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
    assemble_lrt_coefficient(
        base_mean,
        full_fit,
        lrt,
        coefficient,
        RowAssembly::adjusted(gene_names, dispersions),
    )
}

pub(crate) fn build_lrt_results_unadjusted(
    base_mean: &[f64],
    full_fit: &NbinomGlmFit,
    lrt: &LrtOutput,
    coefficient: usize,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    assemble_lrt_coefficient(
        base_mean,
        full_fit,
        lrt,
        coefficient,
        RowAssembly::unadjusted(gene_names, dispersions),
    )
}

fn assemble_lrt_coefficient(
    base_mean: &[f64],
    full_fit: &NbinomGlmFit,
    lrt: &LrtOutput,
    coefficient: usize,
    options: RowAssembly<'_>,
) -> Result<DeseqResults, DeseqError> {
    let RowAssembly {
        gene_names,
        dispersions,
        adjust,
    } = options;
    validate_result_inputs(base_mean, full_fit, gene_names, dispersions)?;
    if coefficient >= full_fit.beta.n_cols() {
        return Err(DeseqError::InvalidDimensions {
            context: "LRT result coefficient index".to_string(),
            expected: full_fit.beta.n_cols().saturating_sub(1),
            actual: coefficient,
        });
    }
    validate_lrt_output(lrt, full_fit.beta.n_rows())?;
    let padj = adjust.then(|| bh_adjust(&lrt.pvalue));
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
            padj: padj.as_ref().and_then(|values| values[gene]),
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
    build_lrt_contrast_results_view(
        base_mean,
        &GlmInferenceView::from(full_fit),
        lrt,
        contrast,
        gene_names,
        dispersions,
    )
}

pub(crate) fn build_lrt_contrast_results_view(
    base_mean: &[f64],
    full_fit: &GlmInferenceView<'_>,
    lrt: &LrtOutput,
    contrast: &WaldContrastOutput,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    validate_wald_contrast_output(contrast, full_fit.beta.n_rows())?;
    build_lrt_effect_results_view(
        base_mean,
        full_fit,
        lrt,
        &contrast.log2_fold_change,
        &contrast.lfc_se,
        RowAssembly::adjusted(gene_names, dispersions),
    )
}

pub(crate) fn build_lrt_projection_results_view(
    base_mean: &[f64],
    full_fit: &GlmInferenceView<'_>,
    lrt: &LrtOutput,
    projection: &crate::glm::ContrastProjection,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    build_lrt_effect_results_view(
        base_mean,
        full_fit,
        lrt,
        &projection.log2_fold_change,
        &projection.lfc_se,
        RowAssembly::adjusted(gene_names, dispersions),
    )
}

pub(crate) fn build_lrt_projection_results_unadjusted_view(
    base_mean: &[f64],
    full_fit: &GlmInferenceView<'_>,
    lrt: &LrtOutput,
    projection: &crate::glm::ContrastProjection,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    build_lrt_effect_results_view(
        base_mean,
        full_fit,
        lrt,
        &projection.log2_fold_change,
        &projection.lfc_se,
        RowAssembly::unadjusted(gene_names, dispersions),
    )
}

fn build_lrt_effect_results_view(
    base_mean: &[f64],
    full_fit: &GlmInferenceView<'_>,
    lrt: &LrtOutput,
    effect: &[Option<f64>],
    se: &[Option<f64>],
    options: RowAssembly<'_>,
) -> Result<DeseqResults, DeseqError> {
    let RowAssembly {
        gene_names,
        dispersions,
        adjust,
    } = options;
    let n_genes = full_fit.beta.n_rows();
    validate_result_view(base_mean, full_fit, gene_names, dispersions)?;
    validate_lrt_output(lrt, n_genes)?;
    if effect.len() != n_genes {
        return Err(invalid_dimensions(
            "LRT contrast effects",
            n_genes,
            effect.len(),
        ));
    }
    if se.len() != n_genes {
        return Err(invalid_dimensions("LRT contrast SEs", n_genes, se.len()));
    }
    let padj = adjust.then(|| bh_adjust(&lrt.pvalue));

    let mut rows = Vec::with_capacity(n_genes);
    for gene in 0..n_genes {
        rows.push(DeseqResultRow {
            gene: gene_names.and_then(|names| names.get(gene)).cloned(),
            base_mean: base_mean[gene],
            log2_fold_change: effect[gene],
            lfc_se: se[gene],
            stat: lrt.deviance[gene],
            pvalue: lrt.pvalue[gene],
            padj: padj.as_ref().and_then(|values| values[gene]),
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
