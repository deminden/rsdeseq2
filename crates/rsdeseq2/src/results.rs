use crate::core::CountMatrix;
use crate::errors::{invalid_dimensions, DeseqError};
use crate::glm::{wald_test_coefficient, LrtOutput, NbinomGlmFit, WaldContrastOutput, WaldOutput};
use crate::independent_filtering::IndependentFilteringOutput;
use crate::matrix::RowMajorMatrix;
use crate::multiple_testing::bh_adjust;
use crate::options::CooksCutoff;
use statrs::distribution::{ContinuousCDF, FisherSnedecor};

/// One row of a future DESeq2-like results table.
#[derive(Clone, Debug, PartialEq)]
pub struct DeseqResultRow {
    /// Gene identifier, if available.
    pub gene: Option<String>,
    /// Mean normalized count across samples.
    pub base_mean: f64,
    /// Log2 fold change.
    pub log2_fold_change: Option<f64>,
    /// Log2 fold-change standard error.
    pub lfc_se: Option<f64>,
    /// Test statistic.
    pub stat: Option<f64>,
    /// Raw p-value.
    pub pvalue: Option<f64>,
    /// Adjusted p-value.
    pub padj: Option<f64>,
    /// Final dispersion.
    pub dispersion: Option<f64>,
    /// Convergence flag.
    pub converged: Option<bool>,
    /// Maximum Cook's distance over eligible samples.
    pub max_cooks: Option<f64>,
    /// Whether Cook's cutoff filtering masked this row's p-value.
    pub cooks_outlier: Option<bool>,
    /// Whether independent filtering removed this row from adjusted p-values.
    pub filtered: Option<bool>,
}

/// Collection of result rows.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DeseqResults {
    /// Rows in output order.
    pub rows: Vec<DeseqResultRow>,
    /// Independent-filtering metadata, when result assembly has run that stage.
    pub independent_filtering: Option<IndependentFilteringOutput>,
}

impl DeseqResults {
    /// Number of result rows.
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    /// Whether the result table has no rows.
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }
}

/// Build DESeq2-shaped Wald result rows for one coefficient.
///
/// This mirrors the non-contrast, no-independent-filtering result assembly:
/// `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, and `padj`.
pub fn build_wald_results(
    base_mean: &[f64],
    fit: &NbinomGlmFit,
    coefficient: usize,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    let wald = wald_test_coefficient(fit, coefficient)?;
    build_wald_results_from_wald(base_mean, fit, coefficient, gene_names, dispersions, &wald)
}

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
    if lrt.deviance.len() != full_fit.beta.n_rows() {
        return Err(invalid_dimensions(
            "LRT statistic rows",
            full_fit.beta.n_rows(),
            lrt.deviance.len(),
        ));
    }
    if lrt.pvalue.len() != full_fit.beta.n_rows() {
        return Err(invalid_dimensions(
            "LRT p-value rows",
            full_fit.beta.n_rows(),
            lrt.pvalue.len(),
        ));
    }
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
    let distribution =
        FisherSnedecor::new(df1, df2).map_err(|error| DeseqError::InvalidDimensions {
            context: format!("Cook's cutoff F distribution: {error}"),
            expected: n_samples,
            actual: n_coefficients,
        })?;
    Ok(Some(distribution.inverse_cdf(0.99)))
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

    for row in &mut results.rows {
        row.cooks_outlier = row.max_cooks.map(|value| value > cutoff);
        if row.cooks_outlier == Some(true) {
            row.pvalue = None;
        }
    }
    recompute_padj(results);
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
        let is_outlier = row.max_cooks.map(|value| value > cutoff);
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
    recompute_padj(results);
    Ok(())
}

/// Recompute BH-adjusted p-values from the current result p-values.
pub fn recompute_padj(results: &mut DeseqResults) {
    let pvalues = results
        .rows
        .iter()
        .map(|row| row.pvalue)
        .collect::<Vec<_>>();
    let padj = bh_adjust(&pvalues);
    for (row, adjusted) in results.rows.iter_mut().zip(padj) {
        row.padj = adjusted;
    }
}

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
    if let Some(names) = gene_names {
        if names.len() != n_genes {
            return Err(invalid_dimensions(
                "result gene names",
                n_genes,
                names.len(),
            ));
        }
    }
    if let Some(values) = dispersions {
        if values.len() != n_genes {
            return Err(invalid_dimensions(
                "result dispersions",
                n_genes,
                values.len(),
            ));
        }
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
    }
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
    validate_wald_output(&contrast.wald, n_genes)
}

fn finite_option(value: f64) -> Option<f64> {
    value.is_finite().then_some(value)
}
