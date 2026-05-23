/// Diagnostic counters for future fitting stages.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct DiagnosticSummary {
    /// Number of genes attempted by a fitting stage.
    pub attempted_genes: usize,
    /// Number of genes that converged.
    pub converged_genes: usize,
    /// Number of genes routed through a fallback.
    pub fallback_genes: usize,
}

use crate::core::DeseqFit;

/// DESeq2-style row diagnostics derived from an inspectable [`DeseqFit`].
///
/// These names mirror the metadata columns DESeq2 stores in `mcols(dds)` after
/// Wald and LRT pipelines. The Rust fit state keeps more explicit field names;
/// this view exists to make R and parity wrappers straightforward without
/// duplicating data in the core structs.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Deseq2McolsDiagnostics {
    /// DESeq2 `dispGeneIter` column from gene-wise dispersion fitting.
    pub disp_gene_iter: Option<Vec<usize>>,
    /// Wald-style `betaConv` column.
    pub beta_conv: Option<Vec<bool>>,
    /// LRT-style `fullBetaConv` column.
    pub full_beta_conv: Option<Vec<bool>>,
    /// LRT-style `reducedBetaConv` column.
    pub reduced_beta_conv: Option<Vec<bool>>,
    /// DESeq2 `betaIter` column for the full model.
    pub beta_iter: Option<Vec<usize>>,
    /// Reduced-model beta iterations retained for Rust/R parity diagnostics.
    pub reduced_beta_iter: Option<Vec<usize>>,
    /// DESeq2 `deviance` column, equal to `-2 * full logLike`.
    pub deviance: Option<Vec<f64>>,
    /// DESeq2 `maxCooks` column.
    pub max_cooks: Option<Vec<Option<f64>>>,
}

impl DeseqFit {
    /// Return a DESeq2-metadata-shaped diagnostic view.
    ///
    /// For LRT fits, beta convergence is exposed as `fullBetaConv` and
    /// `reducedBetaConv`, matching DESeq2's `nbinomLRT` metadata. For non-LRT
    /// GLM fits, convergence is exposed as Wald-style `betaConv`.
    pub fn deseq2_mcols_diagnostics(&self) -> Deseq2McolsDiagnostics {
        let is_lrt = self.lrt.is_some();
        Deseq2McolsDiagnostics {
            disp_gene_iter: self.disp_gene_iter.clone(),
            beta_conv: (!is_lrt).then(|| self.beta_converged.clone()).flatten(),
            full_beta_conv: is_lrt.then(|| self.beta_converged.clone()).flatten(),
            reduced_beta_conv: self
                .reduced_beta_converged
                .clone()
                .or_else(|| self.lrt.as_ref().map(|lrt| lrt.reduced_converged.clone())),
            beta_iter: self.beta_iter.clone(),
            reduced_beta_iter: self.reduced_beta_iter.clone(),
            deviance: self.full_deviance.clone(),
            max_cooks: self.max_cooks.clone(),
        }
    }
}
