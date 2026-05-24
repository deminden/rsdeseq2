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

/// Typed values for one DESeq2-shaped fit diagnostic column.
#[derive(Clone, Debug, PartialEq)]
pub enum Deseq2McolsDiagnosticValues {
    /// Numeric values, preserving DESeq2-style `NaN` rows.
    Numeric(Vec<f64>),
    /// Optional numeric values, used for nullable Cook's diagnostics.
    OptionalNumeric(Vec<Option<f64>>),
    /// Integer values.
    Integer(Vec<usize>),
    /// Logical values.
    Logical(Vec<bool>),
}

/// One DESeq2-shaped fit diagnostic column.
#[derive(Clone, Debug, PartialEq)]
pub struct Deseq2McolsDiagnosticColumn {
    /// Stable DESeq2-shaped column name.
    pub name: &'static str,
    /// Column values in gene-row order.
    pub values: Deseq2McolsDiagnosticValues,
}

/// Typed DESeq2-shaped diagnostic data-frame view for `mcols(dds)`-style fields.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Deseq2McolsDiagnosticsDataFrame {
    /// Ordered diagnostic columns.
    pub columns: Vec<Deseq2McolsDiagnosticColumn>,
}

/// DESeq2-style row diagnostics derived from an inspectable [`DeseqFit`].
///
/// These names mirror the metadata columns DESeq2 stores in `mcols(dds)` after
/// Wald and LRT pipelines. The Rust fit state keeps more explicit field names;
/// this view exists to make R and parity wrappers straightforward without
/// duplicating data in the core structs.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Deseq2McolsDiagnostics {
    /// Stable fit-type label for the fitted dispersion trend, when present.
    pub dispersion_fit_type: Option<&'static str>,
    /// DESeq2 `dispGeneEst` column from gene-wise dispersion fitting.
    pub disp_gene_est: Option<Vec<f64>>,
    /// DESeq2 `dispGeneIter` column from gene-wise dispersion fitting.
    pub disp_gene_iter: Option<Vec<usize>>,
    /// DESeq2 `dispFit` fitted dispersion trend column.
    pub disp_fit: Option<Vec<f64>>,
    /// DESeq2 `dispersion` final dispersion column.
    pub dispersion: Option<Vec<f64>>,
    /// DESeq2 `dispIter` MAP dispersion iteration column.
    pub disp_iter: Option<Vec<usize>>,
    /// DESeq2 `dispOutlier` MAP dispersion outlier column.
    pub disp_outlier: Option<Vec<bool>>,
    /// Rust convergence flags for implemented MAP dispersion fitting.
    pub dispersion_converged: Option<Vec<bool>>,
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

impl Deseq2McolsDiagnostics {
    /// Stable DESeq2-shaped column names available in this diagnostic view.
    ///
    /// The order follows the usual staged fitting flow: dispersion columns,
    /// GLM convergence/iteration columns, deviance, then Cook's diagnostics.
    pub fn present_column_names(&self) -> Vec<&'static str> {
        let mut names = Vec::new();
        if self.disp_gene_est.is_some() {
            names.push("dispGeneEst");
        }
        if self.disp_gene_iter.is_some() {
            names.push("dispGeneIter");
        }
        if self.disp_fit.is_some() {
            names.push("dispFit");
        }
        if self.dispersion.is_some() {
            names.push("dispersion");
        }
        if self.disp_iter.is_some() {
            names.push("dispIter");
        }
        if self.disp_outlier.is_some() {
            names.push("dispOutlier");
        }
        if self.beta_conv.is_some() {
            names.push("betaConv");
        }
        if self.full_beta_conv.is_some() {
            names.push("fullBetaConv");
        }
        if self.reduced_beta_conv.is_some() {
            names.push("reducedBetaConv");
        }
        if self.beta_iter.is_some() {
            names.push("betaIter");
        }
        if self.reduced_beta_iter.is_some() {
            names.push("reducedBetaIter");
        }
        if self.deviance.is_some() {
            names.push("deviance");
        }
        if self.max_cooks.is_some() {
            names.push("maxCooks");
        }
        names
    }

    /// Assemble present diagnostic fields into a typed DESeq2-shaped data frame.
    pub fn data_frame(&self) -> Deseq2McolsDiagnosticsDataFrame {
        let mut columns = Vec::new();
        if let Some(values) = &self.disp_gene_est {
            columns.push(numeric_diagnostic_column("dispGeneEst", values));
        }
        if let Some(values) = &self.disp_gene_iter {
            columns.push(integer_diagnostic_column("dispGeneIter", values));
        }
        if let Some(values) = &self.disp_fit {
            columns.push(numeric_diagnostic_column("dispFit", values));
        }
        if let Some(values) = &self.dispersion {
            columns.push(numeric_diagnostic_column("dispersion", values));
        }
        if let Some(values) = &self.disp_iter {
            columns.push(integer_diagnostic_column("dispIter", values));
        }
        if let Some(values) = &self.disp_outlier {
            columns.push(logical_diagnostic_column("dispOutlier", values));
        }
        if let Some(values) = &self.beta_conv {
            columns.push(logical_diagnostic_column("betaConv", values));
        }
        if let Some(values) = &self.full_beta_conv {
            columns.push(logical_diagnostic_column("fullBetaConv", values));
        }
        if let Some(values) = &self.reduced_beta_conv {
            columns.push(logical_diagnostic_column("reducedBetaConv", values));
        }
        if let Some(values) = &self.beta_iter {
            columns.push(integer_diagnostic_column("betaIter", values));
        }
        if let Some(values) = &self.reduced_beta_iter {
            columns.push(integer_diagnostic_column("reducedBetaIter", values));
        }
        if let Some(values) = &self.deviance {
            columns.push(numeric_diagnostic_column("deviance", values));
        }
        if let Some(values) = &self.max_cooks {
            columns.push(Deseq2McolsDiagnosticColumn {
                name: "maxCooks",
                values: Deseq2McolsDiagnosticValues::OptionalNumeric(values.clone()),
            });
        }
        Deseq2McolsDiagnosticsDataFrame { columns }
    }
}

fn numeric_diagnostic_column(name: &'static str, values: &[f64]) -> Deseq2McolsDiagnosticColumn {
    Deseq2McolsDiagnosticColumn {
        name,
        values: Deseq2McolsDiagnosticValues::Numeric(values.to_vec()),
    }
}

fn integer_diagnostic_column(name: &'static str, values: &[usize]) -> Deseq2McolsDiagnosticColumn {
    Deseq2McolsDiagnosticColumn {
        name,
        values: Deseq2McolsDiagnosticValues::Integer(values.to_vec()),
    }
}

fn logical_diagnostic_column(name: &'static str, values: &[bool]) -> Deseq2McolsDiagnosticColumn {
    Deseq2McolsDiagnosticColumn {
        name,
        values: Deseq2McolsDiagnosticValues::Logical(values.to_vec()),
    }
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
            dispersion_fit_type: self
                .dispersion_trend
                .as_ref()
                .map(|trend| trend.fit_type_label()),
            disp_gene_est: self.disp_gene_est.clone(),
            disp_gene_iter: self.disp_gene_iter.clone(),
            disp_fit: self.disp_fit.clone(),
            dispersion: self.dispersion.clone(),
            disp_iter: self.disp_iter.clone(),
            disp_outlier: self.disp_outlier.clone(),
            dispersion_converged: self.dispersion_converged.clone(),
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
