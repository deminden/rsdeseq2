use crate::cooks::{
    calculate_cooks_distance, prepare_cooks_replacement_refit, CooksOutput, CooksRefitPlan,
    CooksReplacementOptions,
};
use crate::core::CountMatrix;
use crate::design::{
    expanded_factor_design, expanded_formula_design_with_offsets,
    expanded_formula_design_with_offsets_from_model_frame, DesignMatrix,
    ExpandedAdditiveFactorDesign, ExpandedFactorDesign, ExpandedFactorInteractionSpec,
    ExpandedFactorNumericInteractionSpec, ExpandedFactorSpec, ExpandedFormulaDesignWithOffsets,
    ExpandedNumericInteractionSpec, ExpandedNumericSpec, FormulaModelFrame,
};
use crate::errors::{invalid_dimensions, DeseqError};
use crate::glm::{
    collapse_expanded_model_fit,
    fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights,
    fit_expanded_glms_with_estimated_beta_prior_variance_and_weights, wald_test_coefficient,
    wald_test_contrast, BetaPriorNormalizationFactorWeightInput, BetaPriorRefitOptions,
    BetaPriorSizeFactorWeightInput, ExpandedModelBetaPriorDesignInput,
    ExpandedModelBetaPriorGlmFit, LrtOutput, NbinomGlmFit, WaldAlternative, WaldContrastOutput,
    WaldOutput, WaldTestOptions,
};
use crate::independent_filtering::IndependentFilteringOutput;
use crate::matrix::RowMajorMatrix;
use crate::multiple_testing::bh_adjust;
use crate::normalization::{normalized_counts, normalized_counts_with_factors};
use crate::options::{CooksCutoff, TestType};
use statrs::distribution::{ContinuousCDF, FisherSnedecor};

/// Core DESeq2 `results()` column names emitted by the current Rust result rows.
pub const DESEQ2_RESULT_CORE_COLUMNS: [&str; 6] = [
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
];

/// Optional diagnostic columns carried by [`DeseqResultRow`] when present.
pub const RSDESEQ2_RESULT_DIAGNOSTIC_COLUMNS: [&str; 5] = [
    "dispersion",
    "converged",
    "maxCooks",
    "cooksOutlier",
    "filtered",
];

/// Metadata for one DESeq2-style result-table column.
#[derive(Clone, Debug, PartialEq)]
pub struct DeseqResultColumnMetadata {
    /// Column name.
    pub name: String,
    /// DESeq2-style column group. Core statistical columns use `results`.
    pub column_type: String,
    /// Human-readable column description.
    pub description: String,
}

/// Table-level metadata for a primitive DESeq2-shaped result table.
#[derive(Clone, Debug, PartialEq)]
pub struct DeseqResultsTableMetadata {
    /// Statistical test represented by the p-values.
    pub test_type: Option<TestType>,
    /// Name of the reported coefficient or contrast, if known.
    pub result_name: Option<String>,
    /// Free-form comparison description for wrappers or callers.
    pub comparison: Option<String>,
    /// Resolved numeric contrast over the fitted model coefficients, if this
    /// table reports a contrast rather than a single coefficient.
    pub contrast: Option<Vec<f64>>,
    /// Log2 fold-change threshold used for Wald-threshold tests.
    pub lfc_threshold: f64,
    /// Alternative hypothesis name for thresholded Wald tests.
    pub alt_hypothesis: Option<String>,
    /// P-value adjustment method. The current Rust implementation uses BH.
    pub p_adjust_method: String,
}

/// One scalar table-level result metadata entry.
#[derive(Clone, Debug, PartialEq)]
pub struct DeseqResultsTableMetadataEntry {
    /// Metadata key.
    pub name: String,
    /// Metadata value formatted for table/export consumers.
    pub value: String,
}

impl Default for DeseqResultsTableMetadata {
    fn default() -> Self {
        Self {
            test_type: None,
            result_name: None,
            comparison: None,
            contrast: None,
            lfc_threshold: 0.0,
            alt_hypothesis: None,
            p_adjust_method: "BH".to_string(),
        }
    }
}

impl DeseqResultsTableMetadata {
    /// Label used in effect-size column descriptions.
    ///
    /// Wald contrast result tables prefer the comparison label when present,
    /// while LRT result tables keep the reported full-model coefficient as the
    /// effect label and use the model comparison for test-statistic columns.
    pub fn effect_description_label(&self) -> Option<&str> {
        effect_description_label(self)
    }

    /// Label used in test-statistic and p-value column descriptions.
    pub fn test_description_label(&self) -> Option<&str> {
        test_description_label(self)
    }

    /// Assemble scalar table-level metadata entries.
    ///
    /// The names are stable Rust-side keys for wrapper and file exporters,
    /// while values follow DESeq2-facing labels where applicable.
    pub fn scalar_metadata(&self) -> Vec<DeseqResultsTableMetadataEntry> {
        let mut entries = Vec::new();
        if let Some(test_type) = self.test_type {
            entries.push(DeseqResultsTableMetadataEntry {
                name: "testType".to_string(),
                value: test_type_label(test_type).to_string(),
            });
        }
        if let Some(value) = &self.result_name {
            entries.push(DeseqResultsTableMetadataEntry {
                name: "resultName".to_string(),
                value: value.clone(),
            });
        }
        if let Some(value) = &self.comparison {
            entries.push(DeseqResultsTableMetadataEntry {
                name: "comparison".to_string(),
                value: value.clone(),
            });
        }
        if let Some(values) = &self.contrast {
            entries.push(DeseqResultsTableMetadataEntry {
                name: "contrast".to_string(),
                value: values
                    .iter()
                    .map(|value| value.to_string())
                    .collect::<Vec<_>>()
                    .join(","),
            });
        }
        entries.push(DeseqResultsTableMetadataEntry {
            name: "lfcThreshold".to_string(),
            value: self.lfc_threshold.to_string(),
        });
        if let Some(value) = &self.alt_hypothesis {
            entries.push(DeseqResultsTableMetadataEntry {
                name: "altHypothesis".to_string(),
                value: value.clone(),
            });
        }
        entries.push(DeseqResultsTableMetadataEntry {
            name: "pAdjustMethod".to_string(),
            value: self.p_adjust_method.clone(),
        });
        entries
    }
}

/// Combined metadata view for a result table.
#[derive(Clone, Debug, PartialEq)]
pub struct DeseqResultsMetadata {
    /// Table-level metadata.
    pub table: DeseqResultsTableMetadata,
    /// Column metadata for currently represented columns.
    pub columns: Vec<DeseqResultColumnMetadata>,
    /// Independent-filtering metadata, if filtering has run.
    pub independent_filtering: Option<IndependentFilteringOutput>,
}

/// Typed values for one assembled result-table column.
#[derive(Clone, Debug, PartialEq)]
pub enum DeseqResultColumnValues {
    /// Numeric column values. Missing or non-finite values are represented as `None`.
    Numeric(Vec<Option<f64>>),
    /// Logical column values. Missing values are represented as `None`.
    Logical(Vec<Option<bool>>),
}

impl DeseqResultColumnValues {
    /// Number of values in the column.
    pub fn len(&self) -> usize {
        match self {
            Self::Numeric(values) => values.len(),
            Self::Logical(values) => values.len(),
        }
    }

    /// Whether the column has no values.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Borrow numeric values when this is a numeric column.
    pub fn as_numeric(&self) -> Option<&[Option<f64>]> {
        match self {
            Self::Numeric(values) => Some(values),
            Self::Logical(_) => None,
        }
    }

    /// Borrow logical values when this is a logical column.
    pub fn as_logical(&self) -> Option<&[Option<bool>]> {
        match self {
            Self::Numeric(_) => None,
            Self::Logical(values) => Some(values),
        }
    }
}

/// One assembled result-table column with values and DESeq2-style metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct DeseqResultColumn {
    /// Metadata describing the column.
    pub metadata: DeseqResultColumnMetadata,
    /// Column values in result-row order.
    pub values: DeseqResultColumnValues,
}

/// DESeq2-shaped typed data-frame view of a result table.
#[derive(Clone, Debug, PartialEq)]
pub struct DeseqResultsDataFrame {
    /// Optional row names, typically gene identifiers.
    pub row_names: Vec<Option<String>>,
    /// Ordered result columns.
    pub columns: Vec<DeseqResultColumn>,
    /// Table and independent-filtering metadata.
    pub metadata: DeseqResultsMetadata,
}

/// Inputs for a primitive expanded beta-prior Wald fit-and-results workflow.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedBetaPriorWaldResultsInput<'a> {
    /// Raw count matrix.
    pub counts: &'a CountMatrix,
    /// Expanded and reported design matrices plus coefficient collapse groups.
    pub design: ExpandedModelBetaPriorDesignInput<'a>,
    /// Per-sample size factors.
    pub size_factors: &'a [f64],
    /// Optional normalized observation weights.
    pub weights: Option<&'a RowMajorMatrix<f64>>,
    /// Per-gene final dispersions used by the fixed-dispersion GLM.
    pub dispersions: &'a [f64],
    /// Per-gene base means used for result rows and beta-prior weights.
    pub base_mean: &'a [f64],
    /// Per-gene fitted dispersion trend used for beta-prior weights.
    pub disp_fit: &'a [f64],
    /// Optional gene names for result rows.
    pub gene_names: Option<&'a [String]>,
    /// Beta-prior refit options.
    pub options: BetaPriorRefitOptions,
}

/// Inputs for an expanded beta-prior Wald workflow with normalization factors.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedBetaPriorWaldNormalizedResultsInput<'a> {
    /// Raw count matrix.
    pub counts: &'a CountMatrix,
    /// Expanded and reported design matrices plus coefficient collapse groups.
    pub design: ExpandedModelBetaPriorDesignInput<'a>,
    /// Gene x sample normalization-factor matrix.
    pub normalization_factors: &'a RowMajorMatrix<f64>,
    /// Optional normalized observation weights.
    pub weights: Option<&'a RowMajorMatrix<f64>>,
    /// Per-gene final dispersions used by the fixed-dispersion GLM.
    pub dispersions: &'a [f64],
    /// Per-gene base means used for result rows and beta-prior weights.
    pub base_mean: &'a [f64],
    /// Per-gene fitted dispersion trend used for beta-prior weights.
    pub disp_fit: &'a [f64],
    /// Optional gene names for result rows.
    pub gene_names: Option<&'a [String]>,
    /// Beta-prior refit options.
    pub options: BetaPriorRefitOptions,
}

/// Inputs for a one-factor expanded beta-prior Wald workflow.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedFactorBetaPriorWaldResultsInput<'a> {
    /// Raw count matrix.
    pub counts: &'a CountMatrix,
    /// Factor name used to build coefficient names.
    pub factor: &'a str,
    /// Per-sample factor levels in count-column order.
    pub sample_levels: &'a [String],
    /// Reference level for treatment-style reported coefficients.
    pub reference: &'a str,
    /// Per-sample size factors.
    pub size_factors: &'a [f64],
    /// Optional normalized observation weights.
    pub weights: Option<&'a RowMajorMatrix<f64>>,
    /// Per-gene final dispersions used by the fixed-dispersion GLM.
    pub dispersions: &'a [f64],
    /// Per-gene base means used for result rows and beta-prior weights.
    pub base_mean: &'a [f64],
    /// Per-gene fitted dispersion trend used for beta-prior weights.
    pub disp_fit: &'a [f64],
    /// Optional gene names for result rows.
    pub gene_names: Option<&'a [String]>,
    /// Beta-prior refit options.
    pub options: BetaPriorRefitOptions,
}

/// Inputs for a one-factor expanded beta-prior Wald workflow with normalization factors.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedFactorBetaPriorWaldNormalizedResultsInput<'a> {
    /// Raw count matrix.
    pub counts: &'a CountMatrix,
    /// Factor name used to build coefficient names.
    pub factor: &'a str,
    /// Per-sample factor levels in count-column order.
    pub sample_levels: &'a [String],
    /// Reference level for treatment-style reported coefficients.
    pub reference: &'a str,
    /// Gene x sample normalization-factor matrix.
    pub normalization_factors: &'a RowMajorMatrix<f64>,
    /// Optional normalized observation weights.
    pub weights: Option<&'a RowMajorMatrix<f64>>,
    /// Per-gene final dispersions used by the fixed-dispersion GLM.
    pub dispersions: &'a [f64],
    /// Per-gene base means used for result rows and beta-prior weights.
    pub base_mean: &'a [f64],
    /// Per-gene fitted dispersion trend used for beta-prior weights.
    pub disp_fit: &'a [f64],
    /// Optional gene names for result rows.
    pub gene_names: Option<&'a [String]>,
    /// Beta-prior refit options.
    pub options: BetaPriorRefitOptions,
}

/// Inputs for an additive-factor expanded beta-prior Wald workflow.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedAdditiveBetaPriorWaldResultsInput<'a> {
    /// Raw count matrix.
    pub counts: &'a CountMatrix,
    /// Additive factor terms used to build expanded and reported designs.
    pub factors: &'a [ExpandedFactorSpec<'a>],
    /// Additive numeric covariates included unchanged in both design matrices.
    pub numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    /// Factor-by-factor interactions included after main effects.
    pub interactions: &'a [ExpandedFactorInteractionSpec<'a>],
    /// Factor-by-numeric interactions included after factor-by-factor interactions.
    pub factor_numeric_interactions: &'a [ExpandedFactorNumericInteractionSpec<'a>],
    /// Numeric-by-numeric interactions included after factor-by-numeric interactions.
    pub numeric_interactions: &'a [ExpandedNumericInteractionSpec<'a>],
    /// Per-sample size factors.
    pub size_factors: &'a [f64],
    /// Optional normalized observation weights.
    pub weights: Option<&'a RowMajorMatrix<f64>>,
    /// Per-gene final dispersions used by the fixed-dispersion GLM.
    pub dispersions: &'a [f64],
    /// Per-gene base means used for result rows and beta-prior weights.
    pub base_mean: &'a [f64],
    /// Per-gene fitted dispersion trend used for beta-prior weights.
    pub disp_fit: &'a [f64],
    /// Optional gene names for result rows.
    pub gene_names: Option<&'a [String]>,
    /// Beta-prior refit options.
    pub options: BetaPriorRefitOptions,
}

/// Inputs for an additive-factor expanded beta-prior Wald workflow with normalization factors.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedAdditiveBetaPriorWaldNormalizedResultsInput<'a> {
    /// Raw count matrix.
    pub counts: &'a CountMatrix,
    /// Additive factor terms used to build expanded and reported designs.
    pub factors: &'a [ExpandedFactorSpec<'a>],
    /// Additive numeric covariates included unchanged in both design matrices.
    pub numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    /// Factor-by-factor interactions included after main effects.
    pub interactions: &'a [ExpandedFactorInteractionSpec<'a>],
    /// Factor-by-numeric interactions included after factor-by-factor interactions.
    pub factor_numeric_interactions: &'a [ExpandedFactorNumericInteractionSpec<'a>],
    /// Numeric-by-numeric interactions included after factor-by-numeric interactions.
    pub numeric_interactions: &'a [ExpandedNumericInteractionSpec<'a>],
    /// Gene x sample normalization-factor matrix.
    pub normalization_factors: &'a RowMajorMatrix<f64>,
    /// Optional normalized observation weights.
    pub weights: Option<&'a RowMajorMatrix<f64>>,
    /// Per-gene final dispersions used by the fixed-dispersion GLM.
    pub dispersions: &'a [f64],
    /// Per-gene base means used for result rows and beta-prior weights.
    pub base_mean: &'a [f64],
    /// Per-gene fitted dispersion trend used for beta-prior weights.
    pub disp_fit: &'a [f64],
    /// Optional gene names for result rows.
    pub gene_names: Option<&'a [String]>,
    /// Beta-prior refit options.
    pub options: BetaPriorRefitOptions,
}

/// Inputs for a formula-driven expanded beta-prior Wald workflow.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedFormulaBetaPriorWaldResultsInput<'a> {
    /// Raw count matrix.
    pub counts: &'a CountMatrix,
    /// Primitive formula parsed by [`expanded_formula_design`].
    pub formula: &'a str,
    /// Candidate factor metadata referenced by the formula.
    pub factors: &'a [ExpandedFactorSpec<'a>],
    /// Candidate numeric covariates referenced by the formula.
    pub numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    /// Per-sample size factors.
    pub size_factors: &'a [f64],
    /// Optional normalized observation weights.
    pub weights: Option<&'a RowMajorMatrix<f64>>,
    /// Per-gene final dispersions used by the fixed-dispersion GLM.
    pub dispersions: &'a [f64],
    /// Per-gene base means used for result rows and beta-prior weights.
    pub base_mean: &'a [f64],
    /// Per-gene fitted dispersion trend used for beta-prior weights.
    pub disp_fit: &'a [f64],
    /// Optional gene names for result rows.
    pub gene_names: Option<&'a [String]>,
    /// Beta-prior refit options.
    pub options: BetaPriorRefitOptions,
}

/// Inputs for a formula-driven expanded beta-prior Wald workflow using owned
/// model-frame sample metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedFormulaModelFrameBetaPriorWaldResultsInput<'a> {
    /// Raw count matrix.
    pub counts: &'a CountMatrix,
    /// Primitive formula parsed by [`expanded_formula_design_from_model_frame`].
    pub formula: &'a str,
    /// Owned factor and numeric covariates referenced by the formula.
    pub model_frame: &'a FormulaModelFrame,
    /// Per-sample size factors.
    pub size_factors: &'a [f64],
    /// Optional normalized observation weights.
    pub weights: Option<&'a RowMajorMatrix<f64>>,
    /// Per-gene final dispersions used by the fixed-dispersion GLM.
    pub dispersions: &'a [f64],
    /// Per-gene base means used for result rows and beta-prior weights.
    pub base_mean: &'a [f64],
    /// Per-gene fitted dispersion trend used for beta-prior weights.
    pub disp_fit: &'a [f64],
    /// Optional gene names for result rows.
    pub gene_names: Option<&'a [String]>,
    /// Beta-prior refit options.
    pub options: BetaPriorRefitOptions,
}

/// Inputs for a formula-driven expanded beta-prior Wald workflow with normalization factors.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedFormulaBetaPriorWaldNormalizedResultsInput<'a> {
    /// Raw count matrix.
    pub counts: &'a CountMatrix,
    /// Primitive formula parsed by [`expanded_formula_design`].
    pub formula: &'a str,
    /// Candidate factor metadata referenced by the formula.
    pub factors: &'a [ExpandedFactorSpec<'a>],
    /// Candidate numeric covariates referenced by the formula.
    pub numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    /// Gene x sample normalization-factor matrix.
    pub normalization_factors: &'a RowMajorMatrix<f64>,
    /// Optional normalized observation weights.
    pub weights: Option<&'a RowMajorMatrix<f64>>,
    /// Per-gene final dispersions used by the fixed-dispersion GLM.
    pub dispersions: &'a [f64],
    /// Per-gene base means used for result rows and beta-prior weights.
    pub base_mean: &'a [f64],
    /// Per-gene fitted dispersion trend used for beta-prior weights.
    pub disp_fit: &'a [f64],
    /// Optional gene names for result rows.
    pub gene_names: Option<&'a [String]>,
    /// Beta-prior refit options.
    pub options: BetaPriorRefitOptions,
}

/// Inputs for a formula-driven expanded beta-prior Wald workflow with
/// normalization factors and owned model-frame sample metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedFormulaModelFrameBetaPriorWaldNormalizedResultsInput<'a> {
    /// Raw count matrix.
    pub counts: &'a CountMatrix,
    /// Primitive formula parsed by [`expanded_formula_design_from_model_frame`].
    pub formula: &'a str,
    /// Owned factor and numeric covariates referenced by the formula.
    pub model_frame: &'a FormulaModelFrame,
    /// Gene x sample normalization-factor matrix.
    pub normalization_factors: &'a RowMajorMatrix<f64>,
    /// Optional normalized observation weights.
    pub weights: Option<&'a RowMajorMatrix<f64>>,
    /// Per-gene final dispersions used by the fixed-dispersion GLM.
    pub dispersions: &'a [f64],
    /// Per-gene base means used for result rows and beta-prior weights.
    pub base_mean: &'a [f64],
    /// Per-gene fitted dispersion trend used for beta-prior weights.
    pub disp_fit: &'a [f64],
    /// Optional gene names for result rows.
    pub gene_names: Option<&'a [String]>,
    /// Beta-prior refit options.
    pub options: BetaPriorRefitOptions,
}

/// Expanded beta-prior fit plus DESeq2-shaped Wald result rows.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedBetaPriorWaldResults {
    /// Expanded-design beta-prior fit with collapsed standard-design prior fit.
    pub fit: ExpandedModelBetaPriorGlmFit,
    /// Wald result table built from the collapsed prior fit.
    pub results: DeseqResults,
}

/// Expanded beta-prior Wald result workflow with Cook's replacement-refit metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedBetaPriorWaldReplacementResults {
    /// Original beta-prior fit and result rows before count replacement.
    pub original: ExpandedBetaPriorWaldResults,
    /// Cook's distances calculated from the original collapsed prior fit.
    pub cooks: CooksOutput,
    /// Count-replacement and refit plan derived from original Cook's distances.
    pub refit_plan: CooksRefitPlan,
    /// Optional beta-prior refit on replacement counts.
    pub refit: Option<ExpandedBetaPriorWaldResults>,
    /// Final merged result rows after replacement refit and Cook's filtering.
    pub results: DeseqResults,
}

/// One-factor expanded beta-prior design, fit, and DESeq2-shaped Wald rows.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedFactorBetaPriorWaldResults {
    /// Generated expanded and standard one-factor design matrices.
    pub design: ExpandedFactorDesign,
    /// Expanded-design beta-prior fit with collapsed standard-design prior fit.
    pub fit: ExpandedModelBetaPriorGlmFit,
    /// Wald result table built from the collapsed prior fit.
    pub results: DeseqResults,
}

/// One-factor expanded beta-prior Wald replacement workflow with generated design metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedFactorBetaPriorWaldReplacementResults {
    /// Generated expanded and standard one-factor design matrices.
    pub design: ExpandedFactorDesign,
    /// Replacement-refit workflow output for the generated design.
    pub replacement: ExpandedBetaPriorWaldReplacementResults,
}

/// Additive-factor expanded beta-prior design, fit, and DESeq2-shaped Wald rows.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedAdditiveBetaPriorWaldResults {
    /// Generated expanded and standard additive-factor design matrices.
    pub design: ExpandedAdditiveFactorDesign,
    /// Expanded-design beta-prior fit with collapsed standard-design prior fit.
    pub fit: ExpandedModelBetaPriorGlmFit,
    /// Wald result table built from the collapsed prior fit.
    pub results: DeseqResults,
}

/// Additive expanded beta-prior Wald replacement workflow with generated design metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedAdditiveBetaPriorWaldReplacementResults {
    /// Generated expanded and standard additive-factor design matrices.
    pub design: ExpandedAdditiveFactorDesign,
    /// Replacement-refit workflow output for the generated design.
    pub replacement: ExpandedBetaPriorWaldReplacementResults,
}

/// One row of a DESeq2-like results table.
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
    /// Table-level metadata, such as test type and reported coefficient.
    pub metadata: DeseqResultsTableMetadata,
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

    /// Column names represented by this result table.
    ///
    /// The six core names match the simple DESeq2 `results()` table. Optional
    /// Rust diagnostics are included only when at least one row carries that
    /// field. Gene identifiers are represented as row names by R-style
    /// frontends, not as a result column.
    pub fn column_names(&self) -> Vec<&'static str> {
        let mut columns = deseq2_result_core_column_names().to_vec();
        if self.rows.iter().any(|row| row.dispersion.is_some()) {
            columns.push("dispersion");
        }
        if self.rows.iter().any(|row| row.converged.is_some()) {
            columns.push("converged");
        }
        if self.rows.iter().any(|row| row.max_cooks.is_some()) {
            columns.push("maxCooks");
        }
        if self.rows.iter().any(|row| row.cooks_outlier.is_some()) {
            columns.push("cooksOutlier");
        }
        if self.rows.iter().any(|row| row.filtered.is_some()) {
            columns.push("filtered");
        }
        columns
    }

    /// Column metadata represented by this result table.
    pub fn column_metadata(&self) -> Vec<DeseqResultColumnMetadata> {
        self.column_names()
            .into_iter()
            .map(|name| DeseqResultColumnMetadata {
                name: name.to_string(),
                column_type: result_column_type(name).to_string(),
                description: result_column_description(name, &self.metadata),
            })
            .collect()
    }

    /// DESeq2-style metadata view for table and represented columns.
    pub fn deseq2_metadata(&self) -> DeseqResultsMetadata {
        DeseqResultsMetadata {
            table: self.metadata.clone(),
            columns: self.column_metadata(),
            independent_filtering: self.independent_filtering.clone(),
        }
    }

    /// Assemble a typed DESeq2-shaped data-frame view.
    ///
    /// Row names are kept separate from result columns, matching R's
    /// `DataFrame`/data-frame convention. Columns follow [`Self::column_names`]
    /// and carry the same metadata returned by [`Self::deseq2_metadata`].
    pub fn data_frame(&self) -> DeseqResultsDataFrame {
        let column_metadata = self.column_metadata();
        let columns = column_metadata
            .iter()
            .map(|metadata| DeseqResultColumn {
                metadata: metadata.clone(),
                values: self.column_values(&metadata.name),
            })
            .collect();
        DeseqResultsDataFrame {
            row_names: self.rows.iter().map(|row| row.gene.clone()).collect(),
            columns,
            metadata: DeseqResultsMetadata {
                table: self.metadata.clone(),
                columns: column_metadata,
                independent_filtering: self.independent_filtering.clone(),
            },
        }
    }

    /// Return a copy with updated table metadata.
    pub fn with_metadata(mut self, metadata: DeseqResultsTableMetadata) -> Self {
        self.metadata = metadata;
        self
    }

    /// Attach metadata for a resolved numeric contrast.
    pub fn set_resolved_contrast_metadata(
        &mut self,
        result_name: impl Into<String>,
        comparison: impl Into<String>,
        contrast: &[f64],
    ) {
        self.metadata.result_name = Some(result_name.into());
        self.metadata.comparison = Some(comparison.into());
        self.metadata.contrast = Some(contrast.to_vec());
    }

    /// Attach the Wald threshold and alternative used to produce p-values.
    pub fn apply_wald_test_options(&mut self, options: &WaldTestOptions) {
        self.metadata.lfc_threshold = options.lfc_threshold;
        self.metadata.alt_hypothesis = Some(wald_alternative_name(options.alternative).to_string());
    }

    /// Return a copy with Wald threshold metadata attached.
    pub fn with_wald_test_options(mut self, options: &WaldTestOptions) -> Self {
        self.apply_wald_test_options(options);
        self
    }

    fn column_values(&self, name: &str) -> DeseqResultColumnValues {
        match name {
            "baseMean" => numeric_column(&self.rows, |row| finite_option(row.base_mean)),
            "log2FoldChange" => numeric_column(&self.rows, |row| row.log2_fold_change),
            "lfcSE" => numeric_column(&self.rows, |row| row.lfc_se),
            "stat" => numeric_column(&self.rows, |row| row.stat),
            "pvalue" => numeric_column(&self.rows, |row| row.pvalue),
            "padj" => numeric_column(&self.rows, |row| row.padj),
            "dispersion" => numeric_column(&self.rows, |row| row.dispersion),
            "maxCooks" => numeric_column(&self.rows, |row| row.max_cooks),
            "converged" => logical_column(&self.rows, |row| row.converged),
            "cooksOutlier" => logical_column(&self.rows, |row| row.cooks_outlier),
            "filtered" => logical_column(&self.rows, |row| row.filtered),
            _ => numeric_column(&self.rows, |_| None),
        }
    }
}
