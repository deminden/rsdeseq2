use crate::cooks::{
    calculate_cooks_distance, prepare_cooks_replacement_refit, CooksOutput, CooksRefitPlan,
    CooksReplacementOptions,
};
use crate::core::CountMatrix;
use crate::design::{
    expanded_factor_design, expanded_formula_design_with_offsets, DesignMatrix,
    ExpandedAdditiveFactorDesign, ExpandedFactorDesign, ExpandedFactorInteractionSpec,
    ExpandedFactorNumericInteractionSpec, ExpandedFactorSpec, ExpandedNumericInteractionSpec,
    ExpandedNumericSpec,
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
    /// Expanded and reported design surfaces plus coefficient collapse groups.
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
    /// Expanded and reported design surfaces plus coefficient collapse groups.
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
    /// Additive numeric covariates included unchanged in both design surfaces.
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
    /// Additive numeric covariates included unchanged in both design surfaces.
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
    /// Generated expanded/standard one-factor design surfaces.
    pub design: ExpandedFactorDesign,
    /// Expanded-design beta-prior fit with collapsed standard-design prior fit.
    pub fit: ExpandedModelBetaPriorGlmFit,
    /// Wald result table built from the collapsed prior fit.
    pub results: DeseqResults,
}

/// One-factor expanded beta-prior Wald replacement workflow with generated design metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedFactorBetaPriorWaldReplacementResults {
    /// Generated expanded/standard one-factor design surfaces.
    pub design: ExpandedFactorDesign,
    /// Replacement-refit workflow output for the generated design.
    pub replacement: ExpandedBetaPriorWaldReplacementResults,
}

/// Additive-factor expanded beta-prior design, fit, and DESeq2-shaped Wald rows.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedAdditiveBetaPriorWaldResults {
    /// Generated expanded/standard additive-factor design surfaces.
    pub design: ExpandedAdditiveFactorDesign,
    /// Expanded-design beta-prior fit with collapsed standard-design prior fit.
    pub fit: ExpandedModelBetaPriorGlmFit,
    /// Wald result table built from the collapsed prior fit.
    pub results: DeseqResults,
}

/// Additive expanded beta-prior Wald replacement workflow with generated design metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedAdditiveBetaPriorWaldReplacementResults {
    /// Generated expanded/standard additive-factor design surfaces.
    pub design: ExpandedAdditiveFactorDesign,
    /// Replacement-refit workflow output for the generated design.
    pub replacement: ExpandedBetaPriorWaldReplacementResults,
}

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

/// Return the core DESeq2 `results()` column names currently emitted by Rust.
pub fn deseq2_result_core_column_names() -> &'static [&'static str] {
    &DESEQ2_RESULT_CORE_COLUMNS
}

/// Return optional diagnostic result-column names used by Rust result rows.
pub fn rsdeseq2_result_diagnostic_column_names() -> &'static [&'static str] {
    &RSDESEQ2_RESULT_DIAGNOSTIC_COLUMNS
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

/// Collapse an expanded-model fit and build DESeq2-shaped Wald results.
///
/// This is a primitive result-table companion for the beta-prior expanded
/// model workflow. It performs grouped coefficient/covariance collapse, then
/// reports the requested standard-design coefficient with ordinary Wald
/// statistics and BH adjustment.
pub fn build_wald_results_from_expanded_model_fit(
    base_mean: &[f64],
    expanded_fit: &NbinomGlmFit,
    standard_design: &DesignMatrix,
    coefficient_groups: &[Vec<usize>],
    coefficient: usize,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    let collapsed = collapse_expanded_model_fit(expanded_fit, standard_design, coefficient_groups)?;
    build_wald_results(base_mean, &collapsed, coefficient, gene_names, dispersions)
}

/// Collapse an expanded-model fit and build DESeq2-shaped Wald contrast rows.
///
/// This is the contrast companion to
/// [`build_wald_results_from_expanded_model_fit`]. The supplied contrast is on
/// the collapsed standard-design coefficient scale; the helper propagates the
/// expanded covariance through the grouped coefficient average before computing
/// `c' beta` and `sqrt(c' Sigma c)`.
pub fn build_wald_contrast_results_from_expanded_model_fit(
    base_mean: &[f64],
    expanded_fit: &NbinomGlmFit,
    standard_design: &DesignMatrix,
    coefficient_groups: &[Vec<usize>],
    contrast: &[f64],
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    let collapsed = collapse_expanded_model_fit(expanded_fit, standard_design, coefficient_groups)?;
    let contrast = wald_test_contrast(&collapsed, contrast)?;
    build_wald_contrast_results(base_mean, &collapsed, &contrast, gene_names, dispersions)
}

/// Build DESeq2-shaped Wald rows from an expanded beta-prior refit output.
///
/// The helper reports the already-collapsed standard-design prior fit stored in
/// [`ExpandedModelBetaPriorGlmFit`], so callers that use the expanded beta-prior
/// workflow do not need to manually pass the collapsed fit to result assembly.
pub fn build_wald_results_from_expanded_beta_prior_fit(
    base_mean: &[f64],
    fit: &ExpandedModelBetaPriorGlmFit,
    coefficient: usize,
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    build_wald_results(
        base_mean,
        &fit.prior_fit,
        coefficient,
        gene_names,
        dispersions,
    )
}

/// Build DESeq2-shaped Wald contrast rows from an expanded beta-prior refit output.
///
/// The supplied contrast is on the collapsed standard-design coefficient scale.
pub fn build_wald_contrast_results_from_expanded_beta_prior_fit(
    base_mean: &[f64],
    fit: &ExpandedModelBetaPriorGlmFit,
    contrast: &[f64],
    gene_names: Option<&[String]>,
    dispersions: Option<&[f64]>,
) -> Result<DeseqResults, DeseqError> {
    let contrast = wald_test_contrast(&fit.prior_fit, contrast)?;
    build_wald_contrast_results(
        base_mean,
        &fit.prior_fit,
        &contrast,
        gene_names,
        dispersions,
    )
}

/// Fit an expanded beta-prior model and assemble Wald rows for one coefficient.
///
/// This is a primitive all-Rust companion for callers that already provide the
/// expanded design, standard design, and coefficient groups.
pub fn fit_expanded_beta_prior_wald_results(
    input: ExpandedBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedBetaPriorWaldResults, DeseqError> {
    let fit = fit_expanded_glms_with_estimated_beta_prior_variance_and_weights(
        input.counts,
        input.design,
        BetaPriorSizeFactorWeightInput {
            size_factors: input.size_factors,
            weights: input.weights,
        },
        input.dispersions,
        input.base_mean,
        input.disp_fit,
        input.options,
    )?;
    let results = build_wald_results_from_expanded_beta_prior_fit(
        input.base_mean,
        &fit,
        coefficient,
        input.gene_names,
        Some(input.dispersions),
    )?;
    Ok(ExpandedBetaPriorWaldResults { fit, results })
}

/// Fit an expanded beta-prior model and assemble Wald rows for a numeric contrast.
///
/// The contrast is on the collapsed standard-design coefficient scale.
pub fn fit_expanded_beta_prior_wald_contrast_results(
    input: ExpandedBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedBetaPriorWaldResults, DeseqError> {
    let fit = fit_expanded_glms_with_estimated_beta_prior_variance_and_weights(
        input.counts,
        input.design,
        BetaPriorSizeFactorWeightInput {
            size_factors: input.size_factors,
            weights: input.weights,
        },
        input.dispersions,
        input.base_mean,
        input.disp_fit,
        input.options,
    )?;
    let results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        input.base_mean,
        &fit,
        contrast,
        input.gene_names,
        Some(input.dispersions),
    )?;
    Ok(ExpandedBetaPriorWaldResults { fit, results })
}

/// Fit an expanded beta-prior Wald coefficient workflow with Cook's replacement refit.
///
/// Cook's distances are calculated from the collapsed prior fit on the reported
/// standard-design surface. Replacement counts are then refit through the same
/// expanded beta-prior workflow with the original size factors and supplied
/// dispersions.
pub fn fit_expanded_beta_prior_wald_results_with_cooks_replacement(
    input: ExpandedBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedBetaPriorWaldReplacementResults, DeseqError> {
    let original = fit_expanded_beta_prior_wald_results(input.clone(), coefficient)?;
    let cooks = beta_prior_cooks_output(input.counts, input.size_factors, &original.fit)?;
    let normalized = normalized_counts(input.counts, input.size_factors)?;
    let refit_plan = prepare_cooks_replacement_refit(
        input.counts,
        &normalized,
        input.size_factors,
        None,
        &cooks.cooks,
        input.design.standard_design,
        replacement_options,
    )?;

    let refit = if refit_plan.should_refit {
        Some(fit_expanded_beta_prior_wald_results(
            ExpandedBetaPriorWaldResultsInput {
                counts: &refit_plan.replacement.replaced_counts,
                design: input.design,
                size_factors: input.size_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: &refit_plan.replaced_base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            coefficient,
        )?)
    } else {
        None
    };

    let mut original_results = original.results.clone();
    attach_cooks_to_results(&mut original_results, &cooks.max_cooks)?;
    let mut results = merge_beta_prior_replacement_results(
        &original_results,
        refit.as_ref().map(|value| &value.results),
        &refit_plan,
    )?;
    apply_cooks_cutoff(&mut results, Some(replacement_options.cooks_cutoff))?;

    Ok(ExpandedBetaPriorWaldReplacementResults {
        original: ExpandedBetaPriorWaldResults {
            fit: original.fit,
            results: original_results,
        },
        cooks,
        refit_plan,
        refit,
        results,
    })
}

/// Fit an expanded beta-prior Wald contrast workflow with Cook's replacement refit.
pub fn fit_expanded_beta_prior_wald_contrast_results_with_cooks_replacement(
    input: ExpandedBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedBetaPriorWaldReplacementResults, DeseqError> {
    let original = fit_expanded_beta_prior_wald_contrast_results(input.clone(), contrast)?;
    let cooks = beta_prior_cooks_output(input.counts, input.size_factors, &original.fit)?;
    let normalized = normalized_counts(input.counts, input.size_factors)?;
    let refit_plan = prepare_cooks_replacement_refit(
        input.counts,
        &normalized,
        input.size_factors,
        None,
        &cooks.cooks,
        input.design.standard_design,
        replacement_options,
    )?;

    let refit = if refit_plan.should_refit {
        Some(fit_expanded_beta_prior_wald_contrast_results(
            ExpandedBetaPriorWaldResultsInput {
                counts: &refit_plan.replacement.replaced_counts,
                design: input.design,
                size_factors: input.size_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: &refit_plan.replaced_base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            contrast,
        )?)
    } else {
        None
    };

    let mut original_results = original.results.clone();
    attach_cooks_to_results(&mut original_results, &cooks.max_cooks)?;
    let mut results = merge_beta_prior_replacement_results(
        &original_results,
        refit.as_ref().map(|value| &value.results),
        &refit_plan,
    )?;
    apply_cooks_cutoff(&mut results, Some(replacement_options.cooks_cutoff))?;

    Ok(ExpandedBetaPriorWaldReplacementResults {
        original: ExpandedBetaPriorWaldResults {
            fit: original.fit,
            results: original_results,
        },
        cooks,
        refit_plan,
        refit,
        results,
    })
}

/// Fit an expanded beta-prior model with normalization factors and assemble Wald rows.
pub fn fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
    input: ExpandedBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedBetaPriorWaldResults, DeseqError> {
    let fit =
        fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
            input.counts,
            input.design,
            BetaPriorNormalizationFactorWeightInput {
                normalization_factors: input.normalization_factors,
                weights: input.weights,
            },
            input.dispersions,
            input.base_mean,
            input.disp_fit,
            input.options,
        )?;
    let results = build_wald_results_from_expanded_beta_prior_fit(
        input.base_mean,
        &fit,
        coefficient,
        input.gene_names,
        Some(input.dispersions),
    )?;
    Ok(ExpandedBetaPriorWaldResults { fit, results })
}

/// Fit an expanded beta-prior model with normalization factors and assemble contrast rows.
pub fn fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
    input: ExpandedBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedBetaPriorWaldResults, DeseqError> {
    let fit =
        fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
            input.counts,
            input.design,
            BetaPriorNormalizationFactorWeightInput {
                normalization_factors: input.normalization_factors,
                weights: input.weights,
            },
            input.dispersions,
            input.base_mean,
            input.disp_fit,
            input.options,
        )?;
    let results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        input.base_mean,
        &fit,
        contrast,
        input.gene_names,
        Some(input.dispersions),
    )?;
    Ok(ExpandedBetaPriorWaldResults { fit, results })
}

/// Fit a normalization-factor expanded beta-prior Wald coefficient workflow with Cook's replacement refit.
pub fn fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedBetaPriorWaldReplacementResults, DeseqError> {
    let original = fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
        input.clone(),
        coefficient,
    )?;
    let cooks = beta_prior_normalized_cooks_output(
        input.counts,
        input.normalization_factors,
        &original.fit,
    )?;
    let normalized = normalized_counts_with_factors(input.counts, input.normalization_factors)?;
    let replacement_size_factors = vec![1.0; input.counts.n_samples()];
    let refit_plan = prepare_cooks_replacement_refit(
        input.counts,
        &normalized,
        &replacement_size_factors,
        Some(input.normalization_factors),
        &cooks.cooks,
        input.design.standard_design,
        replacement_options,
    )?;

    let refit = if refit_plan.should_refit {
        Some(
            fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
                ExpandedBetaPriorWaldNormalizedResultsInput {
                    counts: &refit_plan.replacement.replaced_counts,
                    design: input.design,
                    normalization_factors: input.normalization_factors,
                    weights: input.weights,
                    dispersions: input.dispersions,
                    base_mean: &refit_plan.replaced_base_mean,
                    disp_fit: input.disp_fit,
                    gene_names: input.gene_names,
                    options: input.options,
                },
                coefficient,
            )?,
        )
    } else {
        None
    };

    let mut original_results = original.results.clone();
    attach_cooks_to_results(&mut original_results, &cooks.max_cooks)?;
    let mut results = merge_beta_prior_replacement_results(
        &original_results,
        refit.as_ref().map(|value| &value.results),
        &refit_plan,
    )?;
    apply_cooks_cutoff(&mut results, Some(replacement_options.cooks_cutoff))?;

    Ok(ExpandedBetaPriorWaldReplacementResults {
        original: ExpandedBetaPriorWaldResults {
            fit: original.fit,
            results: original_results,
        },
        cooks,
        refit_plan,
        refit,
        results,
    })
}

/// Fit a normalization-factor expanded beta-prior Wald contrast workflow with Cook's replacement refit.
pub fn fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedBetaPriorWaldReplacementResults, DeseqError> {
    let original =
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            input.clone(),
            contrast,
        )?;
    let cooks = beta_prior_normalized_cooks_output(
        input.counts,
        input.normalization_factors,
        &original.fit,
    )?;
    let normalized = normalized_counts_with_factors(input.counts, input.normalization_factors)?;
    let replacement_size_factors = vec![1.0; input.counts.n_samples()];
    let refit_plan = prepare_cooks_replacement_refit(
        input.counts,
        &normalized,
        &replacement_size_factors,
        Some(input.normalization_factors),
        &cooks.cooks,
        input.design.standard_design,
        replacement_options,
    )?;

    let refit = if refit_plan.should_refit {
        Some(
            fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
                ExpandedBetaPriorWaldNormalizedResultsInput {
                    counts: &refit_plan.replacement.replaced_counts,
                    design: input.design,
                    normalization_factors: input.normalization_factors,
                    weights: input.weights,
                    dispersions: input.dispersions,
                    base_mean: &refit_plan.replaced_base_mean,
                    disp_fit: input.disp_fit,
                    gene_names: input.gene_names,
                    options: input.options,
                },
                contrast,
            )?,
        )
    } else {
        None
    };

    let mut original_results = original.results.clone();
    attach_cooks_to_results(&mut original_results, &cooks.max_cooks)?;
    let mut results = merge_beta_prior_replacement_results(
        &original_results,
        refit.as_ref().map(|value| &value.results),
        &refit_plan,
    )?;
    apply_cooks_cutoff(&mut results, Some(replacement_options.cooks_cutoff))?;

    Ok(ExpandedBetaPriorWaldReplacementResults {
        original: ExpandedBetaPriorWaldResults {
            fit: original.fit,
            results: original_results,
        },
        cooks,
        refit_plan,
        refit,
        results,
    })
}

/// Build a one-factor expanded design, fit the beta-prior model, and assemble Wald rows.
pub fn fit_expanded_factor_beta_prior_wald_results(
    input: ExpandedFactorBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedFactorBetaPriorWaldResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_results(
            ExpandedBetaPriorWaldResultsInput {
                counts: input.counts,
                design: design_input,
                size_factors: input.size_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            coefficient,
        )?
    };
    Ok(ExpandedFactorBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build a one-factor expanded design, fit the beta-prior model, and assemble contrast rows.
pub fn fit_expanded_factor_beta_prior_wald_contrast_results(
    input: ExpandedFactorBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedFactorBetaPriorWaldResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_contrast_results(
            ExpandedBetaPriorWaldResultsInput {
                counts: input.counts,
                design: design_input,
                size_factors: input.size_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            contrast,
        )?
    };
    Ok(ExpandedFactorBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build a one-factor expanded design and run coefficient beta-prior Wald replacement refit.
pub fn fit_expanded_factor_beta_prior_wald_results_with_cooks_replacement(
    input: ExpandedFactorBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedFactorBetaPriorWaldReplacementResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
    let replacement = fit_expanded_beta_prior_wald_results_with_cooks_replacement(
        ExpandedBetaPriorWaldResultsInput {
            counts: input.counts,
            design: ExpandedModelBetaPriorDesignInput {
                expanded_design: &design.expanded_design,
                standard_design: &design.standard_design,
                coefficient_groups: &design.coefficient_groups,
            },
            size_factors: input.size_factors,
            weights: input.weights,
            dispersions: input.dispersions,
            base_mean: input.base_mean,
            disp_fit: input.disp_fit,
            gene_names: input.gene_names,
            options: input.options,
        },
        coefficient,
        replacement_options,
    )?;
    Ok(ExpandedFactorBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

/// Build a one-factor expanded design and run contrast beta-prior Wald replacement refit.
pub fn fit_expanded_factor_beta_prior_wald_contrast_results_with_cooks_replacement(
    input: ExpandedFactorBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedFactorBetaPriorWaldReplacementResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
    let replacement = fit_expanded_beta_prior_wald_contrast_results_with_cooks_replacement(
        ExpandedBetaPriorWaldResultsInput {
            counts: input.counts,
            design: ExpandedModelBetaPriorDesignInput {
                expanded_design: &design.expanded_design,
                standard_design: &design.standard_design,
                coefficient_groups: &design.coefficient_groups,
            },
            size_factors: input.size_factors,
            weights: input.weights,
            dispersions: input.dispersions,
            base_mean: input.base_mean,
            disp_fit: input.disp_fit,
            gene_names: input.gene_names,
            options: input.options,
        },
        contrast,
        replacement_options,
    )?;
    Ok(ExpandedFactorBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

/// Build a one-factor expanded design, use normalization factors, and assemble Wald rows.
pub fn fit_expanded_factor_beta_prior_wald_results_with_normalization_factors_and_weights(
    input: ExpandedFactorBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedFactorBetaPriorWaldResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: design_input,
                normalization_factors: input.normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            coefficient,
        )?
    };
    Ok(ExpandedFactorBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build a one-factor expanded design, use normalization factors, and assemble contrast rows.
pub fn fit_expanded_factor_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
    input: ExpandedFactorBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedFactorBetaPriorWaldResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: design_input,
                normalization_factors: input.normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            contrast,
        )?
    };
    Ok(ExpandedFactorBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build a one-factor expanded design, use normalization factors, and run coefficient beta-prior Wald replacement refit.
pub fn fit_expanded_factor_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedFactorBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedFactorBetaPriorWaldReplacementResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
    let replacement =
        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: ExpandedModelBetaPriorDesignInput {
                    expanded_design: &design.expanded_design,
                    standard_design: &design.standard_design,
                    coefficient_groups: &design.coefficient_groups,
                },
                normalization_factors: input.normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            coefficient,
            replacement_options,
        )?;
    Ok(ExpandedFactorBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

/// Build a one-factor expanded design, use normalization factors, and run contrast beta-prior Wald replacement refit.
pub fn fit_expanded_factor_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedFactorBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedFactorBetaPriorWaldReplacementResults, DeseqError> {
    let design = expanded_factor_design(input.factor, input.sample_levels, input.reference)?;
    let replacement =
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: ExpandedModelBetaPriorDesignInput {
                    expanded_design: &design.expanded_design,
                    standard_design: &design.standard_design,
                    coefficient_groups: &design.coefficient_groups,
                },
                normalization_factors: input.normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            contrast,
            replacement_options,
        )?;
    Ok(ExpandedFactorBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

/// Build an additive-factor expanded design, fit the beta-prior model, and assemble Wald rows.
pub fn fit_expanded_additive_beta_prior_wald_results(
    input: ExpandedAdditiveBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_results(
            ExpandedBetaPriorWaldResultsInput {
                counts: input.counts,
                design: design_input,
                size_factors: input.size_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            coefficient,
        )?
    };
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build an additive-factor expanded design, fit the beta-prior model, and assemble contrast rows.
pub fn fit_expanded_additive_beta_prior_wald_contrast_results(
    input: ExpandedAdditiveBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_contrast_results(
            ExpandedBetaPriorWaldResultsInput {
                counts: input.counts,
                design: design_input,
                size_factors: input.size_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            contrast,
        )?
    };
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build an additive-factor expanded design and run coefficient beta-prior Wald replacement refit.
pub fn fit_expanded_additive_beta_prior_wald_results_with_cooks_replacement(
    input: ExpandedAdditiveBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let replacement = fit_expanded_beta_prior_wald_results_with_cooks_replacement(
        ExpandedBetaPriorWaldResultsInput {
            counts: input.counts,
            design: ExpandedModelBetaPriorDesignInput {
                expanded_design: &design.expanded_design,
                standard_design: &design.standard_design,
                coefficient_groups: &design.coefficient_groups,
            },
            size_factors: input.size_factors,
            weights: input.weights,
            dispersions: input.dispersions,
            base_mean: input.base_mean,
            disp_fit: input.disp_fit,
            gene_names: input.gene_names,
            options: input.options,
        },
        coefficient,
        replacement_options,
    )?;
    Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

/// Build an additive-factor expanded design and run contrast beta-prior Wald replacement refit.
pub fn fit_expanded_additive_beta_prior_wald_contrast_results_with_cooks_replacement(
    input: ExpandedAdditiveBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let replacement = fit_expanded_beta_prior_wald_contrast_results_with_cooks_replacement(
        ExpandedBetaPriorWaldResultsInput {
            counts: input.counts,
            design: ExpandedModelBetaPriorDesignInput {
                expanded_design: &design.expanded_design,
                standard_design: &design.standard_design,
                coefficient_groups: &design.coefficient_groups,
            },
            size_factors: input.size_factors,
            weights: input.weights,
            dispersions: input.dispersions,
            base_mean: input.base_mean,
            disp_fit: input.disp_fit,
            gene_names: input.gene_names,
            options: input.options,
        },
        contrast,
        replacement_options,
    )?;
    Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

/// Build an additive-factor expanded design, use normalization factors, and assemble Wald rows.
pub fn fit_expanded_additive_beta_prior_wald_results_with_normalization_factors_and_weights(
    input: ExpandedAdditiveBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: design_input,
                normalization_factors: input.normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            coefficient,
        )?
    };
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build an additive-factor expanded design, use normalization factors, and assemble contrast rows.
pub fn fit_expanded_additive_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
    input: ExpandedAdditiveBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: design_input,
                normalization_factors: input.normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            contrast,
        )?
    };
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Build an additive-factor expanded design, use normalization factors, and run coefficient beta-prior Wald replacement refit.
pub fn fit_expanded_additive_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedAdditiveBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let replacement =
        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: ExpandedModelBetaPriorDesignInput {
                    expanded_design: &design.expanded_design,
                    standard_design: &design.standard_design,
                    coefficient_groups: &design.coefficient_groups,
                },
                normalization_factors: input.normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            coefficient,
            replacement_options,
        )?;
    Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

/// Build an additive-factor expanded design, use normalization factors, and run contrast beta-prior Wald replacement refit.
pub fn fit_expanded_additive_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedAdditiveBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let design = crate::design::expanded_additive_design_with_all_interactions(
        input.factors,
        input.numeric_covariates,
        input.interactions,
        input.factor_numeric_interactions,
        input.numeric_interactions,
    )?;
    let replacement =
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: ExpandedModelBetaPriorDesignInput {
                    expanded_design: &design.expanded_design,
                    standard_design: &design.standard_design,
                    coefficient_groups: &design.coefficient_groups,
                },
                normalization_factors: input.normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            contrast,
            replacement_options,
        )?;
    Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

/// Parse a primitive formula, fit the expanded beta-prior model, and assemble Wald rows.
pub fn fit_expanded_formula_beta_prior_wald_results(
    input: ExpandedFormulaBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    let design = formula_design.design;
    let offset_factors =
        formula_size_factor_offsets(input.counts, input.size_factors, &formula_design.offsets)?;
    if let Some(normalization_factors) = offset_factors.as_ref() {
        let fit_and_results =
            fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
                ExpandedBetaPriorWaldNormalizedResultsInput {
                    counts: input.counts,
                    design: ExpandedModelBetaPriorDesignInput {
                        expanded_design: &design.expanded_design,
                        standard_design: &design.standard_design,
                        coefficient_groups: &design.coefficient_groups,
                    },
                    normalization_factors,
                    weights: input.weights,
                    dispersions: input.dispersions,
                    base_mean: input.base_mean,
                    disp_fit: input.disp_fit,
                    gene_names: input.gene_names,
                    options: input.options,
                },
                coefficient,
            )?;
        return Ok(ExpandedAdditiveBetaPriorWaldResults {
            design,
            fit: fit_and_results.fit,
            results: fit_and_results.results,
        });
    }
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_results(
            ExpandedBetaPriorWaldResultsInput {
                counts: input.counts,
                design: design_input,
                size_factors: input.size_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            coefficient,
        )?
    };
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Parse a primitive formula, fit the expanded beta-prior model, and assemble contrast rows.
pub fn fit_expanded_formula_beta_prior_wald_contrast_results(
    input: ExpandedFormulaBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    let design = formula_design.design;
    let offset_factors =
        formula_size_factor_offsets(input.counts, input.size_factors, &formula_design.offsets)?;
    if let Some(normalization_factors) = offset_factors.as_ref() {
        let fit_and_results =
            fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
                ExpandedBetaPriorWaldNormalizedResultsInput {
                    counts: input.counts,
                    design: ExpandedModelBetaPriorDesignInput {
                        expanded_design: &design.expanded_design,
                        standard_design: &design.standard_design,
                        coefficient_groups: &design.coefficient_groups,
                    },
                    normalization_factors,
                    weights: input.weights,
                    dispersions: input.dispersions,
                    base_mean: input.base_mean,
                    disp_fit: input.disp_fit,
                    gene_names: input.gene_names,
                    options: input.options,
                },
                contrast,
            )?;
        return Ok(ExpandedAdditiveBetaPriorWaldResults {
            design,
            fit: fit_and_results.fit,
            results: fit_and_results.results,
        });
    }
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_contrast_results(
            ExpandedBetaPriorWaldResultsInput {
                counts: input.counts,
                design: design_input,
                size_factors: input.size_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            contrast,
        )?
    };
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Parse a primitive formula and run coefficient beta-prior Wald replacement refit.
pub fn fit_expanded_formula_beta_prior_wald_results_with_cooks_replacement(
    input: ExpandedFormulaBetaPriorWaldResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    let design = formula_design.design;
    let offset_factors =
        formula_size_factor_offsets(input.counts, input.size_factors, &formula_design.offsets)?;
    if let Some(normalization_factors) = offset_factors.as_ref() {
        let replacement =
            fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
                ExpandedBetaPriorWaldNormalizedResultsInput {
                    counts: input.counts,
                    design: ExpandedModelBetaPriorDesignInput {
                        expanded_design: &design.expanded_design,
                        standard_design: &design.standard_design,
                        coefficient_groups: &design.coefficient_groups,
                    },
                    normalization_factors,
                    weights: input.weights,
                    dispersions: input.dispersions,
                    base_mean: input.base_mean,
                    disp_fit: input.disp_fit,
                    gene_names: input.gene_names,
                    options: input.options,
                },
                coefficient,
                replacement_options,
            )?;
        return Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
            design,
            replacement,
        });
    }
    let replacement = fit_expanded_beta_prior_wald_results_with_cooks_replacement(
        ExpandedBetaPriorWaldResultsInput {
            counts: input.counts,
            design: ExpandedModelBetaPriorDesignInput {
                expanded_design: &design.expanded_design,
                standard_design: &design.standard_design,
                coefficient_groups: &design.coefficient_groups,
            },
            size_factors: input.size_factors,
            weights: input.weights,
            dispersions: input.dispersions,
            base_mean: input.base_mean,
            disp_fit: input.disp_fit,
            gene_names: input.gene_names,
            options: input.options,
        },
        coefficient,
        replacement_options,
    )?;
    Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

/// Parse a primitive formula and run contrast beta-prior Wald replacement refit.
pub fn fit_expanded_formula_beta_prior_wald_contrast_results_with_cooks_replacement(
    input: ExpandedFormulaBetaPriorWaldResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    let design = formula_design.design;
    let offset_factors =
        formula_size_factor_offsets(input.counts, input.size_factors, &formula_design.offsets)?;
    if let Some(normalization_factors) = offset_factors.as_ref() {
        let replacement =
            fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
                ExpandedBetaPriorWaldNormalizedResultsInput {
                    counts: input.counts,
                    design: ExpandedModelBetaPriorDesignInput {
                        expanded_design: &design.expanded_design,
                        standard_design: &design.standard_design,
                        coefficient_groups: &design.coefficient_groups,
                    },
                    normalization_factors,
                    weights: input.weights,
                    dispersions: input.dispersions,
                    base_mean: input.base_mean,
                    disp_fit: input.disp_fit,
                    gene_names: input.gene_names,
                    options: input.options,
                },
                contrast,
                replacement_options,
            )?;
        return Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
            design,
            replacement,
        });
    }
    let replacement = fit_expanded_beta_prior_wald_contrast_results_with_cooks_replacement(
        ExpandedBetaPriorWaldResultsInput {
            counts: input.counts,
            design: ExpandedModelBetaPriorDesignInput {
                expanded_design: &design.expanded_design,
                standard_design: &design.standard_design,
                coefficient_groups: &design.coefficient_groups,
            },
            size_factors: input.size_factors,
            weights: input.weights,
            dispersions: input.dispersions,
            base_mean: input.base_mean,
            disp_fit: input.disp_fit,
            gene_names: input.gene_names,
            options: input.options,
        },
        contrast,
        replacement_options,
    )?;
    Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

/// Parse a primitive formula, use normalization factors, and assemble Wald rows.
pub fn fit_expanded_formula_beta_prior_wald_results_with_normalization_factors_and_weights(
    input: ExpandedFormulaBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    let design = formula_design.design;
    let offset_normalization_factors = formula_normalization_factor_offsets(
        input.counts,
        input.normalization_factors,
        &formula_design.offsets,
    )?;
    let normalization_factors = offset_normalization_factors
        .as_ref()
        .unwrap_or(input.normalization_factors);
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: design_input,
                normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            coefficient,
        )?
    };
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Parse a primitive formula, use normalization factors, and assemble contrast rows.
pub fn fit_expanded_formula_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
    input: ExpandedFormulaBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
) -> Result<ExpandedAdditiveBetaPriorWaldResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    let design = formula_design.design;
    let offset_normalization_factors = formula_normalization_factor_offsets(
        input.counts,
        input.normalization_factors,
        &formula_design.offsets,
    )?;
    let normalization_factors = offset_normalization_factors
        .as_ref()
        .unwrap_or(input.normalization_factors);
    let fit_and_results = {
        let design_input = ExpandedModelBetaPriorDesignInput {
            expanded_design: &design.expanded_design,
            standard_design: &design.standard_design,
            coefficient_groups: &design.coefficient_groups,
        };
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: design_input,
                normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            contrast,
        )?
    };
    Ok(ExpandedAdditiveBetaPriorWaldResults {
        design,
        fit: fit_and_results.fit,
        results: fit_and_results.results,
    })
}

/// Parse a primitive formula, use normalization factors, and run coefficient beta-prior Wald replacement refit.
pub fn fit_expanded_formula_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedFormulaBetaPriorWaldNormalizedResultsInput<'_>,
    coefficient: usize,
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    let design = formula_design.design;
    let offset_normalization_factors = formula_normalization_factor_offsets(
        input.counts,
        input.normalization_factors,
        &formula_design.offsets,
    )?;
    let normalization_factors = offset_normalization_factors
        .as_ref()
        .unwrap_or(input.normalization_factors);
    let replacement =
        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: ExpandedModelBetaPriorDesignInput {
                    expanded_design: &design.expanded_design,
                    standard_design: &design.standard_design,
                    coefficient_groups: &design.coefficient_groups,
                },
                normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            coefficient,
            replacement_options,
        )?;
    Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

/// Parse a primitive formula, use normalization factors, and run contrast beta-prior Wald replacement refit.
pub fn fit_expanded_formula_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
    input: ExpandedFormulaBetaPriorWaldNormalizedResultsInput<'_>,
    contrast: &[f64],
    replacement_options: &CooksReplacementOptions,
) -> Result<ExpandedAdditiveBetaPriorWaldReplacementResults, DeseqError> {
    let formula_design = expanded_formula_design_with_offsets(
        input.formula,
        input.factors,
        input.numeric_covariates,
    )?;
    let design = formula_design.design;
    let offset_normalization_factors = formula_normalization_factor_offsets(
        input.counts,
        input.normalization_factors,
        &formula_design.offsets,
    )?;
    let normalization_factors = offset_normalization_factors
        .as_ref()
        .unwrap_or(input.normalization_factors);
    let replacement =
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: input.counts,
                design: ExpandedModelBetaPriorDesignInput {
                    expanded_design: &design.expanded_design,
                    standard_design: &design.standard_design,
                    coefficient_groups: &design.coefficient_groups,
                },
                normalization_factors,
                weights: input.weights,
                dispersions: input.dispersions,
                base_mean: input.base_mean,
                disp_fit: input.disp_fit,
                gene_names: input.gene_names,
                options: input.options,
            },
            contrast,
            replacement_options,
        )?;
    Ok(ExpandedAdditiveBetaPriorWaldReplacementResults {
        design,
        replacement,
    })
}

fn formula_size_factor_offsets(
    counts: &CountMatrix,
    size_factors: &[f64],
    offsets: &[f64],
) -> Result<Option<RowMajorMatrix<f64>>, DeseqError> {
    if !formula_offsets_are_active(offsets) {
        return Ok(None);
    }
    if size_factors.len() != counts.n_samples() {
        return Err(invalid_dimensions(
            "formula offset size factors",
            counts.n_samples(),
            size_factors.len(),
        ));
    }
    let offset_scales = formula_offset_scales(offsets, counts.n_samples())?;
    let mut values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for _ in 0..counts.n_genes() {
        for (sample, size_factor) in size_factors.iter().copied().enumerate() {
            if !size_factor.is_finite() || size_factor <= 0.0 {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("size factor at sample {sample} must be finite and positive"),
                });
            }
            let factor = size_factor * offset_scales[sample];
            if !factor.is_finite() || factor <= 0.0 {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula offset normalization factor at sample {sample} must be finite and positive"
                    ),
                });
            }
            values.push(factor);
        }
    }
    RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values).map(Some)
}

fn formula_normalization_factor_offsets(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    offsets: &[f64],
) -> Result<Option<RowMajorMatrix<f64>>, DeseqError> {
    if !formula_offsets_are_active(offsets) {
        return Ok(None);
    }
    if normalization_factors.n_rows() != counts.n_genes()
        || normalization_factors.n_cols() != counts.n_samples()
    {
        return Err(invalid_dimensions(
            "formula offset normalization factors",
            counts.n_genes() * counts.n_samples(),
            normalization_factors.len(),
        ));
    }
    let offset_scales = formula_offset_scales(offsets, counts.n_samples())?;
    let mut values = Vec::with_capacity(normalization_factors.len());
    for gene in 0..normalization_factors.n_rows() {
        for (sample, value) in normalization_factors.row(gene)?.iter().copied().enumerate() {
            if !value.is_finite() || value <= 0.0 {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "normalization factor at gene {gene}, sample {sample} must be finite and positive"
                    ),
                });
            }
            let factor = value * offset_scales[sample];
            if !factor.is_finite() || factor <= 0.0 {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula offset normalization factor at gene {gene}, sample {sample} must be finite and positive"
                    ),
                });
            }
            values.push(factor);
        }
    }
    RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values).map(Some)
}

fn formula_offsets_are_active(offsets: &[f64]) -> bool {
    offsets.iter().any(|value| *value != 0.0)
}

fn beta_prior_cooks_output(
    counts: &CountMatrix,
    size_factors: &[f64],
    fit: &ExpandedModelBetaPriorGlmFit,
) -> Result<CooksOutput, DeseqError> {
    let normalized = normalized_counts(counts, size_factors)?;
    calculate_cooks_distance(
        counts,
        &normalized,
        &fit.prior_fit.mu,
        &fit.prior_fit.hat_diagonal,
        &fit.prior_fit.model_matrix,
    )
}

fn beta_prior_normalized_cooks_output(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    fit: &ExpandedModelBetaPriorGlmFit,
) -> Result<CooksOutput, DeseqError> {
    let normalized = normalized_counts_with_factors(counts, normalization_factors)?;
    calculate_cooks_distance(
        counts,
        &normalized,
        &fit.prior_fit.mu,
        &fit.prior_fit.hat_diagonal,
        &fit.prior_fit.model_matrix,
    )
}

fn attach_cooks_to_results(
    results: &mut DeseqResults,
    max_cooks: &[Option<f64>],
) -> Result<(), DeseqError> {
    if max_cooks.len() != results.rows.len() {
        return Err(invalid_dimensions(
            "Cook's result rows",
            results.rows.len(),
            max_cooks.len(),
        ));
    }
    for (row, max_cook) in results.rows.iter_mut().zip(max_cooks.iter().copied()) {
        row.max_cooks = max_cook;
        row.cooks_outlier = None;
    }
    Ok(())
}

fn merge_beta_prior_replacement_results(
    original_results: &DeseqResults,
    refit_results: Option<&DeseqResults>,
    refit_plan: &CooksRefitPlan,
) -> Result<DeseqResults, DeseqError> {
    if original_results.rows.len() != refit_plan.replacement.replace.len() {
        return Err(invalid_dimensions(
            "beta-prior replacement result rows",
            refit_plan.replacement.replace.len(),
            original_results.rows.len(),
        ));
    }
    if let Some(refit_results) = refit_results {
        if refit_results.rows.len() != original_results.rows.len() {
            return Err(invalid_dimensions(
                "beta-prior replacement refit result rows",
                original_results.rows.len(),
                refit_results.rows.len(),
            ));
        }
    }
    if refit_plan.replaced_base_mean.len() != original_results.rows.len() {
        return Err(invalid_dimensions(
            "beta-prior replacement baseMean rows",
            original_results.rows.len(),
            refit_plan.replaced_base_mean.len(),
        ));
    }
    if refit_plan.post_refit_max_cooks.len() != original_results.rows.len() {
        return Err(invalid_dimensions(
            "beta-prior replacement maxCooks rows",
            original_results.rows.len(),
            refit_plan.post_refit_max_cooks.len(),
        ));
    }

    let mut merged = original_results.clone();
    for (gene, row) in merged.rows.iter_mut().enumerate() {
        row.base_mean = refit_plan.replaced_base_mean[gene];
        if refit_plan.n_refit > 0 && refit_plan.should_refit {
            row.max_cooks = refit_plan.post_refit_max_cooks[gene];
            row.cooks_outlier = None;
            row.filtered = None;
        }
    }

    if let Some(refit_results) = refit_results {
        for gene in refit_plan.refit_rows.iter().copied() {
            merged.rows[gene] = refit_results.rows[gene].clone();
            merged.rows[gene].base_mean = refit_plan.replaced_base_mean[gene];
            merged.rows[gene].max_cooks = refit_plan.post_refit_max_cooks[gene];
            merged.rows[gene].cooks_outlier = None;
            merged.rows[gene].filtered = None;
        }
    }

    for gene in refit_plan.new_all_zero_rows.iter().copied() {
        clear_replacement_all_zero_result(&mut merged.rows[gene]);
        merged.rows[gene].base_mean = refit_plan.replaced_base_mean[gene];
        if refit_plan.n_refit > 0 && refit_plan.should_refit {
            merged.rows[gene].max_cooks = refit_plan.post_refit_max_cooks[gene];
        }
    }

    merged.independent_filtering = None;
    Ok(merged)
}

fn clear_replacement_all_zero_result(row: &mut DeseqResultRow) {
    row.log2_fold_change = None;
    row.lfc_se = None;
    row.stat = None;
    row.pvalue = None;
    row.padj = None;
    row.dispersion = None;
    row.converged = None;
    row.cooks_outlier = None;
    row.filtered = None;
}

fn formula_offset_scales(offsets: &[f64], n_samples: usize) -> Result<Vec<f64>, DeseqError> {
    if offsets.len() != n_samples {
        return Err(invalid_dimensions(
            "formula offsets",
            n_samples,
            offsets.len(),
        ));
    }
    offsets
        .iter()
        .copied()
        .enumerate()
        .map(|(sample, offset)| {
            let scale = offset.exp();
            if !scale.is_finite() || scale <= 0.0 {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula offset scale at sample {sample} must be finite and positive"
                    ),
                });
            }
            Ok(scale)
        })
        .collect()
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
        if let Some(value) = value {
            if !value.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: context.to_string(),
                    index: Some(idx),
                    value,
                });
            }
        }
    }
    Ok(())
}

fn validate_optional_positive_finite(
    values: &[Option<f64>],
    context: &str,
) -> Result<(), DeseqError> {
    for (idx, value) in values.iter().copied().enumerate() {
        if let Some(value) = value {
            if !value.is_finite() || value <= 0.0 {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("{context} at index {idx} must be positive and finite"),
                });
            }
        }
    }
    Ok(())
}

fn validate_optional_probability(values: &[Option<f64>], context: &str) -> Result<(), DeseqError> {
    for (idx, value) in values.iter().copied().enumerate() {
        if let Some(value) = value {
            if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("{context} at index {idx} must be finite and within [0, 1]"),
                });
            }
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
