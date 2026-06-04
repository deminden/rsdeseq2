use std::collections::HashSet;
use std::path::PathBuf;

use clap::{ArgAction, Parser, Subcommand, ValueEnum};

use crate::contrasts::{
    resolve_coefficient_index, resolve_contrast, ContrastSpec, ResultsContrast,
};
use crate::cooks::{CooksRefitPlan, CooksReplacementOptions};
use crate::core::{
    CooksReplacementLrtOutput, CooksReplacementTestOutput, CooksReplacementWaldOutput,
    DeseqBuilder, DeseqFit,
};
use crate::design::{
    expanded_additive_design, expanded_additive_factor_design, expanded_factor_design,
    ExpandedFactorSpec, ExpandedNumericSpec,
};
use crate::errors::DeseqError;
use crate::glm::{BetaPriorRefitOptions, ExpandedModelBetaPriorDesignInput, WaldAlternative};
use crate::io::{
    align_design_matrix_to_samples, align_gene_numeric_values_to_genes,
    align_labeled_assay_matrix_to_counts, align_sample_levels_to_samples,
    align_sample_numeric_values_to_samples, read_count_matrix_tsv, read_labeled_design_matrix_tsv,
    read_labeled_gene_numeric_tsv, read_labeled_geometric_means_tsv,
    read_labeled_normalization_factors_tsv, read_labeled_observation_weights_tsv,
    read_labeled_sample_numeric_tsv, read_labeled_size_factors_tsv,
    read_labeled_wald_t_degrees_of_freedom_tsv, read_sample_levels_tsv, write_base_mean_tsv,
    write_cooks_candidate_replacement_counts_tsv, write_cooks_distance_matrix_tsv,
    write_cooks_outlier_cells_tsv, write_cooks_replaced_counts_tsv,
    write_cooks_replacement_metadata_tsv, write_cooks_replacement_row_metadata_tsv,
    write_deseq_mcols_diagnostics_tsv, write_deseq_result_column_metadata_tsv,
    write_deseq_result_table_metadata_tsv, write_deseq_results_tsv,
    write_independent_filter_lowess_tsv, write_independent_filter_metadata_tsv,
    write_independent_filter_num_rej_tsv, write_normalized_counts_tsv,
    write_optional_numeric_matrix_tsv, write_size_factors_tsv,
};
use crate::matrix::RowMajorMatrix;
use crate::normalization::{
    base_mean, base_mean_with_weights, estimate_size_factors_with_options, normalized_counts,
    normalized_counts_with_factors,
};
use crate::options::{CooksCutoff, FitType, SizeFactorMethod, TestType};
use crate::results::{
    fit_expanded_additive_beta_prior_wald_contrast_results,
    fit_expanded_additive_beta_prior_wald_contrast_results_with_cooks_replacement,
    fit_expanded_additive_beta_prior_wald_contrast_results_with_normalization_factors_and_weights,
    fit_expanded_additive_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement,
    fit_expanded_additive_beta_prior_wald_results,
    fit_expanded_additive_beta_prior_wald_results_with_cooks_replacement,
    fit_expanded_additive_beta_prior_wald_results_with_normalization_factors_and_weights,
    fit_expanded_additive_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement,
    fit_expanded_beta_prior_wald_contrast_results,
    fit_expanded_beta_prior_wald_contrast_results_with_cooks_replacement,
    fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights,
    fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement,
    fit_expanded_beta_prior_wald_results,
    fit_expanded_beta_prior_wald_results_with_cooks_replacement,
    fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights,
    fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement,
    fit_expanded_factor_beta_prior_wald_contrast_results,
    fit_expanded_factor_beta_prior_wald_contrast_results_with_cooks_replacement,
    fit_expanded_factor_beta_prior_wald_contrast_results_with_normalization_factors_and_weights,
    fit_expanded_factor_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement,
    fit_expanded_factor_beta_prior_wald_results,
    fit_expanded_factor_beta_prior_wald_results_with_cooks_replacement,
    fit_expanded_factor_beta_prior_wald_results_with_normalization_factors_and_weights,
    fit_expanded_factor_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement,
    resolve_cooks_cutoff, DeseqResults, ExpandedAdditiveBetaPriorWaldNormalizedResultsInput,
    ExpandedAdditiveBetaPriorWaldReplacementResults, ExpandedAdditiveBetaPriorWaldResults,
    ExpandedAdditiveBetaPriorWaldResultsInput, ExpandedBetaPriorWaldNormalizedResultsInput,
    ExpandedBetaPriorWaldReplacementResults, ExpandedBetaPriorWaldResults,
    ExpandedBetaPriorWaldResultsInput, ExpandedFactorBetaPriorWaldNormalizedResultsInput,
    ExpandedFactorBetaPriorWaldReplacementResults, ExpandedFactorBetaPriorWaldResults,
    ExpandedFactorBetaPriorWaldResultsInput,
};

/// Command-line arguments for the minimal `rsdeseq2` CLI.
#[derive(Debug, Parser)]
#[command(name = "rsdeseq2")]
#[command(about = "Early DESeq2-compatible Rust workflow stages")]
pub struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
#[allow(clippy::large_enum_variant)]
enum Commands {
    /// Estimate sample size factors.
    SizeFactors {
        /// Tab-delimited count matrix with gene IDs in the first column.
        #[arg(long)]
        counts: PathBuf,
        /// Size-factor method.
        #[arg(long, default_value = "ratio")]
        method: SizeFactorMethodArg,
        /// Optional gene x geometric-mean TSV used for frozen size-factor estimation.
        #[arg(long)]
        geometric_means: Option<PathBuf>,
        /// Comma-delimited zero-based row indices used to estimate size factors.
        #[arg(long, value_delimiter = ',')]
        control_genes: Option<Vec<usize>>,
        /// Output TSV path.
        #[arg(long)]
        output: PathBuf,
    },
    /// Estimate size factors, normalized counts, and base means.
    BaseMean {
        /// Tab-delimited count matrix with gene IDs in the first column.
        #[arg(long)]
        counts: PathBuf,
        /// Optional gene x sample normalization-factor TSV.
        #[arg(long)]
        normalization_factors: Option<PathBuf>,
        /// Optional sample-level size-factor TSV.
        #[arg(long)]
        size_factors: Option<PathBuf>,
        /// Optional gene x sample observation-weight TSV.
        #[arg(long)]
        observation_weights: Option<PathBuf>,
        /// Size-factor method.
        #[arg(long, default_value = "ratio")]
        method: SizeFactorMethodArg,
        /// Optional gene x geometric-mean TSV used for frozen size-factor estimation.
        #[arg(long)]
        geometric_means: Option<PathBuf>,
        /// Comma-delimited zero-based row indices used to estimate size factors.
        #[arg(long, value_delimiter = ',')]
        control_genes: Option<Vec<usize>>,
        /// Output TSV path.
        #[arg(long)]
        output: PathBuf,
    },
    /// Write DESeq2-style normalized counts.
    NormalizedCounts {
        /// Tab-delimited count matrix with gene IDs in the first column.
        #[arg(long)]
        counts: PathBuf,
        /// Optional gene x sample normalization-factor TSV.
        #[arg(long)]
        normalization_factors: Option<PathBuf>,
        /// Optional sample-level size-factor TSV.
        #[arg(long)]
        size_factors: Option<PathBuf>,
        /// Size-factor method.
        #[arg(long, default_value = "ratio")]
        method: SizeFactorMethodArg,
        /// Optional gene x geometric-mean TSV used for frozen size-factor estimation.
        #[arg(long)]
        geometric_means: Option<PathBuf>,
        /// Comma-delimited zero-based row indices used to estimate size factors.
        #[arg(long, value_delimiter = ',')]
        control_genes: Option<Vec<usize>>,
        /// Output TSV path.
        #[arg(long)]
        output: PathBuf,
    },
    /// Write a variance-stabilized count matrix.
    Vst {
        /// Tab-delimited count matrix with gene IDs in the first column.
        #[arg(long)]
        counts: PathBuf,
        /// Numeric design matrix TSV. Required with --blind false.
        #[arg(long)]
        design: Option<PathBuf>,
        /// Ignore the design and use an intercept-only trend fit.
        #[arg(long, default_value_t = true, action = ArgAction::Set)]
        blind: bool,
        /// Optional gene x sample normalization-factor TSV.
        #[arg(long)]
        normalization_factors: Option<PathBuf>,
        /// Optional sample-level size-factor TSV.
        #[arg(long)]
        size_factors: Option<PathBuf>,
        /// Optional gene x sample observation-weight TSV.
        #[arg(long)]
        observation_weights: Option<PathBuf>,
        /// Size-factor method.
        #[arg(long, default_value = "ratio")]
        method: SizeFactorMethodArg,
        /// Optional gene x geometric-mean TSV used for frozen size-factor estimation.
        #[arg(long)]
        geometric_means: Option<PathBuf>,
        /// Comma-delimited zero-based row indices used to estimate size factors.
        #[arg(long, value_delimiter = ',')]
        control_genes: Option<Vec<usize>>,
        /// Dispersion trend fit type.
        #[arg(long, default_value = "parametric")]
        fit_type: FitTypeArg,
        /// Fast-VST subset size.
        #[arg(long, default_value_t = 1000)]
        nsub: usize,
        /// Output TSV path.
        #[arg(long)]
        output: PathBuf,
    },
    /// Write a regularized-log count matrix.
    Rlog {
        /// Tab-delimited count matrix with gene IDs in the first column.
        #[arg(long)]
        counts: PathBuf,
        /// Numeric design matrix TSV. Required with --blind false.
        #[arg(long)]
        design: Option<PathBuf>,
        /// Ignore the design and use an intercept-only dispersion workflow.
        #[arg(long, default_value_t = true, action = ArgAction::Set)]
        blind: bool,
        /// Optional gene x sample normalization-factor TSV.
        #[arg(long)]
        normalization_factors: Option<PathBuf>,
        /// Optional sample-level size-factor TSV.
        #[arg(long)]
        size_factors: Option<PathBuf>,
        /// Optional gene x sample observation-weight TSV for dispersion estimation.
        #[arg(long)]
        observation_weights: Option<PathBuf>,
        /// Size-factor method.
        #[arg(long, default_value = "ratio")]
        method: SizeFactorMethodArg,
        /// Optional gene x geometric-mean TSV used for frozen size-factor estimation.
        #[arg(long)]
        geometric_means: Option<PathBuf>,
        /// Comma-delimited zero-based row indices used to estimate size factors.
        #[arg(long, value_delimiter = ',')]
        control_genes: Option<Vec<usize>>,
        /// Dispersion trend fit type.
        #[arg(long, default_value = "parametric")]
        fit_type: FitTypeArg,
        /// Optional gene x frozen rlog intercept TSV for a frozen-intercept transform.
        #[arg(long)]
        frozen_intercept: Option<PathBuf>,
        /// Sample-effect prior variance to use with --frozen-intercept.
        #[arg(long)]
        rlog_prior_variance: Option<f64>,
        /// Output TSV path.
        #[arg(long)]
        output: PathBuf,
    },
    /// Run the implemented GLM-mu Wald workflow and write DESeq2-shaped results.
    Wald {
        /// Tab-delimited count matrix with gene IDs in the first column.
        #[arg(long)]
        counts: PathBuf,
        /// Numeric design matrix TSV with sample IDs in the first column.
        #[arg(long)]
        design: PathBuf,
        /// Optional gene x sample normalization-factor TSV.
        #[arg(long)]
        normalization_factors: Option<PathBuf>,
        /// Optional sample-level size-factor TSV.
        #[arg(long)]
        size_factors: Option<PathBuf>,
        /// Optional gene x sample observation-weight TSV.
        #[arg(long)]
        observation_weights: Option<PathBuf>,
        /// Size-factor method.
        #[arg(long, default_value = "ratio")]
        method: SizeFactorMethodArg,
        /// Optional gene x geometric-mean TSV used for frozen size-factor estimation.
        #[arg(long)]
        geometric_means: Option<PathBuf>,
        /// Comma-delimited zero-based row indices used to estimate size factors.
        #[arg(long, value_delimiter = ',')]
        control_genes: Option<Vec<usize>>,
        /// Dispersion trend fit type.
        #[arg(long, default_value = "parametric")]
        fit_type: FitTypeArg,
        /// Zero-based coefficient index to report. Defaults to the last column.
        #[arg(long)]
        coefficient: Option<usize>,
        /// Design coefficient name to report. Defaults to the last column when no coefficient or contrast is supplied.
        #[arg(long)]
        coefficient_name: Option<String>,
        /// Numeric Wald contrast vector, comma-delimited in design-column order.
        #[arg(long, value_delimiter = ',')]
        contrast: Option<Vec<f64>>,
        /// Expanded design matrix TSV for the supplied-dispersion beta-prior Wald workflow.
        #[arg(long)]
        beta_prior_expanded_design: Option<PathBuf>,
        /// Collapse groups for beta-prior expanded columns, e.g. `0|1,2|3`.
        #[arg(long)]
        beta_prior_coefficient_groups: Option<String>,
        /// Gene x final-dispersion TSV for the beta-prior fixed-dispersion workflow.
        #[arg(long)]
        beta_prior_dispersions: Option<PathBuf>,
        /// Gene x baseMean TSV for beta-prior result rows and prior weights.
        #[arg(long)]
        beta_prior_base_mean: Option<PathBuf>,
        /// Gene x fitted-dispersion TSV for beta-prior prior weights.
        #[arg(long)]
        beta_prior_disp_fit: Option<PathBuf>,
        /// Factor name for internally building a one-factor expanded beta-prior design.
        #[arg(long)]
        beta_prior_factor: Option<String>,
        /// Reference level for internally building a one-factor expanded beta-prior design.
        #[arg(long)]
        beta_prior_reference: Option<String>,
        /// Optional sample x level TSV for internally building a one-factor expanded beta-prior design.
        #[arg(long)]
        beta_prior_sample_levels: Option<PathBuf>,
        /// Comma-delimited factor names for internally building an additive expanded beta-prior design.
        #[arg(long, value_delimiter = ',')]
        beta_prior_additive_factors: Option<Vec<String>>,
        /// Comma-delimited reference levels for additive beta-prior factors.
        #[arg(long, value_delimiter = ',')]
        beta_prior_additive_references: Option<Vec<String>>,
        /// Comma-delimited sample x level TSV paths for additive beta-prior factors.
        #[arg(long, value_delimiter = ',')]
        beta_prior_additive_sample_levels: Option<Vec<PathBuf>>,
        /// Comma-delimited numeric covariate names for additive beta-prior designs.
        #[arg(long, value_delimiter = ',')]
        beta_prior_additive_numeric: Option<Vec<String>>,
        /// Comma-delimited sample x value TSV paths for additive beta-prior numeric covariates.
        #[arg(long, value_delimiter = ',')]
        beta_prior_additive_numeric_values: Option<Vec<PathBuf>>,
        /// Design coefficient name to report as a Wald contrast.
        #[arg(long)]
        contrast_name: Option<String>,
        /// Comma-delimited positive coefficient names for a list contrast.
        #[arg(long, value_delimiter = ',')]
        contrast_positive: Option<Vec<String>>,
        /// Comma-delimited negative coefficient names for a list contrast.
        #[arg(long, value_delimiter = ',')]
        contrast_negative: Option<Vec<String>>,
        /// Weight applied to positive list-contrast coefficients.
        #[arg(long, default_value_t = 1.0)]
        contrast_positive_weight: f64,
        /// Weight applied to negative list-contrast coefficients.
        #[arg(long, default_value_t = -1.0, allow_hyphen_values = true)]
        contrast_negative_weight: f64,
        /// Factor or variable name for a coefficient-name factor-level contrast.
        #[arg(long)]
        contrast_factor: Option<String>,
        /// Numerator level for a coefficient-name factor-level contrast.
        #[arg(long)]
        contrast_numerator: Option<String>,
        /// Denominator level for a coefficient-name factor-level contrast.
        #[arg(long)]
        contrast_denominator: Option<String>,
        /// Optional reference level for a coefficient-name factor-level contrast.
        #[arg(long)]
        contrast_reference: Option<String>,
        /// Sample x level TSV required for DESeq2-style factor-level contrast handling.
        #[arg(long)]
        contrast_sample_levels: Option<PathBuf>,
        /// Non-negative log2 fold-change threshold for Wald p-values.
        #[arg(long, default_value_t = 0.0)]
        lfc_threshold: f64,
        /// Alternative hypothesis for thresholded Wald p-values.
        #[arg(long, default_value = "greater-abs")]
        alternative: WaldAlternativeArg,
        /// Use Student t Wald p-values with residual degrees of freedom.
        #[arg(long, action = ArgAction::SetTrue)]
        use_t: bool,
        /// Use Student t Wald p-values with one scalar degrees-of-freedom value.
        #[arg(long)]
        t_degrees_of_freedom: Option<f64>,
        /// Gene x degrees-of-freedom TSV for Student t Wald p-values.
        #[arg(long)]
        t_degrees_of_freedom_file: Option<PathBuf>,
        /// Disable Cook's distance p-value filtering and replacement/refit.
        #[arg(long, action = ArgAction::SetTrue)]
        disable_cooks_cutoff: bool,
        /// Explicit Cook's distance cutoff. Conflicts with --disable-cooks-cutoff.
        #[arg(long)]
        cooks_cutoff: Option<f64>,
        /// Disable independent filtering and use regular BH adjustment.
        #[arg(long, action = ArgAction::SetTrue)]
        disable_independent_filtering: bool,
        /// Alpha used for independent-filtering threshold selection.
        #[arg(long)]
        independent_filtering_alpha: Option<f64>,
        /// Comma-delimited theta grid for independent filtering.
        #[arg(long, value_delimiter = ',')]
        independent_filtering_theta: Option<Vec<f64>>,
        /// Optional result column metadata TSV output.
        #[arg(long)]
        result_column_metadata_output: Option<PathBuf>,
        /// Optional result table metadata TSV output.
        #[arg(long)]
        result_table_metadata_output: Option<PathBuf>,
        /// Optional independent-filtering scalar metadata TSV output.
        #[arg(long)]
        independent_filter_metadata_output: Option<PathBuf>,
        /// Optional independent-filtering rejection-count curve TSV output.
        #[arg(long)]
        independent_filter_num_rej_output: Option<PathBuf>,
        /// Optional independent-filtering lowess curve TSV output.
        #[arg(long)]
        independent_filter_lowess_output: Option<PathBuf>,
        /// Optional DESeq2-shaped fit diagnostics TSV for the original fit.
        #[arg(long)]
        fit_diagnostics_output: Option<PathBuf>,
        /// Optional DESeq2-shaped fit diagnostics TSV for the replacement refit.
        #[arg(long)]
        refit_diagnostics_output: Option<PathBuf>,
        /// Optional GLM beta matrix TSV for the original fit.
        #[arg(long)]
        fit_beta_output: Option<PathBuf>,
        /// Optional GLM beta standard-error matrix TSV for the original fit.
        #[arg(long)]
        fit_beta_se_output: Option<PathBuf>,
        /// Optional fallback optimizer start beta matrix TSV for the original fit.
        #[arg(long)]
        fit_beta_optim_start_output: Option<PathBuf>,
        /// Optional GLM beta matrix TSV for the replacement refit.
        #[arg(long)]
        refit_beta_output: Option<PathBuf>,
        /// Optional GLM beta standard-error matrix TSV for the replacement refit.
        #[arg(long)]
        refit_beta_se_output: Option<PathBuf>,
        /// Optional fallback optimizer start beta matrix TSV for the replacement refit.
        #[arg(long)]
        refit_beta_optim_start_output: Option<PathBuf>,
        /// Optional Cook's distance matrix TSV output.
        #[arg(long)]
        cooks_distance_output: Option<PathBuf>,
        /// Optional Cook's replacement scalar metadata TSV output.
        #[arg(long)]
        cooks_replacement_metadata_output: Option<PathBuf>,
        /// Optional Cook's replacement row metadata TSV output.
        #[arg(long)]
        cooks_replacement_row_metadata_output: Option<PathBuf>,
        /// Optional Cook's replaced-count assay TSV output.
        #[arg(long)]
        cooks_replaced_counts_output: Option<PathBuf>,
        /// Optional Cook's candidate replacement-count assay TSV output.
        #[arg(long)]
        cooks_candidate_replacement_counts_output: Option<PathBuf>,
        /// Optional Cook's outlier-cell logical assay TSV output.
        #[arg(long)]
        cooks_outlier_cells_output: Option<PathBuf>,
        /// Output TSV path.
        #[arg(long)]
        output: PathBuf,
    },
    /// Run the implemented GLM-mu LRT workflow and write DESeq2-shaped results.
    Lrt {
        /// Tab-delimited count matrix with gene IDs in the first column.
        #[arg(long)]
        counts: PathBuf,
        /// Full numeric design matrix TSV with sample IDs in the first column.
        #[arg(long)]
        design: PathBuf,
        /// Reduced numeric design matrix TSV with sample IDs in the first column.
        #[arg(long)]
        reduced_design: PathBuf,
        /// Optional gene x sample normalization-factor TSV.
        #[arg(long)]
        normalization_factors: Option<PathBuf>,
        /// Optional sample-level size-factor TSV.
        #[arg(long)]
        size_factors: Option<PathBuf>,
        /// Optional gene x sample observation-weight TSV.
        #[arg(long)]
        observation_weights: Option<PathBuf>,
        /// Size-factor method.
        #[arg(long, default_value = "ratio")]
        method: SizeFactorMethodArg,
        /// Optional gene x geometric-mean TSV used for frozen size-factor estimation.
        #[arg(long)]
        geometric_means: Option<PathBuf>,
        /// Comma-delimited zero-based row indices used to estimate size factors.
        #[arg(long, value_delimiter = ',')]
        control_genes: Option<Vec<usize>>,
        /// Dispersion trend fit type.
        #[arg(long, default_value = "parametric")]
        fit_type: FitTypeArg,
        /// Zero-based full-design coefficient index to report. Defaults to the last column.
        #[arg(long)]
        coefficient: Option<usize>,
        /// Full-design coefficient name to report.
        #[arg(long)]
        coefficient_name: Option<String>,
        /// Numeric full-model contrast vector to report, comma-delimited in design-column order.
        #[arg(long, value_delimiter = ',')]
        contrast: Option<Vec<f64>>,
        /// Full-design coefficient name to report as an LRT effect-size contrast.
        #[arg(long)]
        contrast_name: Option<String>,
        /// Comma-delimited positive full-design coefficient names for a list contrast.
        #[arg(long, value_delimiter = ',')]
        contrast_positive: Option<Vec<String>>,
        /// Comma-delimited negative full-design coefficient names for a list contrast.
        #[arg(long, value_delimiter = ',')]
        contrast_negative: Option<Vec<String>>,
        /// Weight applied to positive list-contrast coefficients.
        #[arg(long, default_value_t = 1.0)]
        contrast_positive_weight: f64,
        /// Weight applied to negative list-contrast coefficients.
        #[arg(long, default_value_t = -1.0)]
        contrast_negative_weight: f64,
        /// Factor or variable name for a coefficient-name factor-level contrast.
        #[arg(long)]
        contrast_factor: Option<String>,
        /// Numerator level for a coefficient-name factor-level contrast.
        #[arg(long)]
        contrast_numerator: Option<String>,
        /// Denominator level for a coefficient-name factor-level contrast.
        #[arg(long)]
        contrast_denominator: Option<String>,
        /// Optional reference level for a coefficient-name factor-level contrast.
        #[arg(long)]
        contrast_reference: Option<String>,
        /// Sample x level TSV required for DESeq2-style factor-level contrast handling.
        #[arg(long)]
        contrast_sample_levels: Option<PathBuf>,
        /// Disable Cook's distance p-value filtering and replacement/refit.
        #[arg(long, action = ArgAction::SetTrue)]
        disable_cooks_cutoff: bool,
        /// Explicit Cook's distance cutoff. Conflicts with --disable-cooks-cutoff.
        #[arg(long)]
        cooks_cutoff: Option<f64>,
        /// Disable independent filtering and use regular BH adjustment.
        #[arg(long, action = ArgAction::SetTrue)]
        disable_independent_filtering: bool,
        /// Alpha used for independent-filtering threshold selection.
        #[arg(long)]
        independent_filtering_alpha: Option<f64>,
        /// Comma-delimited theta grid for independent filtering.
        #[arg(long, value_delimiter = ',')]
        independent_filtering_theta: Option<Vec<f64>>,
        /// Optional result column metadata TSV output.
        #[arg(long)]
        result_column_metadata_output: Option<PathBuf>,
        /// Optional result table metadata TSV output.
        #[arg(long)]
        result_table_metadata_output: Option<PathBuf>,
        /// Optional independent-filtering scalar metadata TSV output.
        #[arg(long)]
        independent_filter_metadata_output: Option<PathBuf>,
        /// Optional independent-filtering rejection-count curve TSV output.
        #[arg(long)]
        independent_filter_num_rej_output: Option<PathBuf>,
        /// Optional independent-filtering lowess curve TSV output.
        #[arg(long)]
        independent_filter_lowess_output: Option<PathBuf>,
        /// Optional DESeq2-shaped fit diagnostics TSV for the original fit.
        #[arg(long)]
        fit_diagnostics_output: Option<PathBuf>,
        /// Optional DESeq2-shaped fit diagnostics TSV for the replacement refit.
        #[arg(long)]
        refit_diagnostics_output: Option<PathBuf>,
        /// Optional GLM beta matrix TSV for the original fit.
        #[arg(long)]
        fit_beta_output: Option<PathBuf>,
        /// Optional GLM beta standard-error matrix TSV for the original fit.
        #[arg(long)]
        fit_beta_se_output: Option<PathBuf>,
        /// Optional fallback optimizer start beta matrix TSV for the original fit.
        #[arg(long)]
        fit_beta_optim_start_output: Option<PathBuf>,
        /// Optional GLM beta matrix TSV for the replacement refit.
        #[arg(long)]
        refit_beta_output: Option<PathBuf>,
        /// Optional GLM beta standard-error matrix TSV for the replacement refit.
        #[arg(long)]
        refit_beta_se_output: Option<PathBuf>,
        /// Optional fallback optimizer start beta matrix TSV for the replacement refit.
        #[arg(long)]
        refit_beta_optim_start_output: Option<PathBuf>,
        /// Optional Cook's distance matrix TSV output.
        #[arg(long)]
        cooks_distance_output: Option<PathBuf>,
        /// Optional Cook's replacement scalar metadata TSV output.
        #[arg(long)]
        cooks_replacement_metadata_output: Option<PathBuf>,
        /// Optional Cook's replacement row metadata TSV output.
        #[arg(long)]
        cooks_replacement_row_metadata_output: Option<PathBuf>,
        /// Optional Cook's replaced-count assay TSV output.
        #[arg(long)]
        cooks_replaced_counts_output: Option<PathBuf>,
        /// Optional Cook's candidate replacement-count assay TSV output.
        #[arg(long)]
        cooks_candidate_replacement_counts_output: Option<PathBuf>,
        /// Optional Cook's outlier-cell logical assay TSV output.
        #[arg(long)]
        cooks_outlier_cells_output: Option<PathBuf>,
        /// Output TSV path.
        #[arg(long)]
        output: PathBuf,
    },
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum SizeFactorMethodArg {
    Ratio,
    Poscounts,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum FitTypeArg {
    Parametric,
    Local,
    Mean,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum WaldAlternativeArg {
    GreaterAbs,
    GreaterAbsUpshot,
    GreaterAbs2014,
    LessAbs,
    Greater,
    Less,
}

struct CliAnalysisOutput {
    results: DeseqResults,
    fit: Option<DeseqFit>,
    refit: Option<DeseqFit>,
    cooks: Option<RowMajorMatrix<f64>>,
    refit_plan: Option<CooksRefitPlan>,
}

struct CliCooksOutputPaths {
    cooks_distance: Option<PathBuf>,
    replacement_metadata: Option<PathBuf>,
    replacement_row_metadata: Option<PathBuf>,
    replaced_counts: Option<PathBuf>,
    candidate_replacement_counts: Option<PathBuf>,
    outlier_cells: Option<PathBuf>,
}

struct CliResultSidecarPaths {
    column_metadata: Option<PathBuf>,
    table_metadata: Option<PathBuf>,
    independent_filter_metadata: Option<PathBuf>,
    independent_filter_num_rej: Option<PathBuf>,
    independent_filter_lowess: Option<PathBuf>,
    fit_diagnostics: Option<PathBuf>,
    refit_diagnostics: Option<PathBuf>,
    fit_beta: Option<PathBuf>,
    fit_beta_se: Option<PathBuf>,
    fit_beta_optim_start: Option<PathBuf>,
    refit_beta: Option<PathBuf>,
    refit_beta_se: Option<PathBuf>,
    refit_beta_optim_start: Option<PathBuf>,
}

impl From<SizeFactorMethodArg> for SizeFactorMethod {
    fn from(value: SizeFactorMethodArg) -> Self {
        match value {
            SizeFactorMethodArg::Ratio => Self::Ratio,
            SizeFactorMethodArg::Poscounts => Self::PosCounts,
        }
    }
}

impl From<FitTypeArg> for FitType {
    fn from(value: FitTypeArg) -> Self {
        match value {
            FitTypeArg::Parametric => Self::Parametric,
            FitTypeArg::Local => Self::Local,
            FitTypeArg::Mean => Self::Mean,
        }
    }
}

impl From<WaldAlternativeArg> for WaldAlternative {
    fn from(value: WaldAlternativeArg) -> Self {
        match value {
            WaldAlternativeArg::GreaterAbs => Self::GreaterAbs,
            WaldAlternativeArg::GreaterAbsUpshot => Self::GreaterAbsUpshot,
            WaldAlternativeArg::GreaterAbs2014 => Self::GreaterAbs2014,
            WaldAlternativeArg::LessAbs => Self::LessAbs,
            WaldAlternativeArg::Greater => Self::Greater,
            WaldAlternativeArg::Less => Self::Less,
        }
    }
}
