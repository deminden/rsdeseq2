use std::path::PathBuf;

use clap::{ArgAction, Parser, Subcommand, ValueEnum};

use crate::contrasts::{ContrastSpec, FactorLevelContrast};
use crate::cooks::{CooksRefitPlan, CooksReplacementOptions};
use crate::core::{CooksReplacementLrtOutput, CooksReplacementWaldOutput, DeseqBuilder, DeseqFit};
use crate::errors::DeseqError;
use crate::glm::WaldAlternative;
use crate::io::{
    align_design_matrix_to_samples, align_gene_numeric_values_to_genes,
    align_labeled_assay_matrix_to_counts, align_sample_levels_to_samples,
    align_sample_numeric_values_to_samples, read_count_matrix_tsv, read_labeled_design_matrix_tsv,
    read_labeled_gene_numeric_tsv, read_labeled_geometric_means_tsv,
    read_labeled_normalization_factors_tsv, read_labeled_observation_weights_tsv,
    read_labeled_size_factors_tsv, read_labeled_wald_t_degrees_of_freedom_tsv,
    read_sample_levels_tsv, write_base_mean_tsv, write_cooks_candidate_replacement_counts_tsv,
    write_cooks_distance_matrix_tsv, write_cooks_outlier_cells_tsv,
    write_cooks_replaced_counts_tsv, write_cooks_replacement_metadata_tsv,
    write_cooks_replacement_row_metadata_tsv, write_deseq_result_column_metadata_tsv,
    write_deseq_result_table_metadata_tsv, write_deseq_results_tsv,
    write_independent_filter_lowess_tsv, write_independent_filter_metadata_tsv,
    write_independent_filter_num_rej_tsv, write_normalized_counts_tsv, write_size_factors_tsv,
};
use crate::matrix::RowMajorMatrix;
use crate::normalization::{
    base_mean, base_mean_with_weights, estimate_size_factors_with_options, normalized_counts,
    normalized_counts_with_factors,
};
use crate::options::{CooksCutoff, FitType, SizeFactorMethod};
use crate::results::{resolve_cooks_cutoff, DeseqResults};

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
        /// Optional sample x level TSV for DESeq2-style factor-level all-zero contrast handling.
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
        /// Optional sample x level TSV for DESeq2-style factor-level all-zero contrast handling.
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

/// Parse process arguments and run the CLI.
pub fn run_cli() -> Result<(), DeseqError> {
    run(Cli::parse())
}

fn run(cli: Cli) -> Result<(), DeseqError> {
    match cli.command {
        Commands::SizeFactors {
            counts,
            method,
            geometric_means,
            control_genes,
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let geometric_means = read_cli_geometric_means(geometric_means, &counts)?;
            let size_factors = estimate_size_factors_with_options(
                &counts,
                method.into(),
                geometric_means.as_deref(),
                control_genes.as_deref(),
            )?;
            write_size_factors_tsv(output, counts.sample_names(), &size_factors)
        }
        Commands::BaseMean {
            counts,
            normalization_factors,
            size_factors,
            observation_weights,
            method,
            geometric_means,
            control_genes,
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let normalized = cli_normalized_counts(
                &counts,
                normalization_factors,
                size_factors,
                method,
                geometric_means,
                control_genes,
            )?;
            let base_mean = if let Some(path) = observation_weights {
                let weights = read_cli_observation_weights(path, &counts)?;
                base_mean_with_weights(&normalized, &weights)?
            } else {
                base_mean(&normalized)?
            };
            write_base_mean_tsv(output, counts.gene_names(), &base_mean)
        }
        Commands::NormalizedCounts {
            counts,
            normalization_factors,
            size_factors,
            method,
            geometric_means,
            control_genes,
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let normalized = cli_normalized_counts(
                &counts,
                normalization_factors,
                size_factors,
                method,
                geometric_means,
                control_genes,
            )?;
            write_normalized_counts_tsv(
                output,
                counts.gene_names(),
                counts.sample_names(),
                &normalized,
            )
        }
        Commands::Vst {
            counts,
            design,
            blind,
            normalization_factors,
            size_factors,
            observation_weights,
            method,
            geometric_means,
            control_genes,
            fit_type,
            nsub,
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let mut builder = DeseqBuilder::new()
                .size_factor_method(method.into())
                .fit_type(fit_type.into());
            builder = apply_cli_normalization_inputs(
                builder,
                &counts,
                normalization_factors,
                size_factors,
            )?;
            builder =
                apply_cli_size_factor_controls(builder, &counts, geometric_means, control_genes)?;
            if let Some(path) = observation_weights {
                builder = builder.observation_weights(read_cli_observation_weights(path, &counts)?);
            }
            let transformed = if blind {
                builder.blind_vst_glm_mu_auto(&counts, nsub)?.transformed
            } else {
                let design = design.ok_or_else(|| DeseqError::InvalidDimensions {
                    context: "VST design path".to_string(),
                    expected: 1,
                    actual: 0,
                })?;
                let design = read_cli_design_matrix(design, &counts)?;
                builder.vst_glm_mu_auto(&counts, &design, nsub)?.transformed
            };
            write_normalized_counts_tsv(
                output,
                counts.gene_names(),
                counts.sample_names(),
                &transformed,
            )
        }
        Commands::Rlog {
            counts,
            design,
            blind,
            normalization_factors,
            size_factors,
            observation_weights,
            method,
            geometric_means,
            control_genes,
            fit_type,
            frozen_intercept,
            rlog_prior_variance,
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let mut builder = DeseqBuilder::new()
                .size_factor_method(method.into())
                .fit_type(fit_type.into());
            builder = apply_cli_normalization_inputs(
                builder,
                &counts,
                normalization_factors,
                size_factors,
            )?;
            builder =
                apply_cli_size_factor_controls(builder, &counts, geometric_means, control_genes)?;
            if let Some(path) = observation_weights {
                builder = builder.observation_weights(read_cli_observation_weights(path, &counts)?);
            }
            let frozen_intercept = read_cli_frozen_intercept(frozen_intercept, &counts)?;
            let transformed = if blind {
                if let Some(frozen_intercept) = frozen_intercept {
                    let prior = required_cli_rlog_prior_variance(rlog_prior_variance)?;
                    let fit = builder.fit_map_dispersions_glm_mu(
                        &counts,
                        &crate::design::DesignMatrix::intercept_only(counts.n_samples())?,
                    )?;
                    fit.frozen_rlog(&counts, &frozen_intercept, prior)?
                        .transformed
                } else {
                    if rlog_prior_variance.is_some() {
                        return Err(cli_rlog_prior_without_frozen_intercept());
                    }
                    builder.blind_rlog_glm_mu(&counts)?.transformed
                }
            } else {
                let design = design.ok_or_else(|| DeseqError::InvalidDimensions {
                    context: "rlog design path".to_string(),
                    expected: 1,
                    actual: 0,
                })?;
                let design = read_cli_design_matrix(design, &counts)?;
                if let Some(frozen_intercept) = frozen_intercept {
                    let prior = required_cli_rlog_prior_variance(rlog_prior_variance)?;
                    let fit = builder.fit_map_dispersions_glm_mu(&counts, &design)?;
                    fit.frozen_rlog(&counts, &frozen_intercept, prior)?
                        .transformed
                } else {
                    if rlog_prior_variance.is_some() {
                        return Err(cli_rlog_prior_without_frozen_intercept());
                    }
                    builder.rlog_glm_mu(&counts, &design)?.transformed
                }
            };
            write_normalized_counts_tsv(
                output,
                counts.gene_names(),
                counts.sample_names(),
                &transformed,
            )
        }
        Commands::Wald {
            counts,
            design,
            normalization_factors,
            size_factors,
            observation_weights,
            method,
            geometric_means,
            control_genes,
            fit_type,
            coefficient,
            coefficient_name,
            contrast,
            contrast_name,
            contrast_positive,
            contrast_negative,
            contrast_positive_weight,
            contrast_negative_weight,
            contrast_factor,
            contrast_numerator,
            contrast_denominator,
            contrast_reference,
            contrast_sample_levels,
            lfc_threshold,
            alternative,
            use_t,
            t_degrees_of_freedom,
            t_degrees_of_freedom_file,
            disable_cooks_cutoff,
            cooks_cutoff,
            disable_independent_filtering,
            independent_filtering_alpha,
            independent_filtering_theta,
            result_column_metadata_output,
            result_table_metadata_output,
            independent_filter_metadata_output,
            independent_filter_num_rej_output,
            independent_filter_lowess_output,
            cooks_distance_output,
            cooks_replacement_metadata_output,
            cooks_replacement_row_metadata_output,
            cooks_replaced_counts_output,
            cooks_candidate_replacement_counts_output,
            cooks_outlier_cells_output,
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let design = read_cli_design_matrix(design, &counts)?;
            let mut builder = DeseqBuilder::new()
                .size_factor_method(method.into())
                .fit_type(fit_type.into())
                .wald_lfc_threshold(lfc_threshold, alternative.into());
            builder = apply_cli_wald_t_options(
                builder,
                &counts,
                use_t,
                t_degrees_of_freedom,
                t_degrees_of_freedom_file,
            )?;
            builder = apply_cli_result_options(
                builder,
                disable_cooks_cutoff,
                cooks_cutoff,
                disable_independent_filtering,
                independent_filtering_alpha,
                independent_filtering_theta,
            )?;
            builder = apply_cli_normalization_inputs(
                builder,
                &counts,
                normalization_factors,
                size_factors,
            )?;
            builder =
                apply_cli_size_factor_controls(builder, &counts, geometric_means, control_genes)?;
            if let Some(path) = observation_weights {
                builder = builder.observation_weights(read_cli_observation_weights(path, &counts)?);
            }
            let contrast_inputs = usize::from(coefficient.is_some())
                + usize::from(coefficient_name.is_some())
                + usize::from(contrast.is_some())
                + usize::from(contrast_name.is_some())
                + usize::from(contrast_positive.is_some() || contrast_negative.is_some())
                + usize::from(
                    contrast_factor.is_some()
                        || contrast_numerator.is_some()
                        || contrast_denominator.is_some()
                        || contrast_reference.is_some()
                        || contrast_sample_levels.is_some(),
                );
            if contrast_inputs > 1 {
                return Err(DeseqError::InvalidDimensions {
                    context: "Wald coefficient and contrast inputs".to_string(),
                    expected: 1,
                    actual: contrast_inputs,
                });
            }
            let cutoff = resolve_cooks_cutoff(
                builder.current_cooks_cutoff(),
                design.n_samples(),
                design.n_coefficients(),
            )?;
            let analysis = if let Some(contrast) = contrast {
                if let Some(cutoff) = cutoff {
                    cli_wald_replacement_output(
                        builder.fit_wald_glm_mu_contrast_with_cooks_replacement(
                            &counts,
                            &design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?,
                    )
                } else {
                    cli_fit_output(builder.fit_wald_glm_mu_contrast(&counts, &design, &contrast)?)
                }
            } else if let Some(contrast_name) = contrast_name {
                let contrast = ContrastSpec::coefficient_name(contrast_name);
                if let Some(cutoff) = cutoff {
                    cli_wald_replacement_output(
                        builder.fit_wald_glm_mu_contrast_spec_with_cooks_replacement(
                            &counts,
                            &design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?,
                    )
                } else {
                    cli_fit_output(
                        builder.fit_wald_glm_mu_contrast_spec(&counts, &design, &contrast)?,
                    )
                }
            } else if contrast_positive.is_some() || contrast_negative.is_some() {
                let contrast = ContrastSpec::list_with_values(
                    contrast_positive.unwrap_or_default(),
                    contrast_negative.unwrap_or_default(),
                    contrast_positive_weight,
                    contrast_negative_weight,
                );
                if let Some(cutoff) = cutoff {
                    cli_wald_replacement_output(
                        builder.fit_wald_glm_mu_contrast_spec_with_cooks_replacement(
                            &counts,
                            &design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?,
                    )
                } else {
                    cli_fit_output(
                        builder.fit_wald_glm_mu_contrast_spec(&counts, &design, &contrast)?,
                    )
                }
            } else if contrast_factor.is_some()
                || contrast_numerator.is_some()
                || contrast_denominator.is_some()
                || contrast_reference.is_some()
                || contrast_sample_levels.is_some()
            {
                let contrast = cli_factor_level_contrast(
                    contrast_factor,
                    contrast_numerator,
                    contrast_denominator,
                    contrast_reference.as_deref(),
                )?;
                if let Some(path) = contrast_sample_levels {
                    let levels = align_sample_levels_to_samples(
                        &read_sample_levels_tsv(path)?,
                        counts
                            .sample_names()
                            .ok_or_else(|| DeseqError::InvalidOptions {
                                reason: "count sample names are required to align sample levels"
                                    .to_string(),
                            })?,
                    )?;
                    let contrast = cli_factor_level_contrast_with_samples(&contrast, &levels)?;
                    if let Some(cutoff) = cutoff {
                        cli_wald_replacement_output(
                            builder.fit_wald_glm_mu_factor_level_contrast_with_cooks_replacement(
                                &counts,
                                &design,
                                contrast,
                                &CooksReplacementOptions::new(cutoff),
                            )?,
                        )
                    } else {
                        cli_fit_output(
                            builder.fit_wald_glm_mu_factor_level_contrast(
                                &counts, &design, contrast,
                            )?,
                        )
                    }
                } else if let Some(cutoff) = cutoff {
                    cli_wald_replacement_output(
                        builder.fit_wald_glm_mu_contrast_spec_with_cooks_replacement(
                            &counts,
                            &design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?,
                    )
                } else {
                    cli_fit_output(
                        builder.fit_wald_glm_mu_contrast_spec(&counts, &design, &contrast)?,
                    )
                }
            } else {
                let coefficient = match (coefficient, coefficient_name) {
                    (Some(coefficient), None) => coefficient,
                    (None, Some(name)) => design.coefficient_index(&name)?,
                    (None, None) => default_cli_coefficient(&design)?,
                    (Some(_), Some(_)) => unreachable!("checked above"),
                };
                if let Some(cutoff) = cutoff {
                    cli_wald_replacement_output(builder.fit_wald_glm_mu_with_cooks_replacement(
                        &counts,
                        &design,
                        coefficient,
                        &CooksReplacementOptions::new(cutoff),
                    )?)
                } else {
                    cli_fit_output(builder.fit_wald_glm_mu(&counts, &design, coefficient)?)
                }
            };
            let sidecars = CliCooksOutputPaths {
                cooks_distance: cooks_distance_output,
                replacement_metadata: cooks_replacement_metadata_output,
                replacement_row_metadata: cooks_replacement_row_metadata_output,
                replaced_counts: cooks_replaced_counts_output,
                candidate_replacement_counts: cooks_candidate_replacement_counts_output,
                outlier_cells: cooks_outlier_cells_output,
            };
            let result_sidecars = CliResultSidecarPaths {
                column_metadata: result_column_metadata_output,
                table_metadata: result_table_metadata_output,
                independent_filter_metadata: independent_filter_metadata_output,
                independent_filter_num_rej: independent_filter_num_rej_output,
                independent_filter_lowess: independent_filter_lowess_output,
            };
            write_cli_cooks_outputs(
                &sidecars,
                counts.gene_names(),
                counts.sample_names(),
                &analysis,
            )?;
            write_cli_result_sidecars(&result_sidecars, &analysis.results)?;
            write_deseq_results_tsv(output, &analysis.results)
        }
        Commands::Lrt {
            counts,
            design,
            reduced_design,
            normalization_factors,
            size_factors,
            observation_weights,
            method,
            geometric_means,
            control_genes,
            fit_type,
            coefficient,
            coefficient_name,
            contrast,
            contrast_name,
            contrast_positive,
            contrast_negative,
            contrast_positive_weight,
            contrast_negative_weight,
            contrast_factor,
            contrast_numerator,
            contrast_denominator,
            contrast_reference,
            contrast_sample_levels,
            disable_cooks_cutoff,
            cooks_cutoff,
            disable_independent_filtering,
            independent_filtering_alpha,
            independent_filtering_theta,
            result_column_metadata_output,
            result_table_metadata_output,
            independent_filter_metadata_output,
            independent_filter_num_rej_output,
            independent_filter_lowess_output,
            cooks_distance_output,
            cooks_replacement_metadata_output,
            cooks_replacement_row_metadata_output,
            cooks_replaced_counts_output,
            cooks_candidate_replacement_counts_output,
            cooks_outlier_cells_output,
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let design = read_cli_design_matrix(design, &counts)?;
            let reduced_design = read_cli_design_matrix(reduced_design, &counts)?;
            let mut builder = DeseqBuilder::new()
                .size_factor_method(method.into())
                .fit_type(fit_type.into());
            builder = apply_cli_result_options(
                builder,
                disable_cooks_cutoff,
                cooks_cutoff,
                disable_independent_filtering,
                independent_filtering_alpha,
                independent_filtering_theta,
            )?;
            builder = apply_cli_normalization_inputs(
                builder,
                &counts,
                normalization_factors,
                size_factors,
            )?;
            builder =
                apply_cli_size_factor_controls(builder, &counts, geometric_means, control_genes)?;
            if let Some(path) = observation_weights {
                builder = builder.observation_weights(read_cli_observation_weights(path, &counts)?);
            }
            let contrast_inputs = usize::from(coefficient.is_some())
                + usize::from(coefficient_name.is_some())
                + usize::from(contrast.is_some())
                + usize::from(contrast_name.is_some())
                + usize::from(contrast_positive.is_some() || contrast_negative.is_some())
                + usize::from(
                    contrast_factor.is_some()
                        || contrast_numerator.is_some()
                        || contrast_denominator.is_some()
                        || contrast_reference.is_some()
                        || contrast_sample_levels.is_some(),
                );
            if contrast_inputs > 1 {
                return Err(DeseqError::InvalidDimensions {
                    context: "LRT coefficient and contrast inputs".to_string(),
                    expected: 1,
                    actual: contrast_inputs,
                });
            }
            let cutoff = resolve_cooks_cutoff(
                builder.current_cooks_cutoff(),
                design.n_samples(),
                design.n_coefficients(),
            )?;
            let analysis = if let Some(contrast) = contrast {
                if let Some(cutoff) = cutoff {
                    cli_lrt_replacement_output(
                        builder.fit_lrt_glm_mu_contrast_with_cooks_replacement(
                            &counts,
                            &design,
                            &reduced_design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?,
                    )
                } else {
                    cli_fit_output(builder.fit_lrt_glm_mu_contrast(
                        &counts,
                        &design,
                        &reduced_design,
                        &contrast,
                    )?)
                }
            } else if let Some(contrast_name) = contrast_name {
                let contrast = ContrastSpec::coefficient_name(contrast_name);
                if let Some(cutoff) = cutoff {
                    cli_lrt_replacement_output(
                        builder.fit_lrt_glm_mu_contrast_spec_with_cooks_replacement(
                            &counts,
                            &design,
                            &reduced_design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?,
                    )
                } else {
                    cli_fit_output(builder.fit_lrt_glm_mu_contrast_spec(
                        &counts,
                        &design,
                        &reduced_design,
                        &contrast,
                    )?)
                }
            } else if contrast_positive.is_some() || contrast_negative.is_some() {
                let contrast = ContrastSpec::list_with_values(
                    contrast_positive.unwrap_or_default(),
                    contrast_negative.unwrap_or_default(),
                    contrast_positive_weight,
                    contrast_negative_weight,
                );
                if let Some(cutoff) = cutoff {
                    cli_lrt_replacement_output(
                        builder.fit_lrt_glm_mu_contrast_spec_with_cooks_replacement(
                            &counts,
                            &design,
                            &reduced_design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?,
                    )
                } else {
                    cli_fit_output(builder.fit_lrt_glm_mu_contrast_spec(
                        &counts,
                        &design,
                        &reduced_design,
                        &contrast,
                    )?)
                }
            } else if contrast_factor.is_some()
                || contrast_numerator.is_some()
                || contrast_denominator.is_some()
                || contrast_reference.is_some()
                || contrast_sample_levels.is_some()
            {
                let contrast = cli_factor_level_contrast(
                    contrast_factor,
                    contrast_numerator,
                    contrast_denominator,
                    contrast_reference.as_deref(),
                )?;
                if let Some(path) = contrast_sample_levels {
                    let levels = align_sample_levels_to_samples(
                        &read_sample_levels_tsv(path)?,
                        counts
                            .sample_names()
                            .ok_or_else(|| DeseqError::InvalidOptions {
                                reason: "count sample names are required to align sample levels"
                                    .to_string(),
                            })?,
                    )?;
                    let contrast = cli_factor_level_contrast_with_samples(&contrast, &levels)?;
                    if let Some(cutoff) = cutoff {
                        cli_lrt_replacement_output(
                            builder.fit_lrt_glm_mu_factor_level_contrast_with_cooks_replacement(
                                &counts,
                                &design,
                                &reduced_design,
                                contrast,
                                &CooksReplacementOptions::new(cutoff),
                            )?,
                        )
                    } else {
                        cli_fit_output(builder.fit_lrt_glm_mu_factor_level_contrast(
                            &counts,
                            &design,
                            &reduced_design,
                            contrast,
                        )?)
                    }
                } else if let Some(cutoff) = cutoff {
                    cli_lrt_replacement_output(
                        builder.fit_lrt_glm_mu_contrast_spec_with_cooks_replacement(
                            &counts,
                            &design,
                            &reduced_design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?,
                    )
                } else {
                    cli_fit_output(builder.fit_lrt_glm_mu_contrast_spec(
                        &counts,
                        &design,
                        &reduced_design,
                        &contrast,
                    )?)
                }
            } else {
                let coefficient = match (coefficient, coefficient_name) {
                    (Some(coefficient), None) => coefficient,
                    (None, Some(name)) => design.coefficient_index(&name)?,
                    (None, None) => default_cli_coefficient(&design)?,
                    _ => unreachable!("checked above"),
                };
                if let Some(cutoff) = cutoff {
                    cli_lrt_replacement_output(builder.fit_lrt_glm_mu_with_cooks_replacement(
                        &counts,
                        &design,
                        &reduced_design,
                        coefficient,
                        &CooksReplacementOptions::new(cutoff),
                    )?)
                } else {
                    cli_fit_output(builder.fit_lrt_glm_mu(
                        &counts,
                        &design,
                        &reduced_design,
                        coefficient,
                    )?)
                }
            };
            let sidecars = CliCooksOutputPaths {
                cooks_distance: cooks_distance_output,
                replacement_metadata: cooks_replacement_metadata_output,
                replacement_row_metadata: cooks_replacement_row_metadata_output,
                replaced_counts: cooks_replaced_counts_output,
                candidate_replacement_counts: cooks_candidate_replacement_counts_output,
                outlier_cells: cooks_outlier_cells_output,
            };
            let result_sidecars = CliResultSidecarPaths {
                column_metadata: result_column_metadata_output,
                table_metadata: result_table_metadata_output,
                independent_filter_metadata: independent_filter_metadata_output,
                independent_filter_num_rej: independent_filter_num_rej_output,
                independent_filter_lowess: independent_filter_lowess_output,
            };
            write_cli_cooks_outputs(
                &sidecars,
                counts.gene_names(),
                counts.sample_names(),
                &analysis,
            )?;
            write_cli_result_sidecars(&result_sidecars, &analysis.results)?;
            write_deseq_results_tsv(output, &analysis.results)
        }
    }
}

fn cli_normalized_counts(
    counts: &crate::core::CountMatrix,
    normalization_factors: Option<PathBuf>,
    size_factors: Option<PathBuf>,
    method: SizeFactorMethodArg,
    geometric_means: Option<PathBuf>,
    control_genes: Option<Vec<usize>>,
) -> Result<crate::matrix::RowMajorMatrix<f64>, DeseqError> {
    if let Some(path) = normalization_factors {
        if size_factors.is_some() {
            return Err(cli_conflicting_normalization_inputs());
        }
        let factors = read_cli_normalization_factors(path, counts)?;
        normalized_counts_with_factors(counts, &factors)
    } else if let Some(path) = size_factors {
        let size_factors = read_cli_size_factors(path, counts)?;
        normalized_counts(counts, &size_factors)
    } else {
        let geometric_means = read_cli_geometric_means(geometric_means, counts)?;
        let size_factors = estimate_size_factors_with_options(
            counts,
            method.into(),
            geometric_means.as_deref(),
            control_genes.as_deref(),
        )?;
        normalized_counts(counts, &size_factors)
    }
}

fn cli_fit_output((fit, results): (DeseqFit, DeseqResults)) -> CliAnalysisOutput {
    CliAnalysisOutput {
        results,
        cooks: fit.cooks.clone(),
        refit_plan: None,
    }
}

fn cli_wald_replacement_output(output: CooksReplacementWaldOutput) -> CliAnalysisOutput {
    CliAnalysisOutput {
        results: output.results,
        cooks: output.original_fit.cooks.clone(),
        refit_plan: Some(output.refit_plan),
    }
}

fn cli_lrt_replacement_output(output: CooksReplacementLrtOutput) -> CliAnalysisOutput {
    CliAnalysisOutput {
        results: output.results,
        cooks: output.original_fit.cooks.clone(),
        refit_plan: Some(output.refit_plan),
    }
}

fn write_cli_cooks_outputs(
    paths: &CliCooksOutputPaths,
    gene_names: Option<&[String]>,
    sample_names: Option<&[String]>,
    analysis: &CliAnalysisOutput,
) -> Result<(), DeseqError> {
    if paths.cooks_distance.is_some() {
        let cooks = analysis.cooks.as_ref().ok_or_else(|| DeseqError::InvalidOptions {
            reason: "Cook's diagnostic sidecar output requires a workflow that computes Cook's distances"
                .to_string(),
        })?;
        if let Some(path) = &paths.cooks_distance {
            write_cooks_distance_matrix_tsv(path, gene_names, sample_names, cooks)?;
        }
    }

    if paths.replacement_metadata.is_some()
        || paths.replacement_row_metadata.is_some()
        || paths.replaced_counts.is_some()
        || paths.candidate_replacement_counts.is_some()
        || paths.outlier_cells.is_some()
    {
        let refit_plan =
            analysis
                .refit_plan
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason:
                        "Cook's replacement sidecar output requires Cook's replacement/refit to run"
                            .to_string(),
                })?;
        if let Some(path) = &paths.replacement_metadata {
            write_cooks_replacement_metadata_tsv(path, refit_plan)?;
        }
        if let Some(path) = &paths.replacement_row_metadata {
            write_cooks_replacement_row_metadata_tsv(path, refit_plan)?;
        }
        if let Some(path) = &paths.replaced_counts {
            write_cooks_replaced_counts_tsv(path, refit_plan)?;
        }
        if let Some(path) = &paths.candidate_replacement_counts {
            write_cooks_candidate_replacement_counts_tsv(path, refit_plan)?;
        }
        if let Some(path) = &paths.outlier_cells {
            write_cooks_outlier_cells_tsv(path, refit_plan)?;
        }
    }
    Ok(())
}

fn write_cli_result_sidecars(
    paths: &CliResultSidecarPaths,
    results: &DeseqResults,
) -> Result<(), DeseqError> {
    if let Some(path) = &paths.column_metadata {
        write_deseq_result_column_metadata_tsv(path, results)?;
    }
    if let Some(path) = &paths.table_metadata {
        write_deseq_result_table_metadata_tsv(path, results)?;
    }
    if paths.independent_filter_metadata.is_some()
        || paths.independent_filter_num_rej.is_some()
        || paths.independent_filter_lowess.is_some()
    {
        let filtering =
            results
                .independent_filtering
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason:
                        "independent-filtering sidecar output requires independent filtering to run"
                            .to_string(),
                })?;
        if let Some(path) = &paths.independent_filter_metadata {
            write_independent_filter_metadata_tsv(path, filtering)?;
        }
        if let Some(path) = &paths.independent_filter_num_rej {
            write_independent_filter_num_rej_tsv(path, filtering)?;
        }
        if let Some(path) = &paths.independent_filter_lowess {
            write_independent_filter_lowess_tsv(path, filtering)?;
        }
    }
    Ok(())
}

fn apply_cli_normalization_inputs(
    builder: DeseqBuilder,
    counts: &crate::core::CountMatrix,
    normalization_factors: Option<PathBuf>,
    size_factors: Option<PathBuf>,
) -> Result<DeseqBuilder, DeseqError> {
    match (normalization_factors, size_factors) {
        (Some(normalization_factors), None) => Ok(builder.normalization_factors(
            read_cli_normalization_factors(normalization_factors, counts)?,
        )),
        (None, Some(size_factors)) => {
            Ok(builder.size_factors(read_cli_size_factors(size_factors, counts)?))
        }
        (None, None) => Ok(builder),
        (Some(_), Some(_)) => Err(cli_conflicting_normalization_inputs()),
    }
}

fn apply_cli_size_factor_controls(
    mut builder: DeseqBuilder,
    counts: &crate::core::CountMatrix,
    geometric_means: Option<PathBuf>,
    control_genes: Option<Vec<usize>>,
) -> Result<DeseqBuilder, DeseqError> {
    if let Some(geometric_means) = read_cli_geometric_means(geometric_means, counts)? {
        builder = builder.geometric_means(geometric_means);
    }
    if let Some(control_genes) = control_genes {
        builder = builder.control_genes(control_genes);
    }
    Ok(builder)
}

fn read_cli_geometric_means(
    path: Option<PathBuf>,
    counts: &crate::core::CountMatrix,
) -> Result<Option<Vec<f64>>, DeseqError> {
    path.map(|path| {
        align_gene_numeric_values_to_genes(
            &read_labeled_geometric_means_tsv(path)?,
            counts
                .gene_names()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "count gene names are required to align geometric means".to_string(),
                })?,
            "geometric-mean",
        )
    })
    .transpose()
}

fn read_cli_frozen_intercept(
    path: Option<PathBuf>,
    counts: &crate::core::CountMatrix,
) -> Result<Option<Vec<f64>>, DeseqError> {
    path.map(|path| {
        align_gene_numeric_values_to_genes(
            &read_labeled_gene_numeric_tsv(path, "rlog frozen intercept")?,
            counts
                .gene_names()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "count gene names are required to align rlog frozen intercepts"
                        .to_string(),
                })?,
            "rlog frozen intercept",
        )
    })
    .transpose()
}

fn required_cli_rlog_prior_variance(value: Option<f64>) -> Result<f64, DeseqError> {
    let value = value.ok_or_else(|| DeseqError::InvalidOptions {
        reason: "--rlog-prior-variance is required with --frozen-intercept".to_string(),
    })?;
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(DeseqError::InvalidOptions {
            reason: "--rlog-prior-variance must be positive and finite".to_string(),
        })
    }
}

fn cli_rlog_prior_without_frozen_intercept() -> DeseqError {
    DeseqError::InvalidOptions {
        reason: "--rlog-prior-variance requires --frozen-intercept".to_string(),
    }
}

fn read_cli_size_factors(
    path: impl Into<PathBuf>,
    counts: &crate::core::CountMatrix,
) -> Result<Vec<f64>, DeseqError> {
    align_sample_numeric_values_to_samples(
        &read_labeled_size_factors_tsv(path.into())?,
        counts
            .sample_names()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: "count sample names are required to align size factors".to_string(),
            })?,
        "size-factor",
    )
}

fn read_cli_normalization_factors(
    path: impl Into<PathBuf>,
    counts: &crate::core::CountMatrix,
) -> Result<crate::matrix::RowMajorMatrix<f64>, DeseqError> {
    align_labeled_assay_matrix_to_counts(
        read_labeled_normalization_factors_tsv(path.into())?,
        counts,
        "normalization factor",
    )
}

fn read_cli_observation_weights(
    path: impl Into<PathBuf>,
    counts: &crate::core::CountMatrix,
) -> Result<crate::matrix::RowMajorMatrix<f64>, DeseqError> {
    align_labeled_assay_matrix_to_counts(
        read_labeled_observation_weights_tsv(path.into())?,
        counts,
        "observation weight",
    )
}

fn apply_cli_result_options(
    mut builder: DeseqBuilder,
    disable_cooks_cutoff: bool,
    cooks_cutoff: Option<f64>,
    disable_independent_filtering: bool,
    independent_filtering_alpha: Option<f64>,
    independent_filtering_theta: Option<Vec<f64>>,
) -> Result<DeseqBuilder, DeseqError> {
    if disable_cooks_cutoff {
        if cooks_cutoff.is_some() {
            return Err(DeseqError::InvalidDimensions {
                context: "Cook's cutoff inputs".to_string(),
                expected: 1,
                actual: 2,
            });
        }
        builder = builder.cooks_cutoff(CooksCutoff::Disabled);
    } else if let Some(cutoff) = cooks_cutoff {
        builder = builder.cooks_cutoff_threshold(cutoff);
    }

    if disable_independent_filtering {
        builder = builder.disable_independent_filtering();
    }
    if let Some(alpha) = independent_filtering_alpha {
        builder = builder.independent_filtering_alpha(alpha);
    }
    if let Some(theta) = independent_filtering_theta {
        builder = builder.independent_filtering_theta(theta);
    }

    Ok(builder)
}

fn apply_cli_wald_t_options(
    builder: DeseqBuilder,
    counts: &crate::core::CountMatrix,
    use_t: bool,
    t_degrees_of_freedom: Option<f64>,
    t_degrees_of_freedom_file: Option<PathBuf>,
) -> Result<DeseqBuilder, DeseqError> {
    let requested = usize::from(use_t)
        + usize::from(t_degrees_of_freedom.is_some())
        + usize::from(t_degrees_of_freedom_file.is_some());
    if requested > 1 {
        return Err(DeseqError::InvalidDimensions {
            context: "Wald t p-value inputs".to_string(),
            expected: 1,
            actual: requested,
        });
    }

    if use_t {
        Ok(builder.wald_t_residual_degrees_of_freedom())
    } else if let Some(degrees_of_freedom) = t_degrees_of_freedom {
        Ok(builder.wald_t_degrees_of_freedom(degrees_of_freedom))
    } else if let Some(path) = t_degrees_of_freedom_file {
        Ok(
            builder.wald_t_per_gene_degrees_of_freedom(align_gene_numeric_values_to_genes(
                &read_labeled_wald_t_degrees_of_freedom_tsv(path)?,
                counts
                    .gene_names()
                    .ok_or_else(|| DeseqError::InvalidOptions {
                        reason: "count gene names are required to align Wald t degrees of freedom"
                            .to_string(),
                    })?,
                "Wald t degrees-of-freedom",
            )?),
        )
    } else {
        Ok(builder)
    }
}

fn cli_factor_level_contrast(
    factor: Option<String>,
    numerator: Option<String>,
    denominator: Option<String>,
    reference: Option<&str>,
) -> Result<ContrastSpec, DeseqError> {
    let supplied = usize::from(factor.is_some())
        + usize::from(numerator.is_some())
        + usize::from(denominator.is_some());
    let (Some(factor), Some(numerator), Some(denominator)) = (factor, numerator, denominator)
    else {
        return Err(DeseqError::InvalidDimensions {
            context: "factor-level contrast inputs".to_string(),
            expected: 3,
            actual: supplied,
        });
    };
    Ok(match reference {
        Some(reference) => {
            ContrastSpec::factor_level_with_reference(factor, numerator, denominator, reference)
        }
        None => ContrastSpec::factor_level(factor, numerator, denominator),
    })
}

fn cli_factor_level_contrast_with_samples<'a>(
    contrast: &'a ContrastSpec,
    sample_levels: &'a [String],
) -> Result<FactorLevelContrast<'a>, DeseqError> {
    match contrast {
        ContrastSpec::FactorLevel {
            factor,
            numerator,
            denominator,
            reference,
        } => Ok(FactorLevelContrast {
            factor,
            numerator,
            denominator,
            reference: reference.as_deref(),
            sample_levels,
        }),
        _ => Err(DeseqError::InvalidOptions {
            reason: "sample levels require a factor-level contrast".to_string(),
        }),
    }
}

fn read_cli_design_matrix(
    path: impl Into<PathBuf>,
    counts: &crate::core::CountMatrix,
) -> Result<crate::design::DesignMatrix, DeseqError> {
    align_design_matrix_to_samples(
        read_labeled_design_matrix_tsv(path.into())?,
        counts
            .sample_names()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: "count sample names are required to align design rows".to_string(),
            })?,
    )
}

fn cli_conflicting_normalization_inputs() -> DeseqError {
    DeseqError::InvalidDimensions {
        context: "normalization inputs".to_string(),
        expected: 1,
        actual: 2,
    }
}

fn default_cli_coefficient(design: &crate::design::DesignMatrix) -> Result<usize, DeseqError> {
    design
        .n_coefficients()
        .checked_sub(1)
        .ok_or_else(|| DeseqError::InvalidDimensions {
            context: "design matrix coefficients".to_string(),
            expected: 1,
            actual: 0,
        })
}
