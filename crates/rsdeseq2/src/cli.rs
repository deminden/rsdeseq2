use std::path::PathBuf;

use clap::{ArgAction, Parser, Subcommand, ValueEnum};

use crate::contrasts::ContrastSpec;
use crate::cooks::CooksReplacementOptions;
use crate::core::DeseqBuilder;
use crate::errors::DeseqError;
use crate::glm::WaldAlternative;
use crate::io::{
    read_count_matrix_tsv, read_design_matrix_tsv, read_geometric_means_tsv,
    read_normalization_factors_tsv, read_observation_weights_tsv, read_size_factors_tsv,
    read_wald_t_degrees_of_freedom_tsv, write_base_mean_tsv, write_deseq_results_tsv,
    write_normalized_counts_tsv, write_size_factors_tsv,
};
use crate::normalization::{
    base_mean, base_mean_with_weights, estimate_size_factors_with_options, normalized_counts,
    normalized_counts_with_factors,
};
use crate::options::{CooksCutoff, FitType, SizeFactorMethod};
use crate::results::resolve_cooks_cutoff;

/// Command-line arguments for the minimal `rsdeseq2` CLI.
#[derive(Debug, Parser)]
#[command(name = "rsdeseq2")]
#[command(about = "Early DESeq2-compatible Rust workflow stages")]
pub struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
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
        /// Numeric Wald contrast vector, comma-delimited in design-column order.
        #[arg(long, value_delimiter = ',')]
        contrast: Option<Vec<f64>>,
        /// Design coefficient name to report as a Wald contrast.
        #[arg(long)]
        contrast_name: Option<String>,
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
            let geometric_means = read_cli_geometric_means(geometric_means)?;
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
                let weights = read_observation_weights_tsv(path)?;
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
            builder = apply_cli_normalization_inputs(builder, normalization_factors, size_factors)?;
            builder = apply_cli_size_factor_controls(builder, geometric_means, control_genes)?;
            if let Some(path) = observation_weights {
                builder = builder.observation_weights(read_observation_weights_tsv(path)?);
            }
            let transformed = if blind {
                builder.blind_vst_glm_mu_auto(&counts, nsub)?.transformed
            } else {
                let design = design.ok_or_else(|| DeseqError::InvalidDimensions {
                    context: "VST design path".to_string(),
                    expected: 1,
                    actual: 0,
                })?;
                let design = read_design_matrix_tsv(design)?;
                builder.vst_glm_mu_auto(&counts, &design, nsub)?.transformed
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
            contrast,
            contrast_name,
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
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let design = read_design_matrix_tsv(design)?;
            let mut builder = DeseqBuilder::new()
                .size_factor_method(method.into())
                .fit_type(fit_type.into())
                .wald_lfc_threshold(lfc_threshold, alternative.into());
            builder = apply_cli_wald_t_options(
                builder,
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
            builder = apply_cli_normalization_inputs(builder, normalization_factors, size_factors)?;
            builder = apply_cli_size_factor_controls(builder, geometric_means, control_genes)?;
            if let Some(path) = observation_weights {
                builder = builder.observation_weights(read_observation_weights_tsv(path)?);
            }
            let contrast_inputs = usize::from(coefficient.is_some())
                + usize::from(contrast.is_some())
                + usize::from(contrast_name.is_some());
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
            let results = if let Some(contrast) = contrast {
                if let Some(cutoff) = cutoff {
                    builder
                        .fit_wald_glm_mu_contrast_with_cooks_replacement(
                            &counts,
                            &design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?
                        .results
                } else {
                    builder
                        .fit_wald_glm_mu_contrast(&counts, &design, &contrast)?
                        .1
                }
            } else if let Some(contrast_name) = contrast_name {
                let contrast = ContrastSpec::coefficient_name(contrast_name);
                if let Some(cutoff) = cutoff {
                    builder
                        .fit_wald_glm_mu_contrast_spec_with_cooks_replacement(
                            &counts,
                            &design,
                            &contrast,
                            &CooksReplacementOptions::new(cutoff),
                        )?
                        .results
                } else {
                    builder
                        .fit_wald_glm_mu_contrast_spec(&counts, &design, &contrast)?
                        .1
                }
            } else {
                let coefficient = if let Some(coefficient) = coefficient {
                    coefficient
                } else {
                    default_cli_coefficient(&design)?
                };
                if let Some(cutoff) = cutoff {
                    builder
                        .fit_wald_glm_mu_with_cooks_replacement(
                            &counts,
                            &design,
                            coefficient,
                            &CooksReplacementOptions::new(cutoff),
                        )?
                        .results
                } else {
                    builder.fit_wald_glm_mu(&counts, &design, coefficient)?.1
                }
            };
            write_deseq_results_tsv(output, &results)
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
            disable_cooks_cutoff,
            cooks_cutoff,
            disable_independent_filtering,
            independent_filtering_alpha,
            independent_filtering_theta,
            output,
        } => {
            let counts = read_count_matrix_tsv(counts)?;
            let design = read_design_matrix_tsv(design)?;
            let reduced_design = read_design_matrix_tsv(reduced_design)?;
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
            builder = apply_cli_normalization_inputs(builder, normalization_factors, size_factors)?;
            builder = apply_cli_size_factor_controls(builder, geometric_means, control_genes)?;
            if let Some(path) = observation_weights {
                builder = builder.observation_weights(read_observation_weights_tsv(path)?);
            }
            let coefficient = coefficient.unwrap_or(default_cli_coefficient(&design)?);
            let cutoff = resolve_cooks_cutoff(
                builder.current_cooks_cutoff(),
                design.n_samples(),
                design.n_coefficients(),
            )?;
            let results = if let Some(cutoff) = cutoff {
                builder
                    .fit_lrt_glm_mu_with_cooks_replacement(
                        &counts,
                        &design,
                        &reduced_design,
                        coefficient,
                        &CooksReplacementOptions::new(cutoff),
                    )?
                    .results
            } else {
                builder
                    .fit_lrt_glm_mu(&counts, &design, &reduced_design, coefficient)?
                    .1
            };
            write_deseq_results_tsv(output, &results)
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
        let factors = read_normalization_factors_tsv(path)?;
        normalized_counts_with_factors(counts, &factors)
    } else if let Some(path) = size_factors {
        let size_factors = read_size_factors_tsv(path)?;
        normalized_counts(counts, &size_factors)
    } else {
        let geometric_means = read_cli_geometric_means(geometric_means)?;
        let size_factors = estimate_size_factors_with_options(
            counts,
            method.into(),
            geometric_means.as_deref(),
            control_genes.as_deref(),
        )?;
        normalized_counts(counts, &size_factors)
    }
}

fn apply_cli_normalization_inputs(
    builder: DeseqBuilder,
    normalization_factors: Option<PathBuf>,
    size_factors: Option<PathBuf>,
) -> Result<DeseqBuilder, DeseqError> {
    match (normalization_factors, size_factors) {
        (Some(normalization_factors), None) => {
            Ok(builder
                .normalization_factors(read_normalization_factors_tsv(normalization_factors)?))
        }
        (None, Some(size_factors)) => {
            Ok(builder.size_factors(read_size_factors_tsv(size_factors)?))
        }
        (None, None) => Ok(builder),
        (Some(_), Some(_)) => Err(cli_conflicting_normalization_inputs()),
    }
}

fn apply_cli_size_factor_controls(
    mut builder: DeseqBuilder,
    geometric_means: Option<PathBuf>,
    control_genes: Option<Vec<usize>>,
) -> Result<DeseqBuilder, DeseqError> {
    if let Some(geometric_means) = read_cli_geometric_means(geometric_means)? {
        builder = builder.geometric_means(geometric_means);
    }
    if let Some(control_genes) = control_genes {
        builder = builder.control_genes(control_genes);
    }
    Ok(builder)
}

fn read_cli_geometric_means(path: Option<PathBuf>) -> Result<Option<Vec<f64>>, DeseqError> {
    path.map(read_geometric_means_tsv).transpose()
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
        Ok(builder.wald_t_per_gene_degrees_of_freedom(read_wald_t_degrees_of_freedom_tsv(path)?))
    } else {
        Ok(builder)
    }
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
