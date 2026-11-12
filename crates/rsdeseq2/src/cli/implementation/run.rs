use super::*;

/// Parse process arguments and run the CLI.
pub fn run_cli() -> Result<(), DeseqError> {
    run(Cli::parse())
}

pub(super) fn run(mut cli: Cli) -> Result<(), DeseqError> {
    let (inputs, mut outputs) = cli.command.paths_mut();
    let transaction =
        crate::cli::output_transaction::OutputTransaction::prepare(&inputs, &mut outputs)?;
    run_staged(cli)?;
    transaction.publish()
}

pub(super) fn run_staged(cli: Cli) -> Result<(), DeseqError> {
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
                .automatic_dispersion_means()
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
                .automatic_dispersion_means()
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
        command @ Commands::Wald { .. } => run_wald::run_wald(command),
        command @ Commands::Lrt { .. } => run_lrt::run_lrt(command),
    }
}
