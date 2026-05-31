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
    gene_names: Option<&[String]>,
    analysis: &CliAnalysisOutput,
) -> Result<(), DeseqError> {
    let results = &analysis.results;
    // Sidecars are optional, but each requested file must correspond to data
    // produced by the selected workflow rather than silently exporting empties.
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
    if let Some(path) = &paths.fit_diagnostics {
        let fit = analysis
            .fit
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: "fit diagnostics sidecar output requires a native fit workflow".to_string(),
            })?;
        write_deseq_mcols_diagnostics_tsv(path, gene_names, &fit.deseq2_mcols_diagnostics())?;
    }
    if let Some(path) = &paths.refit_diagnostics {
        let refit = analysis
            .refit
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: "replacement refit diagnostics sidecar output requires rows to be refit"
                    .to_string(),
            })?;
        write_deseq_mcols_diagnostics_tsv(path, gene_names, &refit.deseq2_mcols_diagnostics())?;
    }
    if let Some(path) = &paths.fit_beta {
        let fit = require_cli_fit(analysis, "fit beta sidecar output")?;
        let beta = fit
            .beta
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: "fit beta sidecar output requires GLM beta estimates".to_string(),
            })?;
        write_optional_numeric_matrix_tsv(path, gene_names, fit_coefficient_names(fit), beta)?;
    }
    if let Some(path) = &paths.fit_beta_se {
        let fit = require_cli_fit(analysis, "fit beta standard-error sidecar output")?;
        let beta_se = fit
            .beta_se
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: "fit beta standard-error sidecar output requires GLM beta standard errors"
                    .to_string(),
            })?;
        write_optional_numeric_matrix_tsv(path, gene_names, fit_coefficient_names(fit), beta_se)?;
    }
    if let Some(path) = &paths.fit_beta_optim_start {
        let fit = require_cli_fit(analysis, "fit optimizer-start beta sidecar output")?;
        let beta_optim_start =
            fit.beta_optim_start
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason:
                        "fit optimizer-start beta sidecar output requires GLM fallback diagnostics"
                            .to_string(),
                })?;
        write_optional_numeric_matrix_tsv(
            path,
            gene_names,
            fit_coefficient_names(fit),
            beta_optim_start,
        )?;
    }
    if let Some(path) = &paths.refit_beta {
        let refit = require_cli_refit(analysis, "replacement refit beta sidecar output")?;
        let beta = refit
            .beta
            .as_ref()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: "replacement refit beta sidecar output requires GLM beta estimates"
                    .to_string(),
            })?;
        write_optional_numeric_matrix_tsv(path, gene_names, fit_coefficient_names(refit), beta)?;
    }
    if let Some(path) = &paths.refit_beta_se {
        let refit = require_cli_refit(
            analysis,
            "replacement refit beta standard-error sidecar output",
        )?;
        let beta_se = refit.beta_se.as_ref().ok_or_else(|| DeseqError::InvalidOptions {
            reason:
                "replacement refit beta standard-error sidecar output requires GLM beta standard errors"
                    .to_string(),
        })?;
        write_optional_numeric_matrix_tsv(path, gene_names, fit_coefficient_names(refit), beta_se)?;
    }
    if let Some(path) = &paths.refit_beta_optim_start {
        let refit = require_cli_refit(
            analysis,
            "replacement refit optimizer-start beta sidecar output",
        )?;
        let beta_optim_start =
            refit
                .beta_optim_start
                .as_ref()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "replacement refit optimizer-start beta sidecar output requires GLM fallback diagnostics"
                        .to_string(),
                })?;
        write_optional_numeric_matrix_tsv(
            path,
            gene_names,
            fit_coefficient_names(refit),
            beta_optim_start,
        )?;
    }
    Ok(())
}
