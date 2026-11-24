use super::*;
use crate::errors::{ErrorKind, missing_prerequisite};

pub(super) fn write_cli_analysis_outputs(
    outputs: DiagnosticOutputArgs,
    counts: &crate::matrix::CountMatrix,
    analysis: &CliAnalysisOutput,
) -> Result<(), DeseqError> {
    let DiagnosticOutputArgs {
        result_column_metadata_output,
        result_table_metadata_output,
        independent_filter_metadata_output,
        independent_filter_num_rej_output,
        independent_filter_lowess_output,
        fit_diagnostics_output,
        refit_diagnostics_output,
        fit_diagnostics_summary_output,
        fit_beta_output,
        fit_beta_se_output,
        fit_beta_optim_start_output,
        refit_beta_output,
        refit_beta_se_output,
        refit_beta_optim_start_output,
        cooks_distance_output,
        cooks_replacement_metadata_output,
        cooks_replacement_row_metadata_output,
        cooks_replaced_counts_output,
        cooks_candidate_replacement_counts_output,
        cooks_outlier_cells_output,
        output,
    } = outputs;
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
        fit_diagnostics: fit_diagnostics_output,
        refit_diagnostics: refit_diagnostics_output,
        fit_beta: fit_beta_output,
        fit_beta_se: fit_beta_se_output,
        fit_beta_optim_start: fit_beta_optim_start_output,
        refit_beta: refit_beta_output,
        refit_beta_se: refit_beta_se_output,
        refit_beta_optim_start: refit_beta_optim_start_output,
    };
    write_cli_cooks_outputs(
        &sidecars,
        counts.gene_names(),
        counts.sample_names(),
        analysis,
    )?;
    write_cli_result_sidecars(&result_sidecars, counts.gene_names(), analysis)?;
    if let Some(path) = fit_diagnostics_summary_output {
        write_cli_fit_diagnostics_summary(&path, analysis)?;
    }
    write_deseq_results_tsv(output, &analysis.results)
}

pub(super) fn write_cli_fit_diagnostics_summary(
    path: &std::path::Path,
    analysis: &CliAnalysisOutput,
) -> Result<(), DeseqError> {
    let fit = require_cli_fit(analysis, "fit diagnostics summary output")?;
    let mut writer = crate::io::TsvWriter::from_path(path)?;
    writer.write_record(["fit", "nonfinite_covariance_rows"])?;
    for (name, fit) in [("original", Some(fit)), ("refit", analysis.refit.as_ref())] {
        if let Some(fit) = fit {
            let count = fit
                .nonfinite_covariance_rows()
                .map_or_else(|| "NA".to_string(), |count| count.to_string());
            writer.write_record([name, &count])?;
        }
    }
    writer.flush()?;
    Ok(())
}

pub(super) fn write_cli_cooks_outputs(
    paths: &CliCooksOutputPaths,
    gene_names: Option<&[String]>,
    sample_names: Option<&[String]>,
    analysis: &CliAnalysisOutput,
) -> Result<(), DeseqError> {
    if paths.cooks_distance.is_some() {
        let cooks = analysis
            .cooks
            .as_ref()
            .or_else(|| analysis.fit.as_ref().and_then(|fit| fit.cooks.as_ref()))
            .ok_or_else(|| {
                missing_prerequisite(
                    ErrorKind::Options,
                    "Cook's diagnostic sidecar output",
                    "a workflow that computes Cook's distances",
                )
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
        let refit_plan = analysis.refit_plan.as_ref().ok_or_else(|| {
            missing_prerequisite(
                ErrorKind::Options,
                "Cook's replacement sidecar output",
                "Cook's replacement/refit to run",
            )
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

pub(super) fn write_cli_result_sidecars(
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
        let filtering = results.independent_filtering.as_ref().ok_or_else(|| {
            missing_prerequisite(
                ErrorKind::Options,
                "independent-filtering sidecar output",
                "independent filtering to run",
            )
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
        let fit = analysis.fit.as_ref().ok_or_else(|| {
            missing_prerequisite(
                ErrorKind::Options,
                "fit diagnostics sidecar output",
                "a native fit workflow",
            )
        })?;
        write_deseq_mcols_diagnostics_tsv(path, gene_names, &fit.deseq2_mcols_diagnostics())?;
    }
    if let Some(path) = &paths.refit_diagnostics {
        let refit = analysis.refit.as_ref().ok_or_else(|| {
            missing_prerequisite(
                ErrorKind::Options,
                "replacement refit diagnostics sidecar output",
                "rows to be refit",
            )
        })?;
        write_deseq_mcols_diagnostics_tsv(path, gene_names, &refit.deseq2_mcols_diagnostics())?;
    }
    if let Some(path) = &paths.fit_beta {
        let fit = require_cli_fit(analysis, "fit beta sidecar output")?;
        let beta = fit.beta.as_ref().ok_or_else(|| {
            missing_prerequisite(
                ErrorKind::Options,
                "fit beta sidecar output",
                "GLM beta estimates",
            )
        })?;
        write_optional_numeric_matrix_tsv(path, gene_names, fit_coefficient_names(fit), beta)?;
    }
    if let Some(path) = &paths.fit_beta_se {
        let fit = require_cli_fit(analysis, "fit beta standard-error sidecar output")?;
        let beta_se = fit.beta_se.as_ref().ok_or_else(|| {
            missing_prerequisite(
                ErrorKind::Options,
                "fit beta standard-error sidecar output",
                "GLM beta standard errors",
            )
        })?;
        write_optional_numeric_matrix_tsv(path, gene_names, fit_coefficient_names(fit), beta_se)?;
    }
    if let Some(path) = &paths.fit_beta_optim_start {
        let fit = require_cli_fit(analysis, "fit optimizer-start beta sidecar output")?;
        let beta_optim_start = fit.beta_optim_start.as_ref().ok_or_else(|| {
            missing_prerequisite(
                ErrorKind::Options,
                "fit optimizer-start beta sidecar output",
                "GLM fallback diagnostics",
            )
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
        let beta = refit.beta.as_ref().ok_or_else(|| {
            missing_prerequisite(
                ErrorKind::Options,
                "replacement refit beta sidecar output",
                "GLM beta estimates",
            )
        })?;
        write_optional_numeric_matrix_tsv(path, gene_names, fit_coefficient_names(refit), beta)?;
    }
    if let Some(path) = &paths.refit_beta_se {
        let refit = require_cli_refit(
            analysis,
            "replacement refit beta standard-error sidecar output",
        )?;
        let beta_se = refit.beta_se.as_ref().ok_or_else(|| {
            missing_prerequisite(
                ErrorKind::Options,
                "replacement refit beta standard-error sidecar output",
                "GLM beta standard errors",
            )
        })?;
        write_optional_numeric_matrix_tsv(path, gene_names, fit_coefficient_names(refit), beta_se)?;
    }
    if let Some(path) = &paths.refit_beta_optim_start {
        let refit = require_cli_refit(
            analysis,
            "replacement refit optimizer-start beta sidecar output",
        )?;
        let beta_optim_start = refit.beta_optim_start.as_ref().ok_or_else(|| {
            missing_prerequisite(
                ErrorKind::Options,
                "replacement refit optimizer-start beta sidecar output",
                "GLM fallback diagnostics",
            )
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
