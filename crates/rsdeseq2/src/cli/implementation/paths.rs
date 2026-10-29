// Explicit path inventory for collision checks and output staging.
impl Commands {
    fn paths_mut(&mut self) -> (Vec<&mut PathBuf>, Vec<&mut PathBuf>) {
        match self {
            Self::SizeFactors {
                counts,
                geometric_means,
                output,
                ..
            } => {
                let mut inputs = Vec::new();
                let mut outputs = Vec::new();
                inputs.push(counts);
                inputs.extend(geometric_means.iter_mut());
                outputs.push(output);
                (inputs, outputs)
            }
            Self::BaseMean {
                counts,
                normalization_factors,
                size_factors,
                observation_weights,
                geometric_means,
                output,
                ..
            } => {
                let mut inputs = Vec::new();
                let mut outputs = Vec::new();
                inputs.push(counts);
                inputs.extend(normalization_factors.iter_mut());
                inputs.extend(size_factors.iter_mut());
                inputs.extend(observation_weights.iter_mut());
                inputs.extend(geometric_means.iter_mut());
                outputs.push(output);
                (inputs, outputs)
            }
            Self::NormalizedCounts {
                counts,
                normalization_factors,
                size_factors,
                geometric_means,
                output,
                ..
            } => {
                let mut inputs = Vec::new();
                let mut outputs = Vec::new();
                inputs.push(counts);
                inputs.extend(normalization_factors.iter_mut());
                inputs.extend(size_factors.iter_mut());
                inputs.extend(geometric_means.iter_mut());
                outputs.push(output);
                (inputs, outputs)
            }
            Self::Vst {
                counts,
                design,
                normalization_factors,
                size_factors,
                observation_weights,
                geometric_means,
                output,
                ..
            } => {
                let mut inputs = Vec::new();
                let mut outputs = Vec::new();
                inputs.push(counts);
                inputs.extend(design.iter_mut());
                inputs.extend(normalization_factors.iter_mut());
                inputs.extend(size_factors.iter_mut());
                inputs.extend(observation_weights.iter_mut());
                inputs.extend(geometric_means.iter_mut());
                outputs.push(output);
                (inputs, outputs)
            }
            Self::Rlog {
                counts,
                design,
                normalization_factors,
                size_factors,
                observation_weights,
                geometric_means,
                frozen_intercept,
                output,
                ..
            } => {
                let mut inputs = Vec::new();
                let mut outputs = Vec::new();
                inputs.push(counts);
                inputs.extend(design.iter_mut());
                inputs.extend(normalization_factors.iter_mut());
                inputs.extend(size_factors.iter_mut());
                inputs.extend(observation_weights.iter_mut());
                inputs.extend(geometric_means.iter_mut());
                inputs.extend(frozen_intercept.iter_mut());
                outputs.push(output);
                (inputs, outputs)
            }
            Self::Wald {
                counts,
                design,
                normalization_factors,
                size_factors,
                observation_weights,
                geometric_means,
                beta_prior_expanded_design,
                beta_prior_dispersions,
                beta_prior_base_mean,
                beta_prior_disp_fit,
                beta_prior_sample_levels,
                beta_prior_additive_sample_levels,
                beta_prior_additive_numeric_values,
                contrast_sample_levels,
                t_degrees_of_freedom_file,
                result_column_metadata_output,
                result_table_metadata_output,
                independent_filter_metadata_output,
                independent_filter_num_rej_output,
                independent_filter_lowess_output,
                fit_diagnostics_output,
                refit_diagnostics_output,
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
                ..
            } => {
                let mut inputs = Vec::new();
                let mut outputs = Vec::new();
                inputs.push(counts);
                inputs.push(design);
                inputs.extend(normalization_factors.iter_mut());
                inputs.extend(size_factors.iter_mut());
                inputs.extend(observation_weights.iter_mut());
                inputs.extend(geometric_means.iter_mut());
                inputs.extend(beta_prior_expanded_design.iter_mut());
                inputs.extend(beta_prior_dispersions.iter_mut());
                inputs.extend(beta_prior_base_mean.iter_mut());
                inputs.extend(beta_prior_disp_fit.iter_mut());
                inputs.extend(beta_prior_sample_levels.iter_mut());
                if let Some(paths) = beta_prior_additive_sample_levels {
                    inputs.extend(paths.iter_mut());
                }
                if let Some(paths) = beta_prior_additive_numeric_values {
                    inputs.extend(paths.iter_mut());
                }
                inputs.extend(contrast_sample_levels.iter_mut());
                inputs.extend(t_degrees_of_freedom_file.iter_mut());
                outputs.extend(result_column_metadata_output.iter_mut());
                outputs.extend(result_table_metadata_output.iter_mut());
                outputs.extend(independent_filter_metadata_output.iter_mut());
                outputs.extend(independent_filter_num_rej_output.iter_mut());
                outputs.extend(independent_filter_lowess_output.iter_mut());
                outputs.extend(fit_diagnostics_output.iter_mut());
                outputs.extend(refit_diagnostics_output.iter_mut());
                outputs.extend(fit_beta_output.iter_mut());
                outputs.extend(fit_beta_se_output.iter_mut());
                outputs.extend(fit_beta_optim_start_output.iter_mut());
                outputs.extend(refit_beta_output.iter_mut());
                outputs.extend(refit_beta_se_output.iter_mut());
                outputs.extend(refit_beta_optim_start_output.iter_mut());
                outputs.extend(cooks_distance_output.iter_mut());
                outputs.extend(cooks_replacement_metadata_output.iter_mut());
                outputs.extend(cooks_replacement_row_metadata_output.iter_mut());
                outputs.extend(cooks_replaced_counts_output.iter_mut());
                outputs.extend(cooks_candidate_replacement_counts_output.iter_mut());
                outputs.extend(cooks_outlier_cells_output.iter_mut());
                outputs.push(output);
                (inputs, outputs)
            }
            Self::Lrt {
                counts,
                design,
                reduced_design,
                normalization_factors,
                size_factors,
                observation_weights,
                geometric_means,
                contrast_sample_levels,
                result_column_metadata_output,
                result_table_metadata_output,
                independent_filter_metadata_output,
                independent_filter_num_rej_output,
                independent_filter_lowess_output,
                fit_diagnostics_output,
                refit_diagnostics_output,
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
                ..
            } => {
                let mut inputs = Vec::new();
                let mut outputs = Vec::new();
                inputs.push(counts);
                inputs.push(design);
                inputs.push(reduced_design);
                inputs.extend(normalization_factors.iter_mut());
                inputs.extend(size_factors.iter_mut());
                inputs.extend(observation_weights.iter_mut());
                inputs.extend(geometric_means.iter_mut());
                inputs.extend(contrast_sample_levels.iter_mut());
                outputs.extend(result_column_metadata_output.iter_mut());
                outputs.extend(result_table_metadata_output.iter_mut());
                outputs.extend(independent_filter_metadata_output.iter_mut());
                outputs.extend(independent_filter_num_rej_output.iter_mut());
                outputs.extend(independent_filter_lowess_output.iter_mut());
                outputs.extend(fit_diagnostics_output.iter_mut());
                outputs.extend(refit_diagnostics_output.iter_mut());
                outputs.extend(fit_beta_output.iter_mut());
                outputs.extend(fit_beta_se_output.iter_mut());
                outputs.extend(fit_beta_optim_start_output.iter_mut());
                outputs.extend(refit_beta_output.iter_mut());
                outputs.extend(refit_beta_se_output.iter_mut());
                outputs.extend(refit_beta_optim_start_output.iter_mut());
                outputs.extend(cooks_distance_output.iter_mut());
                outputs.extend(cooks_replacement_metadata_output.iter_mut());
                outputs.extend(cooks_replacement_row_metadata_output.iter_mut());
                outputs.extend(cooks_replaced_counts_output.iter_mut());
                outputs.extend(cooks_candidate_replacement_counts_output.iter_mut());
                outputs.extend(cooks_outlier_cells_output.iter_mut());
                outputs.push(output);
                (inputs, outputs)
            }
        }
    }
}
