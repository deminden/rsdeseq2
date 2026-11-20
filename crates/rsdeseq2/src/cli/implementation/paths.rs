use super::*;
// Explicit path inventory for collision checks and output staging.
impl Commands {
    pub(super) fn paths_mut(&mut self) -> (Vec<&mut PathBuf>, Vec<&mut PathBuf>) {
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
            Self::Wald(args) => {
                let mut inputs = vec![&mut args.counts, &mut args.design];
                inputs.extend(args.normalization.input_paths_mut());
                inputs.extend(args.beta_prior.input_paths_mut());
                inputs.extend(args.contrast.contrast_sample_levels.iter_mut());
                inputs.extend(args.test.t_degrees_of_freedom_file.iter_mut());
                (inputs, args.outputs.output_paths_mut())
            }
            Self::Lrt(args) => {
                let mut inputs = vec![&mut args.counts, &mut args.design, &mut args.reduced_design];
                inputs.extend(args.normalization.input_paths_mut());
                inputs.extend(args.contrast.contrast_sample_levels.iter_mut());
                (inputs, args.outputs.output_paths_mut())
            }
        }
    }
}

impl NormalizationArgs {
    fn input_paths_mut(&mut self) -> impl Iterator<Item = &mut PathBuf> {
        [
            self.normalization_factors.as_mut(),
            self.size_factors.as_mut(),
            self.observation_weights.as_mut(),
            self.geometric_means.as_mut(),
        ]
        .into_iter()
        .flatten()
    }
}

impl BetaPriorArgs {
    fn input_paths_mut(&mut self) -> impl Iterator<Item = &mut PathBuf> {
        [
            self.beta_prior_expanded_design.as_mut(),
            self.beta_prior_dispersions.as_mut(),
            self.beta_prior_base_mean.as_mut(),
            self.beta_prior_disp_fit.as_mut(),
            self.beta_prior_sample_levels.as_mut(),
        ]
        .into_iter()
        .flatten()
        .chain(self.beta_prior_additive_sample_levels.iter_mut().flatten())
        .chain(self.beta_prior_additive_numeric_values.iter_mut().flatten())
    }
}

impl DiagnosticOutputArgs {
    fn output_paths_mut(&mut self) -> Vec<&mut PathBuf> {
        [
            self.result_column_metadata_output.as_mut(),
            self.result_table_metadata_output.as_mut(),
            self.independent_filter_metadata_output.as_mut(),
            self.independent_filter_num_rej_output.as_mut(),
            self.independent_filter_lowess_output.as_mut(),
            self.fit_diagnostics_output.as_mut(),
            self.refit_diagnostics_output.as_mut(),
            self.fit_diagnostics_summary_output.as_mut(),
            self.fit_beta_output.as_mut(),
            self.fit_beta_se_output.as_mut(),
            self.fit_beta_optim_start_output.as_mut(),
            self.refit_beta_output.as_mut(),
            self.refit_beta_se_output.as_mut(),
            self.refit_beta_optim_start_output.as_mut(),
            self.cooks_distance_output.as_mut(),
            self.cooks_replacement_metadata_output.as_mut(),
            self.cooks_replacement_row_metadata_output.as_mut(),
            self.cooks_replaced_counts_output.as_mut(),
            self.cooks_candidate_replacement_counts_output.as_mut(),
            self.cooks_outlier_cells_output.as_mut(),
            Some(&mut self.output),
        ]
        .into_iter()
        .flatten()
        .collect()
    }
}
