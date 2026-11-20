use super::*;

// Execute the Wald command after output staging.
pub(super) fn run_wald(args: WaldArgs) -> Result<(), DeseqError> {
    let WaldArgs {
        counts,
        design,
        normalization,
        target,
        beta_prior,
        contrast,
        test,
        filtering,
        outputs,
    } = args;
    let counts = read_count_matrix_tsv(counts)?;
    let design = read_cli_design_matrix(design, &counts)?;
    let prior_normalization = (beta_prior.supplied_count() > 0).then(|| normalization.clone());
    let mut builder = DeseqBuilder::new()
        .automatic_dispersion_means()
        .size_factor_method(normalization.method.into())
        .fit_type(normalization.fit_type.into())
        .wald_lfc_threshold(test.lfc_threshold, test.alternative.into());
    builder = apply_cli_wald_t_options(
        builder,
        &counts,
        test.use_t,
        test.t_degrees_of_freedom,
        test.t_degrees_of_freedom_file,
    )?;
    builder = configure_cli_results(builder, &counts, normalization, filtering)?;
    validate_cli_result_target("Wald", &target, &contrast)?;
    let cutoff = resolve_cooks_cutoff(
        builder.current_cooks_cutoff(),
        design.n_samples(),
        design.n_coefficients(),
    )?;
    let analysis = if let Some(normalization) = prior_normalization {
        run_supplied_beta_prior(
            &counts,
            &design,
            beta_prior,
            normalization,
            target,
            contrast,
            cutoff,
        )?
    } else {
        run_cli_model(&builder, &counts, &design, None, target, contrast, cutoff)?
    };
    write_cli_analysis_outputs(outputs, &counts, &analysis)
}

fn run_supplied_beta_prior(
    counts: &crate::matrix::CountMatrix,
    design: &crate::design::DesignMatrix,
    beta_prior: BetaPriorArgs,
    normalization: NormalizationArgs,
    target: ResultTargetArgs,
    contrast: ContrastArgs,
    cutoff: Option<f64>,
) -> Result<CliAnalysisOutput, DeseqError> {
    let beta_prior_inputs = beta_prior.supplied_count();
    let BetaPriorArgs {
        beta_prior_expanded_design,
        beta_prior_coefficient_groups,
        beta_prior_dispersions,
        beta_prior_base_mean,
        beta_prior_disp_fit,
        beta_prior_factor,
        beta_prior_reference,
        beta_prior_sample_levels,
        beta_prior_additive_factors,
        beta_prior_additive_references,
        beta_prior_additive_sample_levels,
        beta_prior_additive_numeric,
        beta_prior_additive_numeric_values,
    } = beta_prior;
    let beta_prior_matrix_inputs =
        beta_prior_expanded_design.is_some() || beta_prior_coefficient_groups.is_some();
    let beta_prior_factor_inputs = beta_prior_factor.is_some()
        || beta_prior_reference.is_some()
        || beta_prior_sample_levels.is_some();
    let beta_prior_additive_inputs = beta_prior_additive_factors.is_some()
        || beta_prior_additive_references.is_some()
        || beta_prior_additive_sample_levels.is_some()
        || beta_prior_additive_numeric.is_some()
        || beta_prior_additive_numeric_values.is_some();
    let beta_prior_design_routes = [
        beta_prior_matrix_inputs,
        beta_prior_factor_inputs,
        beta_prior_additive_inputs,
    ]
    .into_iter()
    .filter(|present| *present)
    .count();
    if beta_prior_design_routes > 1 {
        return Err(DeseqError::InvalidOptions {
            reason:
                "beta-prior expanded matrix, one-factor, and additive-factor inputs are mutually exclusive"
                    .to_string(),
        });
    }
    if beta_prior_matrix_inputs && beta_prior_inputs != 5 {
        return Err(DeseqError::InvalidDimensions {
            context: "beta-prior expanded Wald inputs".to_string(),
            expected: 5,
            actual: beta_prior_inputs,
        });
    }
    if beta_prior_factor_inputs && beta_prior_inputs != 6 {
        return Err(DeseqError::InvalidDimensions {
            context: "beta-prior factor Wald inputs".to_string(),
            expected: 6,
            actual: beta_prior_inputs,
        });
    }
    if beta_prior_additive_inputs {
        let additive_common_inputs = [
            beta_prior_dispersions.is_some(),
            beta_prior_base_mean.is_some(),
            beta_prior_disp_fit.is_some(),
        ]
        .into_iter()
        .filter(|present| *present)
        .count();
        let additive_factor_inputs = [
            beta_prior_additive_factors.is_some(),
            beta_prior_additive_references.is_some(),
            beta_prior_additive_sample_levels.is_some(),
        ]
        .into_iter()
        .filter(|present| *present)
        .count();
        let additive_numeric_inputs = [
            beta_prior_additive_numeric.is_some(),
            beta_prior_additive_numeric_values.is_some(),
        ]
        .into_iter()
        .filter(|present| *present)
        .count();
        if additive_common_inputs != 3
            || matches!(additive_factor_inputs, 1 | 2)
            || additive_numeric_inputs == 1
            || additive_factor_inputs + additive_numeric_inputs == 0
        {
            return Err(DeseqError::InvalidDimensions {
                context: "beta-prior additive-factor Wald inputs".to_string(),
                expected: 3,
                actual: additive_common_inputs,
            });
        }
    }
    if beta_prior_additive_inputs && beta_prior_inputs < 5 {
        return Err(DeseqError::InvalidDimensions {
            context: "beta-prior additive-factor Wald inputs".to_string(),
            expected: 5,
            actual: beta_prior_inputs,
        });
    }
    if contrast.has_factor() {
        return Err(DeseqError::InvalidOptions {
            reason:
                "beta-prior expanded Wald CLI currently accepts coefficient, named, list, or numeric contrasts"
                    .to_string(),
        });
    }
    if beta_prior_factor_inputs {
        cli_factor_beta_prior_wald_analysis(
            counts,
            design,
            required(beta_prior_factor, "--beta-prior-factor")?,
            required(beta_prior_reference, "--beta-prior-reference")?,
            required(beta_prior_sample_levels, "--beta-prior-sample-levels")?,
            required(beta_prior_dispersions, "--beta-prior-dispersions")?,
            required(beta_prior_base_mean, "--beta-prior-base-mean")?,
            required(beta_prior_disp_fit, "--beta-prior-disp-fit")?,
            normalization,
            target,
            contrast,
            cutoff,
        )
    } else if beta_prior_additive_inputs {
        cli_additive_beta_prior_wald_analysis(
            counts,
            design,
            beta_prior_additive_factors.unwrap_or_default(),
            beta_prior_additive_references.unwrap_or_default(),
            beta_prior_additive_sample_levels.unwrap_or_default(),
            beta_prior_additive_numeric.unwrap_or_default(),
            beta_prior_additive_numeric_values.unwrap_or_default(),
            required(beta_prior_dispersions, "--beta-prior-dispersions")?,
            required(beta_prior_base_mean, "--beta-prior-base-mean")?,
            required(beta_prior_disp_fit, "--beta-prior-disp-fit")?,
            normalization,
            target,
            contrast,
            cutoff,
        )
    } else {
        cli_expanded_beta_prior_wald_analysis(
            counts,
            design,
            required(beta_prior_expanded_design, "--beta-prior-expanded-design")?,
            &required(
                beta_prior_coefficient_groups,
                "--beta-prior-coefficient-groups",
            )?,
            required(beta_prior_dispersions, "--beta-prior-dispersions")?,
            required(beta_prior_base_mean, "--beta-prior-base-mean")?,
            required(beta_prior_disp_fit, "--beta-prior-disp-fit")?,
            normalization,
            target,
            contrast,
            cutoff,
        )
    }
}
