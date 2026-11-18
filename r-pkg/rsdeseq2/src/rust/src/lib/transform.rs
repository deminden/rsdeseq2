// VST/rlog transformation export and frozen fit inputs.
use super::*;

/// # Safety
///
/// Non-null input pointers must reference aligned, readable buffers of the
/// lengths specified by the dimensions and the C bridge contract. Output
/// pointers must reference writable buffers of those lengths, uniquely
/// borrowed for this call. All buffers must remain valid until return.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_transform(
    n_genes: usize,
    n_samples: usize,
    counts_values: *const f64,
    size_factors_values: *const f64,
    normalization_factors_values: *const f64,
    observation_weights_values: *const f64,
    n_coefficients: usize,
    design_values: *const f64,
    fit_type: i32,
    transform_type: i32,
    blind: u8,
    nsub: usize,
    threads: usize,
    frozen_intercept_values: *const f64,
    sample_prior_variance: f64,
    frozen_trend_type: i32,
    frozen_trend_parameters: *const f64,
    frozen_disp_fit_values: *const f64,
    frozen_local_curve: *const f64,
    transformed_output: *mut f64,
    size_factors_output: *mut f64,
    intercept_output: *mut f64,
    prior_variance_output: *mut f64,
    fit_type_used_output: *mut i32,
) -> i32 {
    run_backend(|| {
        let fit_type_used =
            // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
            unsafe { output_slice(fit_type_used_output, 1, "fitted dispersion trend type")? };
        fit_type_used[0] = fit_type;
        // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
        let counts = unsafe { counts_from_r(n_genes, n_samples, counts_values)? };
        let design_values =
            // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
            unsafe { row_major_from_column_major(n_samples, n_coefficients, design_values)? };
        let design = DesignMatrix::from_row_major(
            n_samples,
            n_coefficients,
            design_values.into_values(),
            None,
        )
        .map_err(|error| error.to_string())?;
        let selected_fit_type = match fit_type {
            0 => FitType::Parametric,
            1 => FitType::Local,
            2 => FitType::Mean,
            _ => return Err("fitType must be parametric, local, or mean".into()),
        };
        let mut builder = DeseqBuilder::new()
            .automatic_dispersion_means()
            .fit_type(selected_fit_type)
            .threads(threads.max(1));
        if !size_factors_values.is_null() {
            builder = builder.size_factors(
                // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
                unsafe { input_slice(size_factors_values, n_samples, "size factors")? }.to_vec(),
            );
        }
        if !normalization_factors_values.is_null() {
            // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
            let factors = unsafe {
                row_major_from_column_major(n_genes, n_samples, normalization_factors_values)?
            };
            builder = builder.normalization_factors(factors);
        }
        if !observation_weights_values.is_null() {
            // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
            let weights = unsafe {
                row_major_from_column_major(n_genes, n_samples, observation_weights_values)?
            };
            builder = builder.observation_weights(weights);
        }

        let transformed = match transform_type {
            0 => {
                if nsub == 0 {
                    return Err("nsub must be a positive integer".into());
                }
                // DESeq2's `vst()` always fits the fast subset and errors if
                // fewer than `nsub` rows have ordinary normalized mean > 5.
                // The builder's `*_auto` helpers intentionally retain a
                // Rust-specific full-data fallback, so use the strict fast
                // path for the Bioconductor-compatible wrapper.
                let output = builder
                    .fast_vst_glm_mu(&counts, &design, nsub)
                    .map_err(|error| error.to_string())?;
                fit_type_used[0] =
                    fit_type_code(output.subset_fit.dispersion_trend.as_ref(), fit_type);
                // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
                unsafe {
                    output_slice(size_factors_output, n_samples, "size factors")?
                        .copy_from_slice(&output.subset_fit.size_factors);
                    output_slice(intercept_output, n_genes, "rlog intercepts")?.fill(f64::NAN);
                    output_slice(prior_variance_output, 1, "rlog prior variance")?[0] = f64::NAN;
                }
                output.transformed
            }
            1 => {
                let (fit, rlog) = if frozen_intercept_values.is_null() {
                    if !frozen_disp_fit_values.is_null() {
                        let mut fit = trend_stage_fit(&builder, &counts, &design)?;
                        // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
                        fit.disp_fit = Some(unsafe {
                            input_slice(frozen_disp_fit_values, n_genes, "frozen rlog dispFit")?
                                .to_vec()
                        });
                        let rlog = builder
                            .with_execution_context(|| fit.regularized_log_transform(&counts))
                            .and_then(|result| result)
                            .map_err(|error| error.to_string())?;
                        (fit, rlog)
                    } else {
                        let output = if blind != 0 {
                            builder.blind_rlog_glm_mu_with_fit(&counts)
                        } else {
                            builder.rlog_glm_mu_with_fit(&counts, &design)
                        }
                        .map_err(|error| error.to_string())?;
                        (output.fit, output.rlog)
                    }
                } else {
                    if !sample_prior_variance.is_finite() || sample_prior_variance <= 0.0 {
                        return Err("sample prior variance must be positive and finite with frozen intercepts".into());
                    }
                    // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
                    let frozen_intercept = unsafe {
                        input_slice(frozen_intercept_values, n_genes, "frozen rlog intercepts")?
                    };
                    let fit = trend_stage_fit(&builder, &counts, &design)?;
                    let rlog = builder
                        .with_execution_context(|| fit.frozen_rlog(&counts, frozen_intercept, sample_prior_variance))
                        .and_then(|result| result)
                        .map_err(|error| error.to_string())?;
                    (fit, rlog)
                };
                fit_type_used[0] = fit_type_code(fit.dispersion_trend.as_ref(), fit_type);
                // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
                unsafe {
                    output_slice(size_factors_output, n_samples, "size factors")?
                        .copy_from_slice(&fit.size_factors);
                    output_slice(intercept_output, n_genes, "rlog intercepts")?
                        .copy_from_slice(&rlog.intercept);
                    output_slice(prior_variance_output, 1, "rlog prior variance")?[0] =
                        rlog.sample_prior_variance;
                }
                rlog.transformed
            }
            2 => {
                let (transformed, output_size_factors) = if frozen_trend_type >= 0 {
                    let size_factors = if !normalization_factors_values.is_null() {
                        if size_factors_values.is_null() {
                            vec![1.0; n_samples]
                        } else {
                            // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
                            unsafe {
                                input_slice(size_factors_values, n_samples, "size factors")?
                                    .to_vec()
                            }
                        }
                    } else if !size_factors_values.is_null() {
                        // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
                        unsafe {
                            input_slice(size_factors_values, n_samples, "size factors")?.to_vec()
                        }
                    } else {
                        estimate_size_factors_with_options(
                            &counts,
                            SizeFactorMethod::Ratio,
                            None,
                            None,
                        )
                        .map_err(|error| error.to_string())?
                    };
                    let normalization_factors = if normalization_factors_values.is_null() {
                        None
                    } else {
                        // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
                        Some(unsafe {
                            row_major_from_column_major(
                                n_genes,
                                n_samples,
                                normalization_factors_values,
                            )?
                        })
                    };
                    let normalized = match &normalization_factors {
                        Some(factors) => normalized_counts_with_factors(&counts, factors),
                        None => normalized_counts(&counts, &size_factors),
                    }
                    .map_err(|error| error.to_string())?;
                    let transformed =
                        match frozen_trend_type {
                            0 => {
                                // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
                                let parameters = unsafe {
                                    input_slice(
                                        frozen_trend_parameters,
                                        2,
                                        "frozen parametric trend parameters",
                                    )?
                                };
                                vst_parametric(
                                    &normalized,
                                    ParametricDispersionTrend {
                                        extra_pois: parameters[0],
                                        asympt_disp: parameters[1],
                                    },
                                )
                            }
                            1 => {
                                // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
                                let curve = unsafe {
                                    input_slice(frozen_local_curve, 1998, "frozen local VST curve")?
                                };
                                let inverse_size_factor_mean = local_vst_inverse_factor_mean(
                                    match &normalization_factors {
                                        Some(factors) => Offsets::PerCell(factors),
                                        None => Offsets::PerSample(&size_factors),
                                    },
                                )
                                .map_err(|error| error.to_string())?;
                                vst_local_from_dispersion_grid(
                                    &normalized,
                                    &curve[..999],
                                    &curve[999..],
                                    inverse_size_factor_mean,
                                )
                            }
                            2 => {
                                // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
                                let parameters = unsafe {
                                    input_slice(
                                        frozen_trend_parameters,
                                        2,
                                        "frozen mean trend parameters",
                                    )?
                                };
                                vst_mean(&normalized, parameters[0])
                            }
                            _ => return Err(
                                "frozen dispersion trend type must be parametric, local, or mean"
                                    .into(),
                            ),
                        }
                        .map_err(|error| error.to_string())?;
                    (transformed, size_factors)
                } else {
                    let fit = trend_stage_fit(&builder, &counts, &design)?;
                    let transformed = fit.vst(&counts).map_err(|error| error.to_string())?;
                    fit_type_used[0] = fit_type_code(fit.dispersion_trend.as_ref(), fit_type);
                    (transformed, fit.size_factors)
                };
                // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
                unsafe {
                    output_slice(size_factors_output, n_samples, "size factors")?
                        .copy_from_slice(&output_size_factors);
                    output_slice(intercept_output, n_genes, "rlog intercepts")?.fill(f64::NAN);
                    output_slice(prior_variance_output, 1, "rlog prior variance")?[0] = f64::NAN;
                }
                transformed
            }
            _ => return Err("unknown transform type".into()),
        };
        // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
        unsafe { copy_matrix_to_r(&transformed, transformed_output) }
    })
}

/// R's `fitType` code (0 parametric, 1 local, 2 mean) of the trend a fit used,
/// which differs from the request when a failed parametric fit fell back to
/// the local trend; the request when no trend was fitted.
pub(crate) fn fit_type_code(trend: Option<&DispersionTrendFit>, requested: i32) -> i32 {
    match trend {
        Some(DispersionTrendFit::Parametric(_)) => 0,
        Some(DispersionTrendFit::Local(_)) => 1,
        Some(DispersionTrendFit::Mean(_)) => 2,
        None => requested,
    }
}

/// Run the pipeline through the dispersion trend with GLM fitted means.
pub(crate) fn trend_stage_fit(
    builder: &DeseqBuilder,
    counts: &CountMatrix,
    design: &DesignMatrix,
) -> Result<DeseqFit, String> {
    let request = Analysis::new(counts, design)
        .means(MeanModel::GlmMu)
        .stop_after(Stage::Trend);
    match builder.run(&request).map_err(|error| error.to_string())? {
        RunOutput::Fit(fit) => Ok(fit),
        _ => Err("the trend-stage run returned no fit".into()),
    }
}
