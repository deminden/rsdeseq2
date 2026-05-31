fn require_cli_fit<'a>(
    analysis: &'a CliAnalysisOutput,
    context: &str,
) -> Result<&'a DeseqFit, DeseqError> {
    analysis
        .fit
        .as_ref()
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!("{context} requires a native fit workflow"),
        })
}

fn require_cli_refit<'a>(
    analysis: &'a CliAnalysisOutput,
    context: &str,
) -> Result<&'a DeseqFit, DeseqError> {
    analysis
        .refit
        .as_ref()
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!("{context} requires rows to be refit"),
        })
}

fn fit_coefficient_names(fit: &DeseqFit) -> Option<&[String]> {
    fit.design
        .as_ref()
        .and_then(|design| design.coefficient_names())
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

fn read_cli_gene_numeric(
    path: impl Into<PathBuf>,
    counts: &crate::core::CountMatrix,
    context: &str,
) -> Result<Vec<f64>, DeseqError> {
    align_gene_numeric_values_to_genes(
        &read_labeled_gene_numeric_tsv(path.into(), context)?,
        counts
            .gene_names()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!("count gene names are required to align {context} values"),
            })?,
        context,
    )
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
