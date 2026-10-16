// Expanded-model coefficient, covariance and contrast collapse.

/// Average expanded-model coefficients into standard model-term coefficients.
///
/// DESeq2's beta-prior workflow can fit an expanded model matrix and then
/// collapse groups of expanded columns back to the reported standard
/// coefficients. This primitive performs that numeric collapse for an existing
/// log2-scale beta matrix: each output column is the per-gene average of the
/// input columns listed in the corresponding group.
pub fn average_expanded_model_coefficients(
    expanded_beta: &RowMajorMatrix<f64>,
    coefficient_groups: &[Vec<usize>],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    validate_expanded_coefficient_groups(
        expanded_beta.n_cols(),
        coefficient_groups,
        "expanded beta coefficient groups",
    )?;

    let mut values = Vec::with_capacity(expanded_beta.n_rows() * coefficient_groups.len());
    for gene in 0..expanded_beta.n_rows() {
        let row = expanded_beta.row(gene)?;
        for group in coefficient_groups {
            let grouped = group.iter().map(|&column| row[column]).collect::<Vec<_>>();
            values.push(checked_scaled_mean(&grouped).ok_or_else(|| {
                DeseqError::NonFiniteValue {
                    context: "expanded beta coefficient average".to_string(),
                    index: Some(gene),
                    value: f64::NAN,
                }
            })?);
        }
    }

    RowMajorMatrix::from_row_major(expanded_beta.n_rows(), coefficient_groups.len(), values)
}

/// Propagate expanded-model covariance through grouped coefficient averaging.
///
/// The input covariance matrix stores each gene's expanded-model coefficient
/// covariance as one row-major `p x p` row. The output stores each gene's
/// collapsed `q x q` covariance matrix, where `q = coefficient_groups.len()`
/// and every output cell is `A Sigma A'` for the averaging matrix implied by
/// `coefficient_groups`.
pub fn average_expanded_model_covariances(
    expanded_covariance: &RowMajorMatrix<f64>,
    n_expanded_coefficients: usize,
    coefficient_groups: &[Vec<usize>],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    validate_expanded_covariance_inputs(
        expanded_covariance,
        n_expanded_coefficients,
        coefficient_groups,
    )?;

    let n_groups = coefficient_groups.len();
    let mut values = Vec::with_capacity(expanded_covariance.n_rows() * n_groups * n_groups);
    for gene in 0..expanded_covariance.n_rows() {
        let covariance = expanded_covariance.row(gene)?;
        for left_group in coefficient_groups {
            for right_group in coefficient_groups {
                values.push(average_covariance_block(
                    covariance,
                    n_expanded_coefficients,
                    left_group,
                    right_group,
                    gene,
                )?);
            }
        }
    }

    RowMajorMatrix::from_row_major(expanded_covariance.n_rows(), n_groups * n_groups, values)
}

/// Collapse an expanded-model GLM fit to standard-design coefficients.
///
/// This keeps gene/sample diagnostics and fitted means from the expanded fit,
/// replaces the reported model matrix with `standard_design`, averages beta
/// coefficients by group, propagates covariance through the same averaging
/// matrix, and recomputes standard errors from the collapsed covariance
/// diagonal. It is a primitive building block for DESeq2's beta-prior expanded
/// model workflow; callers remain responsible for constructing the expanded
/// and standard design matrices from formula metadata.
pub fn collapse_expanded_model_fit(
    expanded_fit: &NbinomGlmFit,
    standard_design: &DesignMatrix,
    coefficient_groups: &[Vec<usize>],
) -> Result<NbinomGlmFit, DeseqError> {
    validate_expanded_fit_collapse_inputs(expanded_fit, standard_design, coefficient_groups)?;
    let covariance =
        expanded_fit
            .beta_covariance
            .as_ref()
            .ok_or_else(|| DeseqError::UnsupportedFeature {
                feature: "expanded model fit collapse requires beta covariance matrices"
                    .to_string(),
            })?;
    let beta = average_expanded_model_coefficients(&expanded_fit.beta, coefficient_groups)?;
    let beta_covariance = average_expanded_model_covariances(
        covariance,
        expanded_fit.beta.n_cols(),
        coefficient_groups,
    )?;
    let beta_se = beta_se_from_covariance(&beta_covariance, beta.n_cols())?;

    Ok(NbinomGlmFit {
        log_like: expanded_fit.log_like.clone(),
        beta_converged: expanded_fit.beta_converged.clone(),
        beta,
        beta_se,
        beta_optim_start: expanded_fit.beta_optim_start.clone(),
        beta_covariance: Some(beta_covariance),
        mu: expanded_fit.mu.clone(),
        beta_iter: expanded_fit.beta_iter.clone(),
        beta_optim_iter: expanded_fit.beta_optim_iter.clone(),
        beta_optim_start_objective: expanded_fit.beta_optim_start_objective.clone(),
        beta_optim_objective: expanded_fit.beta_optim_objective.clone(),
        beta_optim_gradient_norm: expanded_fit.beta_optim_gradient_norm.clone(),
        model_matrix: standard_design.clone(),
        n_terms: standard_design.n_coefficients(),
        hat_diagonal: expanded_fit.hat_diagonal.clone(),
    })
}

/// Build a numeric contrast vector from two expanded-model coefficient groups.
///
/// The returned vector has `+1 / n` weights across numerator columns and
/// `-1 / m` weights across denominator columns, matching the averaging shape
/// needed to compare collapsed expanded-model terms.
pub fn expanded_model_group_contrast(
    n_coefficients: usize,
    numerator_columns: &[usize],
    denominator_columns: &[usize],
) -> Result<Vec<f64>, DeseqError> {
    validate_expanded_contrast_group(
        n_coefficients,
        numerator_columns,
        "expanded contrast numerator columns",
    )?;
    validate_expanded_contrast_group(
        n_coefficients,
        denominator_columns,
        "expanded contrast denominator columns",
    )?;
    for &column in numerator_columns {
        if denominator_columns.contains(&column) {
            return Err(DeseqError::InvalidOptions {
                reason: "expanded contrast numerator and denominator columns must be disjoint"
                    .to_string(),
            });
        }
    }

    let numerator_weight = checked_div2(1.0, numerator_columns.len() as f64).ok_or_else(|| {
        DeseqError::InvalidOptions {
            reason: "expanded contrast numerator weight is non-finite".to_string(),
        }
    })?;
    let denominator_weight =
        checked_div2(-1.0, denominator_columns.len() as f64).ok_or_else(|| {
            DeseqError::InvalidOptions {
                reason: "expanded contrast denominator weight is non-finite".to_string(),
            }
        })?;

    let mut contrast = vec![0.0; n_coefficients];
    for &column in numerator_columns {
        contrast[column] = numerator_weight;
    }
    for &column in denominator_columns {
        contrast[column] = denominator_weight;
    }
    Ok(contrast)
}
