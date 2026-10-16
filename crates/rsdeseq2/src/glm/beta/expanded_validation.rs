// Expanded-fit options, validation and covariance handling.

fn options_with_beta_prior_variance(
    design: &DesignMatrix,
    beta_prior_variance: &[f64],
    options: IrlsOptions,
) -> Result<IrlsOptions, DeseqError> {
    if beta_prior_variance.len() != design.n_coefficients() {
        return Err(invalid_dimensions(
            "beta prior variance coefficients",
            design.n_coefficients(),
            beta_prior_variance.len(),
        ));
    }
    let ridge_lambda = beta_prior_variance_to_ridge_lambda(beta_prior_variance)?;
    Ok(options.ridge_lambda_by_coefficient(ridge_lambda))
}

fn validate_expanded_coefficient_groups(
    n_coefficients: usize,
    coefficient_groups: &[Vec<usize>],
    context: &str,
) -> Result<(), DeseqError> {
    if coefficient_groups.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: context.to_string(),
            expected: 1,
            actual: 0,
        });
    }
    for (group_idx, group) in coefficient_groups.iter().enumerate() {
        validate_expanded_contrast_group(
            n_coefficients,
            group,
            &format!("{context} group {group_idx}"),
        )?;
    }
    Ok(())
}

fn validate_expanded_covariance_inputs(
    expanded_covariance: &RowMajorMatrix<f64>,
    n_expanded_coefficients: usize,
    coefficient_groups: &[Vec<usize>],
) -> Result<(), DeseqError> {
    if n_expanded_coefficients == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "expanded covariance coefficient count".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    let expected = n_expanded_coefficients
        .checked_mul(n_expanded_coefficients)
        .ok_or_else(|| DeseqError::InvalidDimensions {
            context: "expanded covariance columns".to_string(),
            expected: usize::MAX,
            actual: expanded_covariance.n_cols(),
        })?;
    if expanded_covariance.n_cols() != expected {
        return Err(invalid_dimensions(
            "expanded covariance columns",
            expected,
            expanded_covariance.n_cols(),
        ));
    }
    validate_expanded_coefficient_groups(
        n_expanded_coefficients,
        coefficient_groups,
        "expanded covariance coefficient groups",
    )
}

fn validate_expanded_fit_collapse_inputs(
    expanded_fit: &NbinomGlmFit,
    standard_design: &DesignMatrix,
    coefficient_groups: &[Vec<usize>],
) -> Result<(), DeseqError> {
    validate_expanded_coefficient_groups(
        expanded_fit.beta.n_cols(),
        coefficient_groups,
        "expanded fit coefficient groups",
    )?;
    if standard_design.n_coefficients() != coefficient_groups.len() {
        return Err(invalid_dimensions(
            "standard design coefficients",
            coefficient_groups.len(),
            standard_design.n_coefficients(),
        ));
    }
    if standard_design.n_samples() != expanded_fit.mu.n_cols() {
        return Err(invalid_dimensions(
            "standard design samples",
            expanded_fit.mu.n_cols(),
            standard_design.n_samples(),
        ));
    }
    validate_expanded_fit_row_vector(
        "expanded fit log-likelihood rows",
        expanded_fit.beta.n_rows(),
        expanded_fit.log_like.len(),
    )?;
    validate_expanded_fit_row_vector(
        "expanded fit convergence rows",
        expanded_fit.beta.n_rows(),
        expanded_fit.beta_converged.len(),
    )?;
    validate_expanded_fit_row_vector(
        "expanded fit beta-iteration rows",
        expanded_fit.beta.n_rows(),
        expanded_fit.beta_iter.len(),
    )?;
    if expanded_fit.beta_se.n_rows() != expanded_fit.beta.n_rows()
        || expanded_fit.beta_se.n_cols() != expanded_fit.beta.n_cols()
    {
        return Err(DeseqError::InvalidDimensions {
            context: "expanded fit beta SE matrix".to_string(),
            expected: expanded_fit.beta.len(),
            actual: expanded_fit.beta_se.len(),
        });
    }
    if expanded_fit.mu.n_rows() != expanded_fit.beta.n_rows() {
        return Err(invalid_dimensions(
            "expanded fit fitted mean rows",
            expanded_fit.beta.n_rows(),
            expanded_fit.mu.n_rows(),
        ));
    }
    if expanded_fit.hat_diagonal.n_rows() != expanded_fit.beta.n_rows()
        || expanded_fit.hat_diagonal.n_cols() != expanded_fit.mu.n_cols()
    {
        return Err(DeseqError::InvalidDimensions {
            context: "expanded fit hat diagonal matrix".to_string(),
            expected: expanded_fit.beta.n_rows() * expanded_fit.mu.n_cols(),
            actual: expanded_fit.hat_diagonal.len(),
        });
    }
    Ok(())
}

fn validate_expanded_fit_row_vector(
    context: &str,
    expected: usize,
    actual: usize,
) -> Result<(), DeseqError> {
    if actual != expected {
        return Err(invalid_dimensions(context, expected, actual));
    }
    Ok(())
}

fn validate_expanded_contrast_group(
    n_coefficients: usize,
    columns: &[usize],
    context: &str,
) -> Result<(), DeseqError> {
    if columns.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: context.to_string(),
            expected: 1,
            actual: 0,
        });
    }
    let mut seen = vec![false; n_coefficients];
    for &column in columns {
        if column >= n_coefficients {
            return Err(DeseqError::InvalidDimensions {
                context: context.to_string(),
                expected: n_coefficients.saturating_sub(1),
                actual: column,
            });
        }
        if seen[column] {
            return Err(DeseqError::InvalidOptions {
                reason: format!("{context} contains duplicate column {column}"),
            });
        }
        seen[column] = true;
    }
    Ok(())
}

fn average_covariance_block(
    covariance: &[f64],
    n_coefficients: usize,
    left_group: &[usize],
    right_group: &[usize],
    gene: usize,
) -> Result<f64, DeseqError> {
    let terms = left_group
        .iter()
        .flat_map(|&left| {
            right_group
                .iter()
                .map(move |&right| covariance[left * n_coefficients + right])
        })
        .collect::<Vec<_>>();
    let sum = checked_scaled_sum(terms).ok_or_else(|| DeseqError::NonFiniteValue {
        context: "expanded beta covariance average".to_string(),
        index: Some(gene),
        value: f64::NAN,
    })?;
    let denominator = left_group
        .len()
        .checked_mul(right_group.len())
        .ok_or_else(|| DeseqError::InvalidDimensions {
            context: "expanded covariance average denominator".to_string(),
            expected: usize::MAX,
            actual: 0,
        })?;
    checked_div2(sum, denominator as f64).ok_or_else(|| DeseqError::NonFiniteValue {
        context: "expanded beta covariance average".to_string(),
        index: Some(gene),
        value: f64::NAN,
    })
}

fn beta_se_from_covariance(
    covariance: &RowMajorMatrix<f64>,
    n_coefficients: usize,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let expected_cols = n_coefficients.checked_mul(n_coefficients).ok_or_else(|| {
        DeseqError::InvalidDimensions {
            context: "collapsed covariance columns".to_string(),
            expected: usize::MAX,
            actual: covariance.n_cols(),
        }
    })?;
    if covariance.n_cols() != expected_cols {
        return Err(invalid_dimensions(
            "collapsed covariance columns",
            expected_cols,
            covariance.n_cols(),
        ));
    }
    let mut values = Vec::with_capacity(covariance.n_rows() * n_coefficients);
    for gene in 0..covariance.n_rows() {
        let row = covariance.row(gene)?;
        for coefficient in 0..n_coefficients {
            let variance = row[coefficient * n_coefficients + coefficient];
            if !variance.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: "collapsed beta standard error".to_string(),
                    index: Some(gene),
                    value: variance,
                });
            }
            values.push(variance.max(0.0).sqrt());
        }
    }
    RowMajorMatrix::from_row_major(covariance.n_rows(), n_coefficients, values)
}
