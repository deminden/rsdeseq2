fn best_log_alpha(
    objective: DispersionObjectiveInput<'_>,
    grid: &[f64],
) -> Result<(f64, f64), DeseqError> {
    let mut best_log = grid[0];
    let mut best_score = dispersion_log_posterior_objective(objective, best_log)?;
    for log_alpha in grid.iter().copied().skip(1) {
        let score = dispersion_log_posterior_objective(objective, log_alpha)?;
        if score > best_score {
            best_log = log_alpha;
            best_score = score;
        }
    }
    Ok((best_log, best_score))
}

fn linspace(start: f64, end: f64, len: usize) -> Result<Vec<f64>, DeseqError> {
    if len == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion grid points".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    if !start.is_finite() || !end.is_finite() {
        return Err(DeseqError::NonFiniteValue {
            context: "dispersion grid endpoint".to_string(),
            index: None,
            value: if start.is_finite() { end } else { start },
        });
    }
    if len == 1 {
        return Ok(vec![start]);
    }
    let span = checked_sub(end, start, 0, "dispersion grid span")?;
    let step = span / (len as f64 - 1.0);
    if !step.is_finite() {
        return Err(DeseqError::NonFiniteValue {
            context: "dispersion grid step".to_string(),
            index: None,
            value: step,
        });
    }
    (0..len)
        .map(|idx| {
            let offset = checked_mul(step, idx as f64, idx, "dispersion grid offset")?;
            checked_add(start, offset, idx, "dispersion grid value")
        })
        .collect()
}

fn compact_counts_rows(
    counts: &CountMatrix,
    row_indices: &[usize],
) -> Result<CountMatrix, DeseqError> {
    let mut values = Vec::with_capacity(row_indices.len() * counts.n_samples());
    for row in row_indices {
        values.extend_from_slice(counts.row(*row)?);
    }
    let gene_names = counts.gene_names().map(|names| {
        row_indices
            .iter()
            .map(|row| names[*row].clone())
            .collect::<Vec<_>>()
    });
    let sample_names = counts.sample_names().map(<[String]>::to_vec);
    CountMatrix::from_row_major_u32_with_names(
        row_indices.len(),
        counts.n_samples(),
        values,
        gene_names,
        sample_names,
    )
}

fn compact_matrix_rows(
    matrix: &RowMajorMatrix<f64>,
    row_indices: &[usize],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let mut values = Vec::with_capacity(row_indices.len() * matrix.n_cols());
    for row in row_indices {
        values.extend_from_slice(matrix.row(*row)?);
    }
    RowMajorMatrix::from_row_major(row_indices.len(), matrix.n_cols(), values)
}

fn compact_gene_values(values: &[f64], row_indices: &[usize]) -> Result<Vec<f64>, DeseqError> {
    let mut compact = Vec::with_capacity(row_indices.len());
    for row in row_indices {
        let Some(value) = values.get(*row).copied() else {
            return Err(invalid_dimensions(
                "compact gene values",
                row + 1,
                values.len(),
            ));
        };
        compact.push(value);
    }
    Ok(compact)
}

fn validate_gene_est_inputs(
    input: &GeneWiseDispersionInput<'_>,
    options: GeneWiseDispersionOptions,
) -> Result<(), DeseqError> {
    if input.design.n_samples() != input.counts.n_samples() {
        return Err(invalid_dimensions(
            "dispersion design samples",
            input.counts.n_samples(),
            input.design.n_samples(),
        ));
    }
    if input.normalized_counts.n_rows() != input.counts.n_genes()
        || input.normalized_counts.n_cols() != input.counts.n_samples()
    {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion normalized counts".to_string(),
            expected: input.counts.n_genes() * input.counts.n_samples(),
            actual: input.normalized_counts.len(),
        });
    }
    if input.base_mean.len() != input.counts.n_genes() {
        return Err(invalid_dimensions(
            "dispersion baseMean",
            input.counts.n_genes(),
            input.base_mean.len(),
        ));
    }
    if input.base_var.len() != input.counts.n_genes() {
        return Err(invalid_dimensions(
            "dispersion baseVar",
            input.counts.n_genes(),
            input.base_var.len(),
        ));
    }
    if input.all_zero.len() != input.counts.n_genes() {
        return Err(invalid_dimensions(
            "dispersion allZero",
            input.counts.n_genes(),
            input.all_zero.len(),
        ));
    }
    if input.size_factors.len() != input.counts.n_samples() {
        return Err(invalid_dimensions(
            "dispersion size factors",
            input.counts.n_samples(),
            input.size_factors.len(),
        ));
    }
    validate_size_factors(input.size_factors)?;
    if let Some(normalization_factors) = input.normalization_factors {
        validate_normalization_factors(
            normalization_factors,
            input.counts.n_genes(),
            input.counts.n_samples(),
        )?;
    }
    if let Some(observation_weights) = input.observation_weights {
        if observation_weights.n_rows() != input.counts.n_genes()
            || observation_weights.n_cols() != input.counts.n_samples()
        {
            return Err(DeseqError::InvalidDimensions {
                context: "dispersion observation weights".to_string(),
                expected: input.counts.n_genes() * input.counts.n_samples(),
                actual: observation_weights.len(),
            });
        }
        for (idx, weight) in observation_weights.as_slice().iter().copied().enumerate() {
            if !weight.is_finite() || weight < 0.0 {
                return Err(DeseqError::NonFiniteValue {
                    context: "dispersion observation weight".to_string(),
                    index: Some(idx),
                    value: weight,
                });
            }
        }
    }
    validate_gene_est_options(options)?;
    validate_dispersion_bounds(
        options.min_disp,
        max_dispersion(options, input.counts.n_samples()),
    )?;
    Ok(())
}

fn validate_gene_est_options(options: GeneWiseDispersionOptions) -> Result<(), DeseqError> {
    if !options.min_disp.is_finite() || options.min_disp <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "min_disp must be finite and positive".to_string(),
        });
    }
    if (options.min_disp / 10.0).ln() <= -30.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "log(min_disp / 10) must be above -30 for numerical stability".to_string(),
        });
    }
    if !options.min_mu.is_finite() || options.min_mu <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "dispersion min_mu".to_string(),
            index: None,
            value: options.min_mu,
        });
    }
    if options.grid_points < 3 {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion grid points".to_string(),
            expected: 3,
            actual: options.grid_points,
        });
    }
    if !options.kappa_0.is_finite() || options.kappa_0 <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "kappa_0 must be finite and positive".to_string(),
        });
    }
    if !options.disp_tol.is_finite() || options.disp_tol <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "disp_tol must be finite and positive".to_string(),
        });
    }
    if options.maxit == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion maxit".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    if options.niter == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion niter".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    validate_weight_threshold(options.weight_threshold, "dispersion weight_threshold")?;
    Ok(())
}

fn validate_weight_threshold(value: f64, context: &str) -> Result<(), DeseqError> {
    if !value.is_finite() || value < 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("{context} must be finite and non-negative"),
        });
    }
    Ok(())
}

fn validate_dispersion_prior(prior: Option<DispersionPrior>) -> Result<(), DeseqError> {
    if let Some(prior) = prior {
        if !prior.log_mean.is_finite() {
            return Err(DeseqError::InvalidDispersion {
                reason: "dispersion prior log_mean must be finite".to_string(),
            });
        }
        if !prior.variance.is_finite() || prior.variance <= 0.0 {
            return Err(DeseqError::InvalidDispersion {
                reason: "dispersion prior variance must be finite and positive".to_string(),
            });
        }
    }
    Ok(())
}

fn validate_size_factors(size_factors: &[f64]) -> Result<(), DeseqError> {
    for (idx, value) in size_factors.iter().copied().enumerate() {
        if !value.is_finite() || value <= 0.0 {
            return Err(DeseqError::InvalidSizeFactors {
                reason: format!("size factor at sample {idx} must be finite and positive"),
            });
        }
    }
    Ok(())
}

fn validate_normalization_factors(
    normalization_factors: &RowMajorMatrix<f64>,
    n_genes: usize,
    n_samples: usize,
) -> Result<(), DeseqError> {
    if normalization_factors.n_rows() != n_genes || normalization_factors.n_cols() != n_samples {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion normalization factors".to_string(),
            expected: n_genes * n_samples,
            actual: normalization_factors.len(),
        });
    }
    for (idx, value) in normalization_factors.as_slice().iter().copied().enumerate() {
        validate_normalization_factor(value, idx)?;
    }
    Ok(())
}

fn validate_normalization_factor(value: f64, index: usize) -> Result<(), DeseqError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(DeseqError::InvalidSizeFactors {
            reason: format!("normalization factor at index {index} must be finite and positive"),
        });
    }
    Ok(())
}

fn validate_positive_mu(mu: f64, sample: usize) -> Result<(), DeseqError> {
    if !mu.is_finite() || mu <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "dispersion mean".to_string(),
            index: Some(sample),
            value: mu,
        });
    }
    Ok(())
}

fn checked_mul(left: f64, right: f64, index: usize, context: &str) -> Result<f64, DeseqError> {
    let value = left * right;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value,
        })
    }
}

fn checked_add(left: f64, right: f64, index: usize, context: &str) -> Result<f64, DeseqError> {
    let value = left + right;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value,
        })
    }
}

fn checked_sub(left: f64, right: f64, index: usize, context: &str) -> Result<f64, DeseqError> {
    let value = left - right;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value,
        })
    }
}

fn checked_div(left: f64, right: f64, index: usize, context: &str) -> Result<f64, DeseqError> {
    let value = left / right;
    if left.is_finite() && right.is_finite() && right != 0.0 && value.is_finite() {
        Ok(value)
    } else {
        Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value,
        })
    }
}

fn checked_matrix_add_assign(
    sum: &mut f64,
    term: f64,
    index: usize,
    context: &str,
) -> Result<(), DeseqError> {
    let value = *sum + term;
    if value.is_finite() {
        *sum = value;
        Ok(())
    } else {
        Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value,
        })
    }
}

fn checked_sum_indexed(
    values: impl IntoIterator<Item = f64>,
    context: &str,
) -> Result<f64, DeseqError> {
    let mut sum = 0.0;
    for (idx, value) in values.into_iter().enumerate() {
        checked_matrix_add_assign(&mut sum, value, idx, context)?;
    }
    Ok(sum)
}

#[derive(Clone, Copy, Debug)]
struct MuAlphaTerms {
    log1p: f64,
    ratio: f64,
    alpha_over_one_plus: f64,
    mu_squared_alpha_over_one_plus_squared: f64,
    inv_one_plus_squared: f64,
}

fn mu_alpha_terms(
    mu: f64,
    alpha: f64,
    index: usize,
    context: &str,
) -> Result<MuAlphaTerms, DeseqError> {
    let mu_alpha = mu * alpha;
    let terms = if mu_alpha.is_finite() {
        let one_plus = 1.0 + mu_alpha;
        let inv_one_plus = one_plus.recip();
        let ratio = mu_alpha * inv_one_plus;
        let alpha_over_one_plus = alpha * inv_one_plus;
        let inv_one_plus_squared = checked_mul(
            inv_one_plus,
            inv_one_plus,
            index,
            &format!("{context} inverse denominator square"),
        )?;
        let ratio_squared = checked_mul(
            ratio,
            ratio,
            index,
            &format!("{context} mean-dispersion ratio square"),
        )?;
        let mu_squared_alpha_over_one_plus_squared = checked_mul(
            ratio_squared,
            alpha.recip(),
            index,
            &format!("{context} mean curvature term"),
        )?;
        MuAlphaTerms {
            log1p: mu_alpha.ln_1p(),
            ratio,
            alpha_over_one_plus,
            mu_squared_alpha_over_one_plus_squared,
            inv_one_plus_squared,
        }
    } else {
        let log1p = mu.ln() + alpha.ln();
        let alpha_over_one_plus = mu.recip();
        let mu_squared_alpha_over_one_plus_squared = alpha.recip();
        MuAlphaTerms {
            log1p,
            ratio: 1.0,
            alpha_over_one_plus,
            mu_squared_alpha_over_one_plus_squared,
            inv_one_plus_squared: 0.0,
        }
    };
    for value in [
        terms.log1p,
        terms.ratio,
        terms.alpha_over_one_plus,
        terms.mu_squared_alpha_over_one_plus_squared,
        terms.inv_one_plus_squared,
    ] {
        if !value.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: context.to_string(),
                index: Some(index),
                value,
            });
        }
    }
    Ok(terms)
}

#[derive(Clone, Copy, Debug)]
struct CoxReidWeightTerms {
    weight: f64,
    d_weight: f64,
    d2_weight: f64,
}

fn cox_reid_weight_terms(
    mu: f64,
    alpha: f64,
    index: usize,
) -> Result<CoxReidWeightTerms, DeseqError> {
    let mu_alpha = mu_alpha_terms(mu, alpha, index, "Cox-Reid working weight")?;
    let weight = checked_mul(
        mu_alpha.ratio,
        alpha.recip(),
        index,
        "Cox-Reid working weight",
    )?;
    let weight_square = checked_mul(weight, weight, index, "Cox-Reid working weight square")?;
    let d_weight = -weight_square;
    let d2_weight = checked_mul(
        2.0,
        checked_mul(weight_square, weight, index, "Cox-Reid working weight cube")?,
        index,
        "Cox-Reid working second derivative weight",
    )?;
    Ok(CoxReidWeightTerms {
        weight,
        d_weight,
        d2_weight,
    })
}

fn checked_scaled_sum(values: &[f64], context: &str) -> Result<f64, DeseqError> {
    let mut scale = 0.0_f64;
    for value in values.iter().copied() {
        if !value.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: context.to_string(),
                index: None,
                value,
            });
        }
        scale = scale.max(value.abs());
    }
    if scale == 0.0 {
        return Ok(0.0);
    }
    let mut normalized_sum = 0.0;
    for value in values.iter().copied() {
        let term = value / scale;
        let next = normalized_sum + term;
        if !term.is_finite() || !next.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: context.to_string(),
                index: None,
                value: next,
            });
        }
        normalized_sum = next;
    }
    let sum = normalized_sum * scale;
    if sum.is_finite() {
        Ok(sum)
    } else {
        Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: None,
            value: sum,
        })
    }
}

fn checked_log_alpha_first_derivative(
    inv_alpha: f64,
    derivative_alpha: f64,
    context: &str,
) -> Result<f64, DeseqError> {
    checked_mul(inv_alpha, derivative_alpha, 0, context)
}

fn checked_log_alpha_second_derivative(
    second_alpha_sum: f64,
    inv_alpha: f64,
    first_alpha_sum: f64,
    first_log_alpha: f64,
    context: &str,
) -> Result<f64, DeseqError> {
    let alpha_first_scale = checked_mul(2.0, inv_alpha, 0, context)?;
    let alpha_first_term = -checked_mul(alpha_first_scale, first_alpha_sum, 0, context)?;
    checked_scaled_sum(
        &[second_alpha_sum, alpha_first_term, first_log_alpha],
        context,
    )
}

fn checked_cox_reid_log_alpha_second_derivative(
    second_alpha: f64,
    alpha: f64,
    first_log_alpha: f64,
) -> Result<f64, DeseqError> {
    let alpha_squared = checked_mul(
        alpha,
        alpha,
        0,
        "Cox-Reid log-alpha second derivative alpha square",
    )?;
    let alpha_term = checked_mul(
        second_alpha,
        alpha_squared,
        0,
        "Cox-Reid log-alpha second derivative alpha term",
    )?;
    checked_scaled_sum(
        &[alpha_term, first_log_alpha],
        "Cox-Reid log-alpha second derivative",
    )
}

fn validate_observation_weight_slice(
    weights: Option<&[f64]>,
    expected_len: usize,
    context: &str,
) -> Result<(), DeseqError> {
    let Some(weights) = weights else {
        return Ok(());
    };
    if weights.len() != expected_len {
        return Err(invalid_dimensions(context, expected_len, weights.len()));
    }
    for (idx, weight) in weights.iter().copied().enumerate() {
        if !weight.is_finite() || weight < 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: context.to_string(),
                index: Some(idx),
                value: weight,
            });
        }
    }
    Ok(())
}

fn validate_dispersion_bounds(min_disp: f64, max_disp: f64) -> Result<(), DeseqError> {
    if !max_disp.is_finite() || max_disp <= min_disp {
        return Err(DeseqError::InvalidDispersion {
            reason: "max dispersion must be finite and greater than min dispersion".to_string(),
        });
    }
    Ok(())
}

fn max_dispersion(options: GeneWiseDispersionOptions, n_samples: usize) -> f64 {
    options
        .max_disp
        .unwrap_or_else(|| 10.0_f64.max(n_samples as f64))
}
