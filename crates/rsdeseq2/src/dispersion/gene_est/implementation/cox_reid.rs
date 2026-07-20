/// Cox-Reid adjustment term for one gene and design matrix.
pub fn cox_reid_adjustment(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    cox_reid_adjustment_weighted(design, mu, log_alpha, None)
}

/// Cox-Reid adjustment term with optional DESeq2-style weighted sample subset.
pub fn cox_reid_adjustment_weighted(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    cox_reid_adjustment_weighted_with_threshold(
        design,
        mu,
        log_alpha,
        weights,
        GeneWiseDispersionOptions::default().weight_threshold,
    )
}

fn cox_reid_adjustment_weighted_with_threshold(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
) -> Result<f64, DeseqError> {
    if design.n_samples() != mu.len() {
        return Err(invalid_dimensions(
            "Cox-Reid design samples",
            mu.len(),
            design.n_samples(),
        ));
    }
    validate_observation_weight_slice(weights, mu.len(), "Cox-Reid weights")?;
    if !log_alpha.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "log dispersion must be finite".to_string(),
        });
    }
    let alpha = log_alpha.exp();
    if !alpha.is_finite() || alpha <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "dispersion must be finite and positive".to_string(),
        });
    }
    let matrices = cox_reid_weighted_design_matrices_with_threshold(
        design,
        mu,
        alpha,
        weights,
        weight_threshold,
        CoxReidMatrixOrder::Value,
    )?;
    let determinant = matrices.xtwx.determinant();
    if !determinant.is_finite() || determinant <= 0.0 {
        return Err(DeseqError::InvalidDimensions {
            context: "Cox-Reid weighted design determinant".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        });
    }
    Ok(-0.5 * determinant.ln())
}

struct CoxReidDesignMatrices {
    xtwx: DMatrix<f64>,
    d_xtwx: Option<DMatrix<f64>>,
    d2_xtwx: Option<DMatrix<f64>>,
}

#[derive(Clone, Copy)]
enum CoxReidMatrixOrder {
    Value,
    Derivative,
    SecondDerivative,
}

impl CoxReidMatrixOrder {
    fn needs_derivative(self) -> bool {
        matches!(self, Self::Derivative | Self::SecondDerivative)
    }

    fn needs_second_derivative(self) -> bool {
        matches!(self, Self::SecondDerivative)
    }
}

fn cox_reid_weighted_design_matrices_with_threshold(
    design: &DesignMatrix,
    mu: &[f64],
    alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
    order: CoxReidMatrixOrder,
) -> Result<CoxReidDesignMatrices, DeseqError> {
    if design.n_samples() != mu.len() {
        return Err(invalid_dimensions(
            "Cox-Reid design samples",
            mu.len(),
            design.n_samples(),
        ));
    }
    if !alpha.is_finite() || alpha <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "dispersion must be finite and positive".to_string(),
        });
    }
    validate_observation_weight_slice(weights, mu.len(), "Cox-Reid weights")?;
    validate_weight_threshold(weight_threshold, "Cox-Reid weight threshold")?;
    let selected_samples = cox_reid_sample_indices(weights, mu.len(), weight_threshold);
    let selected_columns = match weights {
        Some(_) => cox_reid_column_indices(design, &selected_samples)?,
        None => (0..design.n_coefficients()).collect(),
    };
    if selected_samples.is_empty() || selected_columns.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "Cox-Reid weighted design subset".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        });
    }
    // Weighted Cox-Reid drops samples below threshold, then drops columns that
    // are zero on the retained sample subset before building determinant terms.
    let p = selected_columns.len();
    let mut xtwx = DMatrix::<f64>::zeros(p, p);
    let mut d_xtwx = order
        .needs_derivative()
        .then(|| DMatrix::<f64>::zeros(p, p));
    let mut d2_xtwx = order
        .needs_second_derivative()
        .then(|| DMatrix::<f64>::zeros(p, p));
    for sample in selected_samples {
        let mu = mu[sample];
        validate_positive_mu(mu, sample)?;
        let row = design.matrix().row(sample)?;
        let active_columns = selected_columns
            .iter()
            .copied()
            .enumerate()
            .filter_map(|(idx, col)| (row[col] != 0.0).then_some((idx, row[col])))
            .collect::<Vec<_>>();
        if active_columns.is_empty() {
            continue;
        }
        let weight_terms = cox_reid_weight_terms(mu, alpha, sample)?;
        for (left_pos, (left_idx, left_value)) in active_columns.iter().copied().enumerate() {
            for (right_idx, right_value) in active_columns.iter().copied().skip(left_pos) {
                let x_product = checked_mul(
                    left_value,
                    right_value,
                    sample,
                    "Cox-Reid weighted design product",
                )?;
                checked_matrix_add_assign(
                    &mut xtwx[(left_idx, right_idx)],
                    checked_mul(
                        x_product,
                        weight_terms.weight,
                        sample,
                        "Cox-Reid weighted design xtwx",
                    )?,
                    sample,
                    "Cox-Reid weighted design xtwx",
                )?;
                if left_idx != right_idx {
                    xtwx[(right_idx, left_idx)] = xtwx[(left_idx, right_idx)];
                }
                if let Some(d_xtwx) = &mut d_xtwx {
                    checked_matrix_add_assign(
                        &mut d_xtwx[(left_idx, right_idx)],
                        checked_mul(
                            x_product,
                            weight_terms.d_weight,
                            sample,
                            "Cox-Reid weighted design derivative",
                        )?,
                        sample,
                        "Cox-Reid weighted design derivative",
                    )?;
                    if left_idx != right_idx {
                        d_xtwx[(right_idx, left_idx)] = d_xtwx[(left_idx, right_idx)];
                    }
                }
                if let Some(d2_xtwx) = &mut d2_xtwx {
                    checked_matrix_add_assign(
                        &mut d2_xtwx[(left_idx, right_idx)],
                        checked_mul(
                            x_product,
                            weight_terms.d2_weight,
                            sample,
                            "Cox-Reid weighted design second derivative",
                        )?,
                        sample,
                        "Cox-Reid weighted design second derivative",
                    )?;
                    if left_idx != right_idx {
                        d2_xtwx[(right_idx, left_idx)] = d2_xtwx[(left_idx, right_idx)];
                    }
                }
            }
        }
    }
    Ok(CoxReidDesignMatrices {
        xtwx,
        d_xtwx,
        d2_xtwx,
    })
}

fn cox_reid_sample_indices(
    weights: Option<&[f64]>,
    n_samples: usize,
    weight_threshold: f64,
) -> Vec<usize> {
    match weights {
        Some(weights) => weights
            .iter()
            .copied()
            .enumerate()
            .filter_map(|(sample, weight)| (weight > weight_threshold).then_some(sample))
            .collect(),
        None => (0..n_samples).collect(),
    }
}

fn cox_reid_column_indices(
    design: &DesignMatrix,
    selected_samples: &[usize],
) -> Result<Vec<usize>, DeseqError> {
    let mut selected = Vec::with_capacity(design.n_coefficients());
    for column in 0..design.n_coefficients() {
        let mut sum_abs = 0.0;
        for sample in selected_samples {
            checked_matrix_add_assign(
                &mut sum_abs,
                design.matrix().row(*sample)?[column].abs(),
                *sample,
                "Cox-Reid selected design column sum",
            )?;
        }
        if sum_abs > 0.0 {
            selected.push(column);
        }
    }
    Ok(selected)
}

fn trace_product(left: &DMatrix<f64>, right: &DMatrix<f64>) -> Result<f64, DeseqError> {
    let product = left * right;
    checked_sum_indexed(product.diagonal().iter().copied(), "Cox-Reid trace product")
}

/// Derivative of the Cox-Reid adjustment with respect to log alpha.
pub fn cox_reid_adjustment_derivative(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    cox_reid_adjustment_derivative_weighted(design, mu, log_alpha, None)
}

/// Derivative of the weighted Cox-Reid adjustment with respect to log alpha.
///
/// Observation weights define the DESeq2 `weightThreshold` sample subset for
/// the determinant; they do not multiply the Cox-Reid working weights.
pub fn cox_reid_adjustment_derivative_weighted(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    cox_reid_adjustment_derivative_weighted_with_threshold(
        design,
        mu,
        log_alpha,
        weights,
        GeneWiseDispersionOptions::default().weight_threshold,
    )
}

fn cox_reid_adjustment_derivative_weighted_with_threshold(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
) -> Result<f64, DeseqError> {
    let alpha = log_alpha.exp();
    let matrices = cox_reid_weighted_design_matrices_with_threshold(
        design,
        mu,
        alpha,
        weights,
        weight_threshold,
        CoxReidMatrixOrder::Derivative,
    )?;
    let Some(inverse) = matrices.xtwx.try_inverse() else {
        return Err(DeseqError::InvalidDimensions {
            context: "Cox-Reid weighted design inverse".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        });
    };
    let d_xtwx = matrices.d_xtwx.ok_or(DeseqError::InvalidDimensions {
        context: "Cox-Reid weighted design derivative".to_string(),
        expected: design.n_coefficients(),
        actual: 0,
    })?;
    Ok(-0.5 * trace_product(&inverse, &d_xtwx)? * alpha)
}

/// Second derivative of the Cox-Reid adjustment with respect to log alpha.
pub fn cox_reid_adjustment_second_derivative(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    cox_reid_adjustment_second_derivative_weighted(design, mu, log_alpha, None)
}

/// Second derivative of the weighted Cox-Reid adjustment with respect to log alpha.
///
/// Observation weights define the DESeq2 `weightThreshold` sample subset for
/// the determinant; they do not multiply the Cox-Reid working weights.
pub fn cox_reid_adjustment_second_derivative_weighted(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    cox_reid_adjustment_second_derivative_weighted_with_threshold(
        design,
        mu,
        log_alpha,
        weights,
        GeneWiseDispersionOptions::default().weight_threshold,
    )
}

fn cox_reid_adjustment_second_derivative_weighted_with_threshold(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
) -> Result<f64, DeseqError> {
    let alpha = log_alpha.exp();
    let matrices = cox_reid_weighted_design_matrices_with_threshold(
        design,
        mu,
        alpha,
        weights,
        weight_threshold,
        CoxReidMatrixOrder::SecondDerivative,
    )?;
    let Some(inverse) = matrices.xtwx.try_inverse() else {
        return Err(DeseqError::InvalidDimensions {
            context: "Cox-Reid weighted design inverse".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        });
    };
    let d_xtwx = matrices.d_xtwx.ok_or(DeseqError::InvalidDimensions {
        context: "Cox-Reid weighted design derivative".to_string(),
        expected: design.n_coefficients(),
        actual: 0,
    })?;
    let d2_xtwx = matrices.d2_xtwx.ok_or(DeseqError::InvalidDimensions {
        context: "Cox-Reid weighted design second derivative".to_string(),
        expected: design.n_coefficients(),
        actual: 0,
    })?;
    let second_trace_product = &inverse * &d_xtwx * &inverse * &d_xtwx;
    let second_trace = checked_sum_indexed(
        second_trace_product.diagonal().iter().copied(),
        "Cox-Reid second trace product",
    )?;
    let trace_bi_d2b = trace_product(&inverse, &d2_xtwx)?;
    let second_alpha = 0.5 * (second_trace - trace_bi_d2b);
    let first_log_alpha = -0.5 * trace_product(&inverse, &d_xtwx)? * alpha;
    checked_cox_reid_log_alpha_second_derivative(second_alpha, alpha, first_log_alpha)
}
