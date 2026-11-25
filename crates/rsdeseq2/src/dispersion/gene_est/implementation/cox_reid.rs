use super::*;
use crate::errors::{ErrorKind, Requirement, invalid_parameter};
use crate::matrix::{try_alloc_matrix, try_zero_matrix};
use crate::numeric::math::ProfileMath;
// Fixed design products for the unweighted optimizer path. Products are
// checked at their original evaluation point, after the mean/weight checks.
pub(super) type PreparedCoxReidDesign = crate::design::SharedDesignProducts;

pub(super) fn prepare_unweighted_cox_reid_design(
    design: Option<&DesignMatrix>,
    weights: Option<&[f64]>,
    use_cox_reid: bool,
) -> Option<std::sync::Arc<PreparedCoxReidDesign>> {
    if weights.is_some() || !use_cox_reid {
        return None;
    }
    let design = design?;
    design.cached_cox_products(|| prepare_cox_reid_products(design))
}

pub(super) fn prepare_cox_reid_products(design: &DesignMatrix) -> Option<PreparedCoxReidDesign> {
    let p = design.n_coefficients();
    // Bound scratch memory; large designs keep the original streaming path.
    if p > 16 {
        return None;
    }
    let capacity = design.n_samples().checked_mul(p * (p + 1) / 2)?;
    if capacity > 1_000_000 {
        return None;
    }
    let mut offsets =
        try_alloc_matrix(design.n_samples() + 1, 1, "Cox-Reid design cache offsets").ok()?;
    let mut products = try_alloc_matrix(capacity, 1, "Cox-Reid design cache products").ok()?;
    let mut active = try_alloc_matrix(p, 1, "Cox-Reid design cache active columns").ok()?;
    for sample in 0..design.n_samples() {
        offsets.push(products.len());
        let row = design.matrix().row(sample).ok()?;
        active.clear();
        active.extend(
            row.iter()
                .copied()
                .enumerate()
                .filter(|(_, value)| *value != 0.0),
        );
        for (left_pos, (left, left_value)) in active.iter().copied().enumerate() {
            for (right, right_value) in active.iter().copied().skip(left_pos) {
                products.push((left, right, left_value * right_value));
            }
        }
    }
    offsets.push(products.len());
    Some(PreparedCoxReidDesign { offsets, products })
}

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

pub(super) fn cox_reid_adjustment_weighted_with_threshold(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
) -> Result<f64, DeseqError> {
    cox_reid_adjustment_weighted_with_threshold_prepared(
        design,
        mu,
        log_alpha,
        weights,
        weight_threshold,
        None,
    )
}

pub(super) fn cox_reid_adjustment_weighted_with_threshold_prepared(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
    prepared: Option<&PreparedCoxReidDesign>,
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
        return Err(invalid_parameter(
            ErrorKind::Dispersion,
            "log dispersion",
            Requirement::Finite,
        ));
    }
    let alpha = log_alpha.exp_p();
    if !alpha.is_finite() || alpha <= 0.0 {
        return Err(invalid_parameter(
            ErrorKind::Dispersion,
            "dispersion",
            Requirement::FiniteAndPositive,
        ));
    }
    let matrices = cox_reid_weighted_design_matrices_prepared(
        design,
        mu,
        alpha,
        weights,
        weight_threshold,
        CoxReidMatrixOrder::Value,
        prepared,
    )?;
    let determinant = determinant(&matrices.xtwx)?;
    if !determinant.is_finite() || determinant <= 0.0 {
        return Err(DeseqError::InvalidDimensions {
            context: "Cox-Reid weighted design determinant".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        });
    }
    Ok(-0.5 * determinant.ln_p())
}

pub(super) struct CoxReidDesignMatrices {
    xtwx: DenseMatrix,
    d_xtwx: Option<DenseMatrix>,
    d2_xtwx: Option<DenseMatrix>,
}

impl CoxReidDesignMatrices {
    fn new<const DERIVATIVE: bool, const SECOND_DERIVATIVE: bool>(
        columns: usize,
    ) -> Result<Self, DeseqError> {
        Ok(Self {
            xtwx: try_zero_matrix(columns, columns, "Cox-Reid crossproduct")?,
            d_xtwx: DERIVATIVE
                .then(|| try_zero_matrix(columns, columns, "Cox-Reid derivative"))
                .transpose()?,
            d2_xtwx: SECOND_DERIVATIVE
                .then(|| try_zero_matrix(columns, columns, "Cox-Reid second derivative"))
                .transpose()?,
        })
    }
}

#[derive(Clone, Copy)]
pub(super) enum CoxReidMatrixOrder {
    Value,
    Derivative,
    SecondDerivative,
}

pub(super) fn cox_reid_weighted_design_matrices_prepared(
    design: &DesignMatrix,
    mu: &[f64],
    alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
    order: CoxReidMatrixOrder,
    prepared: Option<&PreparedCoxReidDesign>,
) -> Result<CoxReidDesignMatrices, DeseqError> {
    match order {
        CoxReidMatrixOrder::Value => cox_reid_design_matrices_impl::<false, false>(
            design,
            mu,
            alpha,
            weights,
            weight_threshold,
            prepared,
        ),
        CoxReidMatrixOrder::Derivative => cox_reid_design_matrices_impl::<true, false>(
            design,
            mu,
            alpha,
            weights,
            weight_threshold,
            prepared,
        ),
        CoxReidMatrixOrder::SecondDerivative => cox_reid_design_matrices_impl::<true, true>(
            design,
            mu,
            alpha,
            weights,
            weight_threshold,
            prepared,
        ),
    }
}

pub(super) fn cox_reid_design_matrices_impl<
    const DERIVATIVE: bool,
    const SECOND_DERIVATIVE: bool,
>(
    design: &DesignMatrix,
    mu: &[f64],
    alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
    prepared: Option<&PreparedCoxReidDesign>,
) -> Result<CoxReidDesignMatrices, DeseqError> {
    if design.n_samples() != mu.len() {
        return Err(invalid_dimensions(
            "Cox-Reid design samples",
            mu.len(),
            design.n_samples(),
        ));
    }
    if !alpha.is_finite() || alpha <= 0.0 {
        return Err(invalid_parameter(
            ErrorKind::Dispersion,
            "dispersion",
            Requirement::FiniteAndPositive,
        ));
    }
    validate_observation_weight_slice(weights, mu.len(), "Cox-Reid weights")?;
    validate_weight_threshold(weight_threshold, "Cox-Reid weight threshold")?;
    if let Some(prepared) = prepared {
        return cox_reid_prepared_design_matrices::<DERIVATIVE, SECOND_DERIVATIVE>(
            design, mu, alpha, prepared,
        );
    }
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
    let mut matrices = CoxReidDesignMatrices::new::<DERIVATIVE, SECOND_DERIVATIVE>(p)?;
    let CoxReidDesignMatrices {
        xtwx,
        d_xtwx,
        d2_xtwx,
    } = &mut matrices;
    let mut active_columns = try_alloc_matrix(p, 1, "Cox-Reid active columns")?;
    for sample in selected_samples {
        let mu = mu[sample];
        validate_positive_mu(mu, sample)?;

        let row = design.matrix().row(sample)?;
        active_columns.clear();
        active_columns.extend(
            selected_columns
                .iter()
                .copied()
                .enumerate()
                .filter_map(|(idx, col)| (row[col] != 0.0).then_some((idx, row[col]))),
        );
        if active_columns.is_empty() {
            continue;
        }
        let weight_terms =
            cox_reid_weight_terms_for_order::<DERIVATIVE, SECOND_DERIVATIVE>(mu, alpha, sample)?;
        for (left_pos, (left_idx, left_value)) in active_columns.iter().copied().enumerate() {
            for (right_idx, right_value) in active_columns.iter().copied().skip(left_pos) {
                let x_product = checked::try_mul(left_value, right_value)
                    .nonfinite_at("Cox-Reid weighted design product", Some(sample))?;
                checked_matrix_add_assign(
                    &mut xtwx[(left_idx, right_idx)],
                    checked::try_mul(x_product, weight_terms.weight)
                        .nonfinite_at("Cox-Reid weighted design xtwx", Some(sample))?,
                    sample,
                    "Cox-Reid weighted design xtwx",
                )?;

                if let Some(d_xtwx) = d_xtwx {
                    checked_matrix_add_assign(
                        &mut d_xtwx[(left_idx, right_idx)],
                        checked::try_mul(x_product, weight_terms.d_weight)
                            .nonfinite_at("Cox-Reid weighted design derivative", Some(sample))?,
                        sample,
                        "Cox-Reid weighted design derivative",
                    )?;
                }
                if let Some(d2_xtwx) = d2_xtwx {
                    checked_matrix_add_assign(
                        &mut d2_xtwx[(left_idx, right_idx)],
                        checked::try_mul(x_product, weight_terms.d2_weight).nonfinite_at(
                            "Cox-Reid weighted design second derivative",
                            Some(sample),
                        )?,
                        sample,
                        "Cox-Reid weighted design second derivative",
                    )?;
                }
            }
        }
    }
    // Lower-triangle entries are read only after the sample reductions.
    for left in 0..p {
        for right in (left + 1)..p {
            xtwx[(right, left)] = xtwx[(left, right)];
            if let Some(matrix) = d_xtwx {
                matrix[(right, left)] = matrix[(left, right)];
            }
            if let Some(matrix) = d2_xtwx {
                matrix[(right, left)] = matrix[(left, right)];
            }
        }
    }
    Ok(matrices)
}

pub(super) fn cox_reid_prepared_design_matrices<
    const DERIVATIVE: bool,
    const SECOND_DERIVATIVE: bool,
>(
    design: &DesignMatrix,
    mu: &[f64],
    alpha: f64,
    prepared: &PreparedCoxReidDesign,
) -> Result<CoxReidDesignMatrices, DeseqError> {
    // The prepared unweighted subset is exactly all samples and columns.
    let p = design.n_coefficients();
    if mu.is_empty() || p == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "Cox-Reid weighted design subset".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        });
    }
    let mut matrices = CoxReidDesignMatrices::new::<DERIVATIVE, SECOND_DERIVATIVE>(p)?;
    let CoxReidDesignMatrices {
        xtwx,
        d_xtwx,
        d2_xtwx,
    } = &mut matrices;
    for (sample, &mu) in mu.iter().enumerate() {
        validate_positive_mu(mu, sample)?;
        let products = &prepared.products[prepared.offsets[sample]..prepared.offsets[sample + 1]];
        if products.is_empty() {
            continue;
        }
        let weight_terms =
            cox_reid_weight_terms_for_order::<DERIVATIVE, SECOND_DERIVATIVE>(mu, alpha, sample)?;
        for &(left_idx, right_idx, x_product) in products {
            if !x_product.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: "Cox-Reid weighted design product".to_string(),
                    index: Some(sample),
                    value: x_product,
                });
            }
            checked_matrix_add_assign(
                &mut xtwx[(left_idx, right_idx)],
                checked::try_mul(x_product, weight_terms.weight)
                    .nonfinite_at("Cox-Reid weighted design xtwx", Some(sample))?,
                sample,
                "Cox-Reid weighted design xtwx",
            )?;

            if let Some(d_xtwx) = d_xtwx {
                checked_matrix_add_assign(
                    &mut d_xtwx[(left_idx, right_idx)],
                    checked::try_mul(x_product, weight_terms.d_weight)
                        .nonfinite_at("Cox-Reid weighted design derivative", Some(sample))?,
                    sample,
                    "Cox-Reid weighted design derivative",
                )?;
            }
            if let Some(d2_xtwx) = d2_xtwx {
                checked_matrix_add_assign(
                    &mut d2_xtwx[(left_idx, right_idx)],
                    checked::try_mul(x_product, weight_terms.d2_weight)
                        .nonfinite_at("Cox-Reid weighted design second derivative", Some(sample))?,
                    sample,
                    "Cox-Reid weighted design second derivative",
                )?;
            }
        }
    }
    // Lower-triangle entries are read only after the sample reductions.
    for left in 0..p {
        for right in (left + 1)..p {
            xtwx[(right, left)] = xtwx[(left, right)];
            if let Some(matrix) = d_xtwx {
                matrix[(right, left)] = matrix[(left, right)];
            }
            if let Some(matrix) = d2_xtwx {
                matrix[(right, left)] = matrix[(left, right)];
            }
        }
    }
    Ok(matrices)
}

pub(super) fn cox_reid_sample_indices(
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

pub(super) fn cox_reid_column_indices(
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

pub(super) fn trace_product(left: &DenseMatrix, right: &DenseMatrix) -> Result<f64, DeseqError> {
    let product = matmul(left, right)?;
    checked_sum_indexed(
        (0..product.nrows().min(product.ncols())).map(|i| product[(i, i)]),
        "Cox-Reid trace product",
    )
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

pub(super) fn cox_reid_adjustment_derivative_weighted_with_threshold(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
) -> Result<f64, DeseqError> {
    cox_reid_adjustment_derivative_weighted_with_threshold_prepared(
        design,
        mu,
        log_alpha,
        weights,
        weight_threshold,
        None,
    )
}

pub(super) fn cox_reid_adjustment_derivative_weighted_with_threshold_prepared(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
    prepared: Option<&PreparedCoxReidDesign>,
) -> Result<f64, DeseqError> {
    let alpha = log_alpha.exp_p();
    let matrices = cox_reid_weighted_design_matrices_prepared(
        design,
        mu,
        alpha,
        weights,
        weight_threshold,
        CoxReidMatrixOrder::Derivative,
        prepared,
    )?;
    let mut inverse = matrices.xtwx;
    if !try_inverse_in_place(&mut inverse)? {
        return Err(DeseqError::InvalidDimensions {
            context: "Cox-Reid weighted design inverse".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        });
    }
    let d_xtwx = matrices
        .d_xtwx
        .ok_or_else(|| DeseqError::InvalidDimensions {
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

pub(super) fn cox_reid_adjustment_second_derivative_weighted_with_threshold(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
) -> Result<f64, DeseqError> {
    cox_reid_adjustment_second_derivative_weighted_with_threshold_prepared(
        design,
        mu,
        log_alpha,
        weights,
        weight_threshold,
        None,
    )
}

pub(super) fn cox_reid_adjustment_second_derivative_weighted_with_threshold_prepared(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
    prepared: Option<&PreparedCoxReidDesign>,
) -> Result<f64, DeseqError> {
    let alpha = log_alpha.exp_p();
    let matrices = cox_reid_weighted_design_matrices_prepared(
        design,
        mu,
        alpha,
        weights,
        weight_threshold,
        CoxReidMatrixOrder::SecondDerivative,
        prepared,
    )?;
    let mut inverse = matrices.xtwx;
    if !try_inverse_in_place(&mut inverse)? {
        return Err(DeseqError::InvalidDimensions {
            context: "Cox-Reid weighted design inverse".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        });
    }
    let d_xtwx = matrices
        .d_xtwx
        .ok_or_else(|| DeseqError::InvalidDimensions {
            context: "Cox-Reid weighted design derivative".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        })?;
    let d2_xtwx = matrices
        .d2_xtwx
        .ok_or_else(|| DeseqError::InvalidDimensions {
            context: "Cox-Reid weighted design second derivative".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        })?;
    let second_trace_product = matmul(&matmul(&matmul(&inverse, &d_xtwx)?, &inverse)?, &d_xtwx)?;
    let second_trace = checked_sum_indexed(
        (0..second_trace_product
            .nrows()
            .min(second_trace_product.ncols()))
            .map(|i| second_trace_product[(i, i)]),
        "Cox-Reid second trace product",
    )?;
    let trace_bi_d2b = trace_product(&inverse, &d2_xtwx)?;
    let second_alpha = 0.5 * (second_trace - trace_bi_d2b);
    let first_log_alpha = -0.5 * trace_product(&inverse, &d_xtwx)? * alpha;
    checked_cox_reid_log_alpha_second_derivative(second_alpha, alpha, first_log_alpha)
}

#[inline(always)]
pub(super) fn cox_reid_weight_terms_for_order<
    const DERIVATIVE: bool,
    const SECOND_DERIVATIVE: bool,
>(
    mu: f64,
    alpha: f64,
    index: usize,
) -> Result<CoxReidWeightTerms, DeseqError> {
    // Both matrix loops validate mu first; their prefix validates alpha.
    // For alpha >= 1e-100, 1/alpha < 2e100 and the positive ratio < 2.
    // The omitted curvature is < 8e100, weight square < 16e200,
    // and twice the weight cube < 128e300: all are below f64::MAX.
    // Inverse denominator squared is <= 1 and alpha/denominator <= alpha.
    // Thus every omitted check would succeed. Extreme dispersions retain the
    // complete original checked helper, including its error precedence.
    if alpha < 1e-100 {
        return cox_reid_weight_terms(mu, alpha, index);
    }
    let product = mu * alpha;
    let ratio = if product.is_finite() {
        let one_plus = 1.0 + product;
        let inverse = one_plus.recip();
        product * inverse
    } else {
        1.0
    };
    let weight = checked::try_mul(ratio, alpha.recip())
        .nonfinite_at("Cox-Reid working weight", Some(index))?;
    let square = if DERIVATIVE || SECOND_DERIVATIVE {
        checked::try_mul(weight, weight)
            .nonfinite_at("Cox-Reid working weight square", Some(index))?
    } else {
        0.0
    };
    let d2_weight = if SECOND_DERIVATIVE {
        checked::try_mul(
            2.0,
            checked::try_mul(square, weight)
                .nonfinite_at("Cox-Reid working weight cube", Some(index))?,
        )
        .nonfinite_at("Cox-Reid working second derivative weight", Some(index))?
    } else {
        0.0
    };
    Ok(CoxReidWeightTerms {
        weight,
        d_weight: if DERIVATIVE { -square } else { 0.0 },
        d2_weight,
    })
}
