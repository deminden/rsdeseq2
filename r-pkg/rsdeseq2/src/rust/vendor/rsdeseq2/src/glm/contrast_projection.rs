//! Borrowed coefficient/covariance projection shared by every interface.
use super::inference_primitives::wald_from_estimate_rows;
use super::wald::{checked_product3_ordered_by_magnitude, checked_scaled_sum};
use super::{WaldContrastOutput, WaldTestOptions};
use crate::errors::{DeseqError, invalid_dimensions};
use crate::matrix::{MatrixOrder, MatrixView};

/// Shape-checked per-gene covariance matrices with explicit tensor layout.
#[derive(Clone, Copy, Debug)]
pub struct CovarianceView<'a> {
    values: MatrixView<'a, f64>,
    coefficients: usize,
}
impl<'a> CovarianceView<'a> {
    /// Wrap G x P² storage. Row-major storage flattens each P x P matrix by
    /// rows; column-major storage is a column-major G x P x P array.
    pub fn new(coefficients: usize, values: MatrixView<'a, f64>) -> Result<Self, DeseqError> {
        let expected = coefficients.checked_mul(coefficients).ok_or_else(|| {
            invalid_dimensions(
                "beta covariance shape overflow",
                usize::MAX,
                values.n_cols(),
            )
        })?;
        if coefficients == 0 || values.n_cols() != expected {
            return Err(invalid_dimensions(
                "beta covariance columns",
                expected,
                values.n_cols(),
            ));
        }
        Ok(Self {
            values,
            coefficients,
        })
    }
    fn at(&self, gene: usize, row: usize, column: usize) -> f64 {
        let index = match self.values.order() {
            MatrixOrder::RowMajor => row * self.coefficients + column,
            MatrixOrder::ColumnMajor => row + column * self.coefficients,
        };
        *self.values.at(gene, index)
    }
}

/// Independently available effect and standard-error projections.
#[derive(Clone, Debug, PartialEq)]
pub struct ContrastProjection {
    /// Estimated log2 fold changes; missing estimates do not discard available SEs.
    pub log2_fold_change: Vec<Option<f64>>,
    /// Standard errors; zero is retained even though it cannot form a Wald test.
    pub lfc_se: Vec<Option<f64>>,
}

/// Borrowed inference arrays, validated once and reusable across requests.
#[derive(Clone, Copy, Debug)]
pub struct InferenceEstimates<'a> {
    beta: MatrixView<'a, f64>,
    beta_se: Option<MatrixView<'a, f64>>,
    covariance: Option<CovarianceView<'a>>,
}
impl<'a> InferenceEstimates<'a> {
    /// Validate aligned shapes without copying beta, SE or covariance storage.
    /// Missing optional variance inputs are diagnosed by the actual request.
    pub fn new(
        beta: MatrixView<'a, f64>,
        beta_se: Option<MatrixView<'a, f64>>,
        covariance: Option<CovarianceView<'a>>,
    ) -> Result<Self, DeseqError> {
        if let Some(se) = beta_se {
            if se.n_rows() != beta.n_rows() {
                return Err(invalid_dimensions(
                    "Wald betaSE rows",
                    beta.n_rows(),
                    se.n_rows(),
                ));
            }
            if se.n_cols() != beta.n_cols() {
                return Err(invalid_dimensions(
                    "Wald betaSE coefficients",
                    beta.n_cols(),
                    se.n_cols(),
                ));
            }
        }
        if let Some(cov) = covariance {
            if cov.values.n_rows() != beta.n_rows() {
                return Err(invalid_dimensions(
                    "beta covariance rows",
                    beta.n_rows(),
                    cov.values.n_rows(),
                ));
            }
            if cov.coefficients != beta.n_cols() {
                return Err(invalid_dimensions(
                    "beta covariance coefficients",
                    beta.n_cols(),
                    cov.coefficients,
                ));
            }
        }
        Ok(Self {
            beta,
            beta_se,
            covariance,
        })
    }
    /// Number of gene rows.
    pub fn n_genes(&self) -> usize {
        self.beta.n_rows()
    }
    /// Number of coefficient columns.
    pub fn n_coefficients(&self) -> usize {
        self.beta.n_cols()
    }

    /// Project one finite nonzero linear contrast using shared arithmetic.
    /// A single active coefficient uses supplied marginal SEs when available;
    /// other requests require covariance. Reuse one P² scratch buffer per call.
    pub fn project_contrast(&self, contrast: &[f64]) -> Result<ContrastProjection, DeseqError> {
        if contrast.len() != self.n_coefficients() {
            return Err(invalid_dimensions(
                "Wald contrast coefficients",
                self.n_coefficients(),
                contrast.len(),
            ));
        }
        let mut active = None;
        let mut nonzero = 0;
        for (i, &weight) in contrast.iter().enumerate() {
            if !weight.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: "Wald contrast".into(),
                    index: Some(i),
                    value: weight,
                });
            }
            if weight != 0.0 {
                nonzero += 1;
                active = Some(i);
            }
        }
        if nonzero == 0 {
            return Err(DeseqError::InvalidOptions {
                reason: "Wald contrast must contain at least one non-zero coefficient".into(),
            });
        }
        let marginal = if nonzero == 1 {
            self.beta_se.map(|se| (active.unwrap(), se))
        } else {
            None
        };
        if marginal.is_none() && self.covariance.is_none() {
            return Err(DeseqError::UnsupportedFeature {feature:"betaCovariance (beta covariance matrices) is required for multi-coefficient or numeric contrast SEs".into()});
        }
        let capacity = if marginal.is_some() {
            contrast.len()
        } else {
            self.covariance.unwrap().values.n_cols()
        };
        let mut scratch = Vec::with_capacity(capacity);
        let mut log2_fold_change = Vec::with_capacity(self.n_genes());
        let mut lfc_se = Vec::with_capacity(self.n_genes());
        for gene in 0..self.n_genes() {
            scratch.clear();
            for (column, &weight) in contrast.iter().enumerate() {
                scratch.push(weight * self.beta.at(gene, column));
            }
            log2_fold_change.push(checked_scaled_sum(&scratch));
            let se = if let Some((column, se)) = marginal {
                let value = contrast[column].abs() * se.at(gene, column);
                (value.is_finite() && value >= 0.0).then_some(value)
            } else {
                scratch.clear();
                let cov = self.covariance.unwrap();
                let mut valid = true;
                'covariance: for (row, &left) in contrast.iter().enumerate() {
                    for (column, &right) in contrast.iter().enumerate() {
                        if let Some(term) = checked_product3_ordered_by_magnitude(
                            left,
                            cov.at(gene, row, column),
                            right,
                        ) {
                            scratch.push(term);
                        } else {
                            valid = false;
                            break 'covariance;
                        }
                    }
                }
                if valid {
                    checked_scaled_sum(&scratch)
                        .filter(|variance| *variance >= 0.0)
                        .map(f64::sqrt)
                        .filter(|se| se.is_finite())
                } else {
                    None
                }
            };
            lfc_se.push(se);
        }
        Ok(ContrastProjection {
            log2_fold_change,
            lfc_se,
        })
    }
}

/// Compute a Wald contrast directly from borrowed primitive arrays.
/// Explicit scalar/per-gene df is required without fitted design context.
/// Invalid inference rows mask the effect and SE as well as the statistic;
/// missing t df alone mask only the p-value. Use `project_contrast` when raw
/// effect and SE availability must remain independent of the Wald test.
pub fn wald_test_contrast_from_estimates(
    estimates: &InferenceEstimates<'_>,
    contrast: &[f64],
    options: &WaldTestOptions,
) -> Result<WaldContrastOutput, DeseqError> {
    wald_contrast_with_context(estimates, contrast, options, None)
}

pub(super) fn wald_contrast_with_context(
    estimates: &InferenceEstimates<'_>,
    contrast: &[f64],
    options: &WaldTestOptions,
    residual_df: Option<f64>,
) -> Result<WaldContrastOutput, DeseqError> {
    let mut projected = estimates.project_contrast(contrast)?;
    let rows = projected
        .log2_fold_change
        .iter()
        .zip(&projected.lfc_se)
        .map(|(&estimate, &se)| (estimate.unwrap_or(f64::NAN), se.unwrap_or(f64::NAN)));
    let wald = wald_from_estimate_rows(rows, options, residual_df)?;
    // A Wald contrast is a joint inference row. Preserve the standalone rule:
    // unavailable estimates, variance or statistic mask both displayed fields.
    // Missing t df retain the statistic, and consequently the effect and SE.
    for (gene, stat) in wald.stat.iter().enumerate() {
        if stat.is_none() {
            projected.log2_fold_change[gene] = None;
            projected.lfc_se[gene] = None;
        }
    }
    Ok(WaldContrastOutput {
        log2_fold_change: projected.log2_fold_change,
        lfc_se: projected.lfc_se,
        wald,
    })
}
