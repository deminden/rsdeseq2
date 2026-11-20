use super::*;
use crate::numeric::math::ProfileMath;
use nalgebra::{DVectorView, Dyn, linalg::LU};
// Starting coefficients and stable fitted means.

pub(super) fn initial_beta(
    x: &DMatrix<f64>,
    y: &[f64],
    nf: &[f64],
) -> Result<DVector<f64>, DeseqError> {
    let mut beta = try_zero_dvector(x.ncols(), "IRLS starting coefficients")?;
    let mut response = try_alloc_matrix(y.len(), 1, "IRLS starting response")?;
    InitialBeta::least_squares(x)?.initialize(y, nf, &mut response, &mut beta)?;
    Ok(beta)
}

/// How `fitNbinomGLMs` (DESeq2 1.52, R/fitNbinomGLMs.R) starts IRLS.
pub(super) enum InitialBeta {
    /// Full-rank design (`qr(x)$rank == ncol(x)`): least squares on
    /// `log(normalized count + 0.1)`.
    LeastSquares {
        transpose: DMatrix<f64>,
        lu: LU<f64, Dyn, Dyn>,
    },
    /// Rank-deficient design with an `Intercept` column: that coefficient
    /// starts at `log(rowMeans(normalized counts))`, the others at 0.
    LogBaseMeanIntercept(usize),
    /// Rank-deficient design without an intercept: every coefficient is 1.
    Ones,
}

impl InitialBeta {
    pub(super) fn for_design(design: &DesignMatrix, x: &DMatrix<f64>) -> Result<Self, DeseqError> {
        if design.rank()? == design.n_coefficients() {
            return Self::least_squares(x);
        }
        let intercept = design.coefficient_names().and_then(|names| {
            names
                .iter()
                .position(|name| name == "Intercept" || name == "(Intercept)")
        });
        Ok(intercept.map_or(Self::Ones, Self::LogBaseMeanIntercept))
    }

    fn least_squares(x: &DMatrix<f64>) -> Result<Self, DeseqError> {
        let coefficients = x.ncols();
        let mut transpose =
            try_zero_dmatrix(coefficients, x.nrows(), "IRLS starting design transpose")?;
        x.transpose_to(&mut transpose);
        let mut crossproduct = try_zero_dmatrix(
            coefficients,
            coefficients,
            "IRLS starting design crossproduct",
        )?;
        transpose.mul_to(x, &mut crossproduct);
        let lu = crossproduct.lu();
        Ok(Self::LeastSquares { transpose, lu })
    }

    pub(super) fn initialize(
        &self,
        y: &[f64],
        nf: &[f64],
        response: &mut Vec<f64>,
        beta: &mut DVector<f64>,
    ) -> Result<(), DeseqError> {
        match self {
            Self::LeastSquares { transpose, lu } => {
                response.clear();
                for (count, factor) in y.iter().copied().zip(nf.iter().copied()) {
                    validate_positive_finite(factor, "normalization factor", 0)?;
                    response.push((count / factor + 0.1).ln_p());
                }
                transpose.mul_to(&DVectorView::from_slice(response, response.len()), beta);
                if !lu.solve_mut(beta) {
                    beta.fill(0.0);
                }
            }
            Self::Ones => beta.fill(1.0),
            Self::LogBaseMeanIntercept(intercept) => {
                let mut total = 0.0;
                for (count, factor) in y.iter().copied().zip(nf.iter().copied()) {
                    validate_positive_finite(factor, "normalization factor", 0)?;
                    total += count / factor;
                }
                beta.fill(0.0);
                beta[*intercept] = (total / y.len() as f64).ln_p();
            }
        }
        Ok(())
    }
}

#[cfg(test)]
pub(super) fn fitted_mu(
    x: &DMatrix<f64>,
    beta: &DVector<f64>,
    nf: &[f64],
    min_mu: f64,
) -> Result<Vec<f64>, DeseqError> {
    fitted_mu_impl(x, beta, nf, Some(min_mu))
}

#[cfg(test)]
pub(super) fn fitted_mu_log2_unfloored(
    x: &DMatrix<f64>,
    beta: &[f64],
    nf: &[f64],
) -> Result<Vec<f64>, DeseqError> {
    let mut output = vec![0.0; x.nrows()];
    fitted_mu_log2_unfloored_into(x, beta, nf, &mut output)?;
    Ok(output)
}

pub(super) fn fitted_mu_log2_unfloored_into(
    x: &DMatrix<f64>,
    beta: &[f64],
    nf: &[f64],
    output: &mut [f64],
) -> Result<(), DeseqError> {
    if beta.len() != x.ncols() {
        return Err(invalid_dimensions(
            "log2 beta coefficients",
            x.ncols(),
            beta.len(),
        ));
    }
    for (sample, value) in output.iter_mut().enumerate() {
        validate_positive_finite(nf[sample], "normalization factor", sample)?;
        let eta = checked_row_dot_slice(x, sample, beta)
            .nonfinite_at("optim fallback linear predictor", Some(sample))?;
        let mu = nf[sample] * 2.0_f64.powf_p(eta);
        if !mu.is_finite() || mu <= 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: "optim fallback fitted mean".to_string(),
                index: Some(sample),
                value: mu,
            });
        }
        *value = mu;
    }
    Ok(())
}

#[cfg(test)]
pub(super) fn fitted_mu_impl(
    x: &DMatrix<f64>,
    beta: &DVector<f64>,
    nf: &[f64],
    min_mu: Option<f64>,
) -> Result<Vec<f64>, DeseqError> {
    let mut output = Vec::with_capacity(x.nrows());
    fitted_mu_into(x, beta, nf, min_mu, None, &mut output)?;
    Ok(output)
}

pub(super) fn fitted_mu_into(
    x: &DMatrix<f64>,
    beta: &DVector<f64>,
    nf: &[f64],
    min_mu: Option<f64>,
    log_nf: Option<&[f64]>,
    output: &mut Vec<f64>,
) -> Result<(), DeseqError> {
    output.clear();
    let values = (0..x.nrows()).zip(nf.iter().copied()).enumerate().map(
        |(sample, (row, factor))| -> Result<f64, DeseqError> {
            validate_positive_finite(factor, "normalization factor", sample)?;
            let eta = checked_row_dot_vector(x, row, beta)
                .nonfinite_at("IRLS linear predictor", Some(sample))?;
            let log_factor = log_nf.map(|values| values[sample]);
            let mu = finite_scaled_exp_with_log(factor, eta, log_factor)
                .nonfinite_at("IRLS fitted mean", Some(sample))?;
            if mu <= 0.0 {
                return Err(DeseqError::NonFiniteValue {
                    context: "IRLS fitted mean".to_string(),
                    index: Some(sample),
                    value: mu,
                });
            }
            Ok(min_mu.map_or(mu, |min_mu| mu.max(min_mu)))
        },
    );
    for value in values {
        output.push(value?);
    }
    Ok(())
}

pub(super) fn checked_row_dot_slice(
    x: &DMatrix<f64>,
    row: usize,
    beta: &[f64],
) -> checked::Checked {
    let mut sum = 0.0;
    for col in 0..x.ncols() {
        sum = x[(row, col)].mul_add(beta[col], sum);
        if !sum.is_finite() {
            return checked::reject(sum);
        }
    }
    Ok(sum)
}

pub(super) fn checked_row_dot_vector(
    x: &DMatrix<f64>,
    row: usize,
    beta: &DVector<f64>,
) -> checked::Checked {
    let mut sum = ScaledSum::default();
    for col in 0..x.ncols() {
        sum.add(checked::try_mul(x[(row, col)], beta[col])?)?;
    }
    sum.finish()
}

#[derive(Default)]
pub(super) struct ScaledSum {
    scale: f64,
    normalized_sum: f64,
}

impl ScaledSum {
    fn add(&mut self, value: f64) -> Result<(), checked::NonFinite> {
        if !value.is_finite() {
            checked::reject(value)?;
        }
        let abs = value.abs();
        if abs == 0.0 {
            return Ok(());
        }
        if self.scale == 0.0 {
            self.scale = abs;
            self.normalized_sum = value / abs;
            return Ok(());
        }
        if abs > self.scale {
            let rescaled = self.normalized_sum * (self.scale / abs);
            self.normalized_sum = checked::try_add(rescaled, value / abs)?;
            self.scale = abs;
        } else {
            self.normalized_sum = checked::try_add(self.normalized_sum, value / self.scale)?;
        }
        Ok(())
    }

    fn finish(self) -> checked::Checked {
        if self.scale == 0.0 {
            return Ok(0.0);
        }
        checked::try_mul(self.normalized_sum, self.scale)
    }
}

pub(super) fn finite_scaled_exp_with_log(
    factor: f64,
    eta: f64,
    log_factor: Option<f64>,
) -> checked::Checked {
    if !factor.is_finite() || factor <= 0.0 {
        return checked::reject(factor);
    }
    if !eta.is_finite() {
        return checked::reject(eta);
    }
    let log_mu = log_factor.unwrap_or_else(|| factor.ln_p()) + eta;
    if log_mu >= f64::MAX.ln_p() {
        return Ok(f64::MAX);
    }
    if log_mu <= f64::MIN_POSITIVE.ln_p() {
        return Ok(f64::MIN_POSITIVE);
    }
    let mu = log_mu.exp_p();
    if mu.is_finite() && mu > 0.0 {
        Ok(mu)
    } else {
        checked::reject(mu)
    }
}
