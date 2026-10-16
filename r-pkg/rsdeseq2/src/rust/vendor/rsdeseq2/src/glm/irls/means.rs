// Starting coefficients and stable fitted means.

fn initial_beta(x: &DMatrix<f64>, y: &[f64], nf: &[f64]) -> Result<DVector<f64>, DeseqError> {
    initial_beta_impl(x, y, nf, None)
}

fn initial_beta_impl(
    x: &DMatrix<f64>,
    y: &[f64],
    nf: &[f64],
    initial_lu: Option<&LU<f64, Dyn, Dyn>>,
) -> Result<DVector<f64>, DeseqError> {
    let response = y
        .iter()
        .copied()
        .zip(nf.iter().copied())
        .map(|(count, factor)| {
            validate_positive_finite(factor, "normalization factor", 0)?;
            Ok((count / factor + 0.1).ln())
        })
        .collect::<Result<Vec<_>, DeseqError>>()?;
    let xty = x.transpose() * DVector::from_vec(response);
    let beta = match initial_lu {
        Some(lu) => lu.solve(&xty),
        None => (x.transpose() * x).lu().solve(&xty),
    };
    if let Some(beta) = beta {
        Ok(beta)
    } else {
        Ok(DVector::from_element(x.ncols(), 0.0))
    }
}

#[cfg(test)]
fn fitted_mu(
    x: &DMatrix<f64>,
    beta: &DVector<f64>,
    nf: &[f64],
    min_mu: f64,
) -> Result<Vec<f64>, DeseqError> {
    fitted_mu_impl(x, beta, nf, Some(min_mu))
}

fn fitted_mu_log2_unfloored(
    x: &DMatrix<f64>,
    beta: &[f64],
    nf: &[f64],
) -> Result<Vec<f64>, DeseqError> {
    if beta.len() != x.ncols() {
        return Err(invalid_dimensions(
            "log2 beta coefficients",
            x.ncols(),
            beta.len(),
        ));
    }
    (0..x.nrows())
        .map(|sample| {
            validate_positive_finite(nf[sample], "normalization factor", sample)?;
            let eta = checked_row_dot_slice(x, sample, beta).ok_or_else(|| {
                DeseqError::NonFiniteValue {
                    context: "optim fallback linear predictor".to_string(),
                    index: Some(sample),
                    value: f64::NAN,
                }
            })?;
            let mu = nf[sample] * 2.0_f64.powf(eta);
            if !mu.is_finite() || mu <= 0.0 {
                return Err(DeseqError::NonFiniteValue {
                    context: "optim fallback fitted mean".to_string(),
                    index: Some(sample),
                    value: mu,
                });
            }
            Ok(mu)
        })
        .collect()
}

#[cfg(test)]
fn fitted_mu_impl(
    x: &DMatrix<f64>,
    beta: &DVector<f64>,
    nf: &[f64],
    min_mu: Option<f64>,
) -> Result<Vec<f64>, DeseqError> {
    fitted_mu_impl_with_logs(x, beta, nf, min_mu, None)
}

fn fitted_mu_impl_with_logs(
    x: &DMatrix<f64>,
    beta: &DVector<f64>,
    nf: &[f64],
    min_mu: Option<f64>,
    log_nf: Option<&[f64]>,
) -> Result<Vec<f64>, DeseqError> {
    let mut output = Vec::with_capacity(x.nrows());
    fitted_mu_into(x, beta, nf, min_mu, log_nf, &mut output)?;
    Ok(output)
}

fn fitted_mu_into(
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
            let eta =
                checked_row_dot_vector(x, row, beta).ok_or_else(|| DeseqError::NonFiniteValue {
                    context: "IRLS linear predictor".to_string(),
                    index: Some(sample),
                    value: f64::NAN,
                })?;
            let log_factor = log_nf.map(|values| values[sample]);
            let mu = finite_scaled_exp_with_log(factor, eta, log_factor).ok_or_else(|| {
                DeseqError::NonFiniteValue {
                    context: "IRLS fitted mean".to_string(),
                    index: Some(sample),
                    value: f64::NAN,
                }
            })?;
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

fn checked_row_dot_slice(x: &DMatrix<f64>, row: usize, beta: &[f64]) -> Option<f64> {
    let mut sum = 0.0;
    for col in 0..x.ncols() {
        sum = x[(row, col)].mul_add(beta[col], sum);
        if !sum.is_finite() {
            return None;
        }
    }
    Some(sum)
}

fn checked_row_dot_vector(x: &DMatrix<f64>, row: usize, beta: &DVector<f64>) -> Option<f64> {
    let mut sum = ScaledSum::default();
    for col in 0..x.ncols() {
        sum.add(checked_product2(x[(row, col)], beta[col])?)?;
    }
    sum.finish()
}

#[derive(Default)]
struct ScaledSum {
    scale: f64,
    normalized_sum: f64,
}

impl ScaledSum {
    fn add(&mut self, value: f64) -> Option<()> {
        if !value.is_finite() {
            return None;
        }
        let abs = value.abs();
        if abs == 0.0 {
            return Some(());
        }
        if self.scale == 0.0 {
            self.scale = abs;
            self.normalized_sum = value / abs;
            return Some(());
        }
        if abs > self.scale {
            let rescaled = self.normalized_sum * (self.scale / abs);
            self.normalized_sum = checked_sum2(rescaled, value / abs)?;
            self.scale = abs;
        } else {
            self.normalized_sum = checked_sum2(self.normalized_sum, value / self.scale)?;
        }
        Some(())
    }

    fn finish(self) -> Option<f64> {
        if self.scale == 0.0 {
            return Some(0.0);
        }
        checked_product2(self.normalized_sum, self.scale)
    }
}

fn finite_scaled_exp_with_log(factor: f64, eta: f64, log_factor: Option<f64>) -> Option<f64> {
    if !factor.is_finite() || factor <= 0.0 || !eta.is_finite() {
        return None;
    }
    let log_mu = log_factor.unwrap_or_else(|| factor.ln()) + eta;
    if log_mu >= f64::MAX.ln() {
        return Some(f64::MAX);
    }
    if log_mu <= f64::MIN_POSITIVE.ln() {
        return Some(f64::MIN_POSITIVE);
    }
    let mu = log_mu.exp();
    (mu.is_finite() && mu > 0.0).then_some(mu)
}
