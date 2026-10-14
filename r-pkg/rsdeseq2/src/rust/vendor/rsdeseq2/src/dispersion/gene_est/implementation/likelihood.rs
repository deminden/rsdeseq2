/// DESeq2's alpha-dependent NB log-likelihood kernel.
///
/// Terms independent of alpha are omitted, matching the objective used inside
/// DESeq2's dispersion optimizer.
pub fn dispersion_nb_log_likelihood_kernel(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    dispersion_nb_log_likelihood_kernel_weighted(counts, mu, log_alpha, None)
}

/// DESeq2's alpha-dependent NB log-likelihood kernel with optional observation weights.
///
/// Terms independent of alpha are omitted, matching the objective used inside
/// DESeq2's dispersion optimizer. When supplied, observation weights multiply
/// the per-sample terms.
pub fn dispersion_nb_log_likelihood_kernel_weighted(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    if counts.len() != mu.len() {
        return Err(invalid_dimensions(
            "dispersion objective mu",
            counts.len(),
            mu.len(),
        ));
    }
    validate_observation_weight_slice(weights, counts.len(), "dispersion objective weights")?;
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
    let inv_alpha = alpha.recip();
    let mut total = 0.0;
    for (sample, (count, mu)) in counts.iter().copied().zip(mu.iter().copied()).enumerate() {
        validate_positive_mu(mu, sample)?;
        let observation_weight = weights.map(|values| values[sample]).unwrap_or(1.0);
        let y = f64::from(count);
        let mu_alpha = mu * alpha;
        let mu_plus_inv_alpha = checked_add(
            mu,
            inv_alpha,
            sample,
            "dispersion objective mean plus inverse alpha",
        )?;
        let term = dispersion_lgamma(y + inv_alpha)
            - dispersion_lgamma(inv_alpha)
            - y * mu_plus_inv_alpha.ln()
            - inv_alpha * (1.0 + mu_alpha).ln();
        if !term.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "dispersion objective likelihood term".to_string(),
                index: Some(sample),
                value: term,
            });
        }
        let weighted_term = checked_mul(
            observation_weight,
            term,
            sample,
            "dispersion objective weighted likelihood term",
        )?;
        total = checked_add(
            total,
            weighted_term,
            sample,
            "dispersion objective likelihood sum",
        )?;
    }
    Ok(total)
}

fn dispersion_lgamma(x: f64) -> f64 {
    if x > 10.0 {
        stirling_lgamma_positive(x)
    } else {
        ln_gamma(x)
    }
}

fn stirling_lgamma_positive(x: f64) -> f64 {
    const LN_SQRT_2PI: f64 = 0.918_938_533_204_672_7;
    if x > 1.0e17 {
        return x * (x.ln() - 1.0);
    }
    let base = LN_SQRT_2PI + (x - 0.5) * x.ln() - x;
    if x > 4_934_720.0 {
        base
    } else {
        base + stirling_log_gamma_correction(x)
    }
}

fn stirling_log_gamma_correction(x: f64) -> f64 {
    const LGAMMA_COR_COEFFS: [f64; 5] = [
        0.166_638_948_045_186_3,
        -0.000_013_849_481_760_675_638,
        0.000_000_009_810_825_646_924_73,
        -0.000_000_000_018_091_294_755_724_94,
        0.000_000_000_000_062_210_980_418_926_05,
    ];
    let tmp = 10.0 / x;
    chebyshev_eval(tmp * tmp * 2.0 - 1.0, &LGAMMA_COR_COEFFS) / x
}

fn chebyshev_eval(x: f64, coeffs: &[f64]) -> f64 {
    let twox = x * 2.0;
    let mut b2 = 0.0;
    let mut b1 = 0.0;
    let mut b0 = 0.0;
    for coeff in coeffs.iter().rev().copied() {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + coeff;
    }
    (b0 - b2) * 0.5
}

fn dispersion_digamma(x: f64) -> f64 {
    if x <= 0.0 || !x.is_finite() {
        return digamma(x);
    }
    let xln = x.ln();
    if x * xln > 1.0 / (2.0 * f64::EPSILON) {
        return xln;
    }

    const WDTOL: f64 = 1.110_223_024_625_156_5e-16;
    if x < WDTOL {
        return -x.recip();
    }

    const BVALUES: [f64; 22] = [
        1.0,
        -0.5,
        0.166_666_666_666_666_66,
        -0.033_333_333_333_333_33,
        0.023_809_523_809_523_808,
        -0.033_333_333_333_333_33,
        0.075_757_575_757_575_76,
        -0.253_113_553_113_553_1,
        1.166_666_666_666_666_7,
        -7.092_156_862_745_098,
        54.971_177_944_862_156,
        -529.124_242_424_242_4,
        6_192.123_188_405_797,
        -86_580.253_113_553_12,
        1_425_517.166_666_666_7,
        -27_298_231.067_816_09,
        601_580_873.900_642_4,
        -15_116_315_767.092_157,
        429_614_643_061.166_7,
        -13_711_655_205_088.332,
        488_332_318_973_593.2,
        -19_296_579_341_940_068.0,
    ];
    let rln = (std::f64::consts::LOG10_2 * 53.0).min(18.06);
    let fln = rln.max(3.0) - 3.0;
    let xmin = (3.50 + 0.40 * fln) as i32 + 1;
    let (xdmy, xdmln, xinc) = if x < f64::from(xmin) {
        let lower = x as i32;
        let increment = f64::from(xmin - lower);
        let shifted = x + increment;
        (shifted, shifted.ln(), increment as i32)
    } else {
        (x, xln, 0)
    };

    let tt = 0.5 / xdmy;
    let tst = WDTOL * tt;
    let rxsq = 1.0 / (xdmy * xdmy);
    let ta = 0.5 * rxsq;
    let mut term = ta;
    let mut sum = term * BVALUES[2];
    if sum.abs() >= tst {
        let mut tk = 2.0;
        for coeff in BVALUES.iter().take(22).skip(3).copied() {
            term = term * (tk / (tk + 2.0)) * rxsq;
            let contribution = term * coeff;
            if contribution.abs() < tst {
                break;
            }
            sum += contribution;
            tk += 2.0;
        }
    }
    sum += tt;
    if xinc != 0 {
        for i in 1..=xinc {
            sum += 1.0 / (x + f64::from(xinc - i));
        }
    }
    xdmln - sum
}

/// Derivative of DESeq2's alpha-dependent NB likelihood kernel with respect to log alpha.
pub fn dispersion_nb_log_likelihood_kernel_derivative(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    dispersion_nb_log_likelihood_kernel_derivative_weighted(counts, mu, log_alpha, None)
}

/// Derivative of the weighted NB likelihood kernel with respect to log alpha.
pub fn dispersion_nb_log_likelihood_kernel_derivative_weighted(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    if counts.len() != mu.len() {
        return Err(invalid_dimensions(
            "dispersion objective derivative mu",
            counts.len(),
            mu.len(),
        ));
    }
    validate_observation_weight_slice(
        weights,
        counts.len(),
        "dispersion objective derivative weights",
    )?;
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
    let inv_alpha = alpha.recip();
    let mut derivative_alpha = 0.0;
    for (sample, (count, mu)) in counts.iter().copied().zip(mu.iter().copied()).enumerate() {
        validate_positive_mu(mu, sample)?;
        let observation_weight = weights.map(|values| values[sample]).unwrap_or(1.0);
        let y = f64::from(count);
        let mu_alpha = mu_alpha_terms(mu, alpha, sample, "dispersion objective derivative")?;
        let term = dispersion_digamma(inv_alpha)
            + mu_alpha.log1p
            - mu_alpha.ratio
            - dispersion_digamma(y + inv_alpha)
            + y * mu_alpha.alpha_over_one_plus;
        if !term.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "dispersion objective derivative term".to_string(),
                index: Some(sample),
                value: term,
            });
        }
        let weighted_term = checked_mul(
            observation_weight,
            term,
            sample,
            "dispersion objective weighted derivative term",
        )?;
        derivative_alpha = checked_add(
            derivative_alpha,
            weighted_term,
            sample,
            "dispersion objective derivative sum",
        )?;
    }
    checked_log_alpha_first_derivative(
        inv_alpha,
        derivative_alpha,
        "dispersion objective log-alpha derivative",
    )
}

/// Second derivative of DESeq2's NB likelihood kernel with respect to log alpha.
pub fn dispersion_nb_log_likelihood_kernel_second_derivative(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    dispersion_nb_log_likelihood_kernel_second_derivative_weighted(counts, mu, log_alpha, None)
}

/// Second derivative of the weighted NB likelihood kernel with respect to log alpha.
pub fn dispersion_nb_log_likelihood_kernel_second_derivative_weighted(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    if counts.len() != mu.len() {
        return Err(invalid_dimensions(
            "dispersion objective second derivative mu",
            counts.len(),
            mu.len(),
        ));
    }
    validate_observation_weight_slice(
        weights,
        counts.len(),
        "dispersion objective second derivative weights",
    )?;
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
    let inv_alpha = alpha.recip();
    let inv_alpha_squared = inv_alpha * inv_alpha;
    let mut first_alpha_sum = 0.0;
    let mut second_alpha_sum = 0.0;
    for (sample, (count, mu)) in counts.iter().copied().zip(mu.iter().copied()).enumerate() {
        validate_positive_mu(mu, sample)?;
        let observation_weight = weights.map(|values| values[sample]).unwrap_or(1.0);
        let y = f64::from(count);
        let mu_alpha = mu_alpha_terms(mu, alpha, sample, "dispersion objective second derivative")?;
        let first_term =
            dispersion_digamma(inv_alpha) + mu_alpha.log1p
                - mu_alpha.ratio
                - dispersion_digamma(y + inv_alpha)
                + y * mu_alpha.alpha_over_one_plus;
        let second_term = -inv_alpha_squared * trigamma(inv_alpha)?
            + mu_alpha.mu_squared_alpha_over_one_plus_squared
            + inv_alpha_squared * trigamma(y + inv_alpha)?
            + y * mu_alpha.inv_one_plus_squared;
        if !first_term.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "dispersion objective second derivative first term".to_string(),
                index: Some(sample),
                value: first_term,
            });
        }
        if !second_term.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "dispersion objective second derivative term".to_string(),
                index: Some(sample),
                value: second_term,
            });
        }
        checked_matrix_add_assign(
            &mut first_alpha_sum,
            checked_mul(
                observation_weight,
                first_term,
                sample,
                "dispersion objective weighted first derivative term",
            )?,
            sample,
            "dispersion objective first derivative sum",
        )?;
        checked_matrix_add_assign(
            &mut second_alpha_sum,
            checked_mul(
                observation_weight,
                second_term,
                sample,
                "dispersion objective weighted second derivative term",
            )?,
            sample,
            "dispersion objective second derivative sum",
        )?;
    }
    let first_log_alpha =
        dispersion_nb_log_likelihood_kernel_derivative_weighted(counts, mu, log_alpha, weights)?;
    checked_log_alpha_second_derivative(
        second_alpha_sum,
        inv_alpha,
        first_alpha_sum,
        first_log_alpha,
        "dispersion objective log-alpha second derivative",
    )
}
