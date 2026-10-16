// Projected-gradient convergence and checked arithmetic.

fn projected_gradient_norm(
    parameters: &[f64],
    gradient: &[f64],
    lower: f64,
    upper: f64,
) -> Option<f64> {
    let scale = parameters
        .iter()
        .copied()
        .zip(gradient.iter().copied())
        .map(|(parameter, gradient)| {
            if (parameter <= lower && gradient > 0.0) || (parameter >= upper && gradient < 0.0) {
                0.0
            } else {
                gradient.abs()
            }
        })
        .try_fold(0.0_f64, |scale, value| {
            value.is_finite().then_some(scale.max(value))
        })?;
    if scale == 0.0 {
        return Some(0.0);
    }
    let mut sum = 0.0;
    for (parameter, gradient) in parameters.iter().copied().zip(gradient.iter().copied()) {
        let value =
            if (parameter <= lower && gradient > 0.0) || (parameter >= upper && gradient < 0.0) {
                0.0
            } else {
                gradient
            };
        let scaled = checked_div2(value, scale)?;
        let term = checked_product2(scaled, scaled)?;
        let next = sum + term;
        if !term.is_finite() || !next.is_finite() {
            return None;
        }
        sum = next;
    }
    let norm = scale * sum.sqrt();
    norm.is_finite().then_some(norm)
}

fn checked_sum2(left: f64, right: f64) -> Option<f64> {
    let sum = left + right;
    (left.is_finite() && right.is_finite() && sum.is_finite()).then_some(sum)
}

fn checked_product2(left: f64, right: f64) -> Option<f64> {
    let product = left * right;
    (left.is_finite() && right.is_finite() && product.is_finite()).then_some(product)
}

fn checked_div2(left: f64, right: f64) -> Option<f64> {
    let quotient = left / right;
    (left.is_finite() && right.is_finite() && right != 0.0 && quotient.is_finite())
        .then_some(quotient)
}

fn checked_product3(left: f64, middle: f64, right: f64) -> Option<f64> {
    checked_product2(checked_product2(left, middle)?, right)
}

fn checked_product4(first: f64, second: f64, third: f64, fourth: f64) -> Option<f64> {
    checked_product2(
        checked_product2(first, second)?,
        checked_product2(third, fourth)?,
    )
}

fn checked_product6(
    first: f64,
    second: f64,
    third: f64,
    fourth: f64,
    fifth: f64,
    sixth: f64,
) -> Option<f64> {
    checked_product2(
        checked_product3(first, second, third)?,
        checked_product3(fourth, fifth, sixth)?,
    )
}

fn irls_deviance_convergence_stat(dev: f64, dev_old: f64) -> Option<f64> {
    let delta = checked_sum2(dev, -dev_old)?.abs();
    let scale = checked_sum2(dev.abs(), 0.1)?;
    checked_div2(delta, scale)
}
