/// Resolve a Cook's cutoff option for a model with `m` samples and `p` columns.
pub fn resolve_cooks_cutoff(
    cutoff: CooksCutoff,
    n_samples: usize,
    n_coefficients: usize,
) -> Result<Option<f64>, DeseqError> {
    match cutoff {
        CooksCutoff::Disabled => Ok(None),
        CooksCutoff::Threshold(value) => {
            if !value.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: "Cook's cutoff".to_string(),
                    index: None,
                    value,
                });
            }
            Ok(Some(value))
        }
        CooksCutoff::Default => default_cooks_cutoff(n_samples, n_coefficients),
    }
}

/// DESeq2 default Cook's cutoff, `qf(.99, p, m - p)`.
///
/// Includes R's established extreme-df chi-square limit policy, evaluated by
/// the same standalone Rust probability implementation for both interfaces.
pub fn default_cooks_cutoff(
    n_samples: usize,
    n_coefficients: usize,
) -> Result<Option<f64>, DeseqError> {
    if n_samples <= n_coefficients {
        return Ok(None);
    }
    let df1 = n_coefficients as f64;
    let df2 = (n_samples - n_coefficients) as f64;
    if !df1.is_finite() || df1 <= 0.0 || !df2.is_finite() || df2 <= 0.0 {
        return Err(DeseqError::InvalidDimensions {
            context: "Cook's cutoff degrees of freedom".to_string(),
            expected: n_samples,
            actual: n_coefficients,
        });
    }
    let cutoff = crate::probability::f_99_quantile(df1, df2);
    if !cutoff.is_finite() {
        return Err(DeseqError::NonFiniteValue {
            context: "Cook's default cutoff".to_string(),
            index: None,
            value: cutoff,
        });
    }
    Ok(Some(cutoff))
}

/// Apply Cook's p-value filtering and recompute BH-adjusted p-values.
///
/// This mirrors the default `results()` behavior where rows with `maxCooks`
/// above the selected cutoff have their p-value set to missing before
/// p-value adjustment. Count replacement is handled separately.
pub fn apply_cooks_cutoff(
    results: &mut DeseqResults,
    cutoff: Option<f64>,
) -> Result<(), DeseqError> {
    if cutoff.is_none() {
        return Ok(());
    }
    let policy = crate::results::CooksPolicy::new(cutoff, None)?;
    apply_cooks_policy(results, policy)
}

/// Apply Cook's p-value filtering with DESeq2's two-group low-count heuristic.
///
/// DESeq2's `results()` applies this heuristic only for formula designs with a
/// single two-level factor. The Rust core cannot infer that from primitive
/// matrices, so callers should use this helper only after establishing that
/// condition. For rows above the Cook's cutoff, the helper finds the sample
/// with the maximum Cook's distance. If at least three counts in that row are
/// larger than the count in that sample, the row is not masked.
pub fn apply_cooks_cutoff_with_low_count_heuristic(
    results: &mut DeseqResults,
    cutoff: Option<f64>,
    counts: &CountMatrix,
    cooks: &RowMajorMatrix<f64>,
) -> Result<(), DeseqError> {
    if cutoff.is_none() {
        return Ok(());
    }
    // Validate threshold before matrix context, as on the legacy path.
    crate::results::CooksPolicy::new(cutoff, None)?;
    let heuristic = crate::results::LowCountHeuristic::new(
        crate::matrix::CountView::from(counts),
        cooks.into(),
    )?;
    let policy = crate::results::CooksPolicy::new(cutoff, Some(heuristic))?;
    apply_cooks_policy(results, policy)
}

fn apply_cooks_policy(
    results: &mut DeseqResults,
    policy: crate::results::CooksPolicy<'_>,
) -> Result<(), DeseqError> {
    update_cooks_policy_rows(results, policy, true)
}

pub(crate) fn mask_cooks_policy_rows(
    results: &mut DeseqResults,
    policy: crate::results::CooksPolicy<'_>,
) -> Result<(), DeseqError> {
    update_cooks_policy_rows(results, policy, false)
}

fn update_cooks_policy_rows(
    results: &mut DeseqResults,
    policy: crate::results::CooksPolicy<'_>,
    adjust: bool,
) -> Result<(), DeseqError> {
    let mask = policy.evaluate(
        results.rows.iter().map(|r| r.pvalue),
        results.rows.iter().map(|r| r.max_cooks),
    )?;
    let adjusted = adjust.then(|| bh_adjust(&mask.pvalue));
    for (index, ((row, pvalue), outlier)) in results
        .rows
        .iter_mut()
        .zip(mask.pvalue)
        .zip(mask.outlier)
        .enumerate()
    {
        row.pvalue = pvalue;
        row.cooks_outlier = outlier;
        row.padj = adjusted.as_ref().and_then(|v| v[index]);
    }
    Ok(())
}

/// Recompute BH-adjusted p-values from the current result p-values.
pub fn recompute_padj(results: &mut DeseqResults) -> Result<(), DeseqError> {
    let pvalues = results
        .rows
        .iter()
        .map(|row| row.pvalue)
        .collect::<Vec<_>>();
    validate_optional_probability(&pvalues, "result p-value")?;
    let padj = bh_adjust(&pvalues);
    for (row, adjusted) in results.rows.iter_mut().zip(padj) {
        row.padj = adjusted;
    }
    Ok(())
}
