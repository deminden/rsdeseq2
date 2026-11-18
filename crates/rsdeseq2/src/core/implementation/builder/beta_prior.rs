//! DESeq2's `DESeq(betaPrior = TRUE)` for [`DeseqBuilder::run`].
//!
//! This follows `nbinomWaldTest(betaPrior = TRUE)` with
//! `modelMatrixType = "expanded"` (`fitGLMsWithPrior`,
//! `estimateBetaPriorVar`, `addAllContrasts`, `averagePriorsOverLevels`):
//!
//! 1. Fit the MLE GLM on the standard (treatment-coded) design. Its fitted
//!    means and hat diagonals give the Cook's distances.
//! 2. Take the MLE log2 fold changes, add every pairwise difference between
//!    the non-reference levels of each factor, and match a weighted upper
//!    quantile per column to a zero-centred normal variance.
//! 3. Average those variances over each factor's columns and assign the
//!    average to every level column of the expanded design (one indicator
//!    per level); the intercept keeps the wide prior.
//! 4. Refit on the expanded design with ridge `1 / variance` and report the
//!    requested standard-design column as the equivalent expanded contrast
//!    (for a treatment column `f_X_vs_R`: level `X` minus level `R`).

use super::*;
use crate::all_zero::nonzero_gene_indices;
use crate::design::ExpandedFormulaDesignWithOffsets;
use crate::errors::{ErrorKind, invalid_options, invariant, missing_prerequisite};
use crate::glm::{GlmInput, estimate_beta_prior_variance, fit_glms_with_beta_prior};
use crate::normalization::Offsets;

/// How the standard and expanded designs of a formula relate.
pub(super) struct ExpandedPriorPlan {
    /// Expanded column of each standard column.
    standard_to_expanded: Vec<usize>,
    /// Per factor: its expanded level columns, the expanded column of its
    /// reference level and its standard columns.
    factors: Vec<FactorColumns>,
    intercept: Option<usize>,
}

struct FactorColumns {
    expanded: Vec<usize>,
    reference: usize,
    standard: Vec<usize>,
}

impl ExpandedPriorPlan {
    pub(super) fn new(formula: &ExpandedFormulaDesignWithOffsets) -> Result<Self, DeseqError> {
        let design = &formula.design;
        let interactions = design.interactions.len()
            + design.factor_numeric_interactions.len()
            + design.numeric_interactions.len()
            + design.higher_order_interactions.len();
        if interactions > 0 {
            // DESeq2's designAndArgChecker rejects the same combination.
            return Err(DeseqError::UnsupportedFeature {
                feature: "the beta prior with interaction terms (use betaPrior = FALSE)".into(),
            });
        }
        let names = design.expanded_design.coefficient_names().ok_or_else(|| {
            invalid_options("beta prior: the expanded design has no column names")
        })?;
        let mut standard_to_expanded = Vec::with_capacity(design.coefficient_groups.len());
        for group in &design.coefficient_groups {
            match group.as_slice() {
                [column] => standard_to_expanded.push(*column),
                _ => {
                    return Err(DeseqError::UnsupportedFeature {
                        feature: "the beta prior with a standard column spanning several expanded columns".into(),
                    });
                }
            }
        }
        let mut factors = Vec::with_capacity(design.factor_levels.len());
        for levels in &design.factor_levels {
            let first = levels
                .first()
                .ok_or_else(|| invariant("a formula factor without levels"))?;
            // The factor's prefix is the one for which every level has a column.
            let prefix = names
                .iter()
                .filter_map(|name| name.strip_suffix(first.as_str()))
                .find(|prefix| {
                    levels
                        .iter()
                        .all(|level| names.iter().any(|name| *name == format!("{prefix}{level}")))
                })
                .ok_or_else(|| invariant("expanded columns for a formula factor"))?;
            let expanded: Vec<usize> = levels
                .iter()
                .filter_map(|level| {
                    names
                        .iter()
                        .position(|name| *name == format!("{prefix}{level}"))
                })
                .collect();
            let standard: Vec<usize> = standard_to_expanded
                .iter()
                .enumerate()
                .filter(|(_, column)| expanded.contains(column))
                .map(|(standard, _)| standard)
                .collect();
            let mut unmapped = expanded
                .iter()
                .copied()
                .filter(|column| !standard_to_expanded.contains(column));
            let reference = match (unmapped.next(), unmapped.next()) {
                (Some(reference), None) => reference,
                _ => return Err(invariant("one reference level per formula factor")),
            };
            factors.push(FactorColumns {
                expanded,
                reference,
                standard,
            });
        }
        let intercept = names.iter().position(|name| name == "Intercept");
        Ok(Self {
            standard_to_expanded,
            factors,
            intercept,
        })
    }

    /// The expanded contrast equal to a standard-design contrast.
    pub(super) fn expanded_contrast(&self, standard: &[f64], expanded_columns: usize) -> Vec<f64> {
        let mut contrast = vec![0.0; expanded_columns];
        for (column, weight) in standard.iter().copied().enumerate() {
            if weight == 0.0 {
                continue;
            }
            let target = self.standard_to_expanded[column];
            contrast[target] += weight;
            if Some(target) == self.intercept {
                // The standard intercept is the mean at every reference level.
                for factor in &self.factors {
                    contrast[factor.reference] += weight;
                }
            } else if let Some(factor) = self.factors.iter().find(|f| f.standard.contains(&column))
            {
                contrast[factor.reference] -= weight;
            }
        }
        contrast
    }

    /// `estimateBetaPriorVar` with `addAllContrasts` and
    /// `averagePriorsOverLevels`, from the MLE log2 betas of the standard
    /// design on the rows that are not all zero.
    pub(super) fn prior_variance(
        &self,
        mle_beta: &RowMajorMatrix<f64>,
        rows: &[usize],
        base_mean: &[f64],
        disp_fit: &[f64],
        expanded_names: &[String],
    ) -> Result<Vec<f64>, DeseqError> {
        let mut columns: Vec<Vec<f64>> = (0..mle_beta.n_cols())
            .map(|column| {
                rows.iter()
                    .map(|&row| mle_beta.get(row, column).copied().unwrap_or(f64::NAN))
                    .collect()
            })
            .collect();
        let mut names: Vec<String> = self
            .standard_to_expanded
            .iter()
            .map(|&column| expanded_names[column].clone())
            .collect();
        let mut factor_columns = Vec::with_capacity(self.factors.len());
        for factor in &self.factors {
            let first_contrast = columns.len();
            for (i, j) in pairwise_level_order(factor.standard.len()) {
                let difference = columns[factor.standard[i]]
                    .iter()
                    .zip(&columns[factor.standard[j]])
                    .map(|(a, b)| a - b)
                    .collect();
                columns.push(difference);
                names.push(String::new());
            }
            let mut all = factor.standard.clone();
            all.extend(first_contrast..columns.len());
            factor_columns.push(all);
        }
        let matrix = RowMajorMatrix::from_row_major(
            rows.len(),
            columns.len(),
            (0..rows.len())
                .flat_map(|row| columns.iter().map(move |column| column[row]))
                .collect(),
        )?;
        let base_mean: Vec<f64> = rows.iter().map(|&row| base_mean[row]).collect();
        let disp_fit: Vec<f64> = rows.iter().map(|&row| disp_fit[row]).collect();
        let prior = estimate_beta_prior_variance(
            &matrix,
            &base_mean,
            &disp_fit,
            Some(&names),
            Default::default(),
        )?;
        let mut expanded = vec![f64::NAN; expanded_names.len()];
        for (standard, &column) in self.standard_to_expanded.iter().enumerate() {
            expanded[column] = prior[standard];
        }
        for (factor, all) in self.factors.iter().zip(&factor_columns) {
            // R's mean(), as averagePriorsOverLevels uses it.
            let values: Vec<f64> = all.iter().map(|&column| prior[column]).collect();
            let mean = crate::design::numeric_transform::scalar_mean(&values);
            for &column in &factor.expanded {
                expanded[column] = mean;
            }
        }
        if let Some(bad) = expanded.iter().position(|v| !(v.is_finite() && *v > 0.0)) {
            return Err(invalid_options(format!(
                "beta prior for {} is not a positive finite variance",
                expanded_names[bad]
            )));
        }
        Ok(expanded)
    }
}

/// `addAllContrasts`' column pairs `(i, j)` for `M[, i] - M[, j]`.
fn pairwise_level_order(n: usize) -> Vec<(usize, usize)> {
    match n {
        0 | 1 => Vec::new(),
        2 => vec![(1, 0)],
        _ => {
            let is = (2..=n).flat_map(|k| k..=n).map(|i| i - 1);
            let js = (1..n)
                .zip((1..n).rev())
                .flat_map(|(j, count)| std::iter::repeat_n(j - 1, count));
            is.zip(js).collect()
        }
    }
}

impl DeseqBuilder {
    /// The Wald results of a beta-prior run, given the MAP-dispersion fit.
    pub(super) fn wald_beta_prior(
        &self,
        counts: &CountMatrix,
        formula: &ExpandedFormulaDesignWithOffsets,
        mut fit: DeseqFit,
        standard_contrast: &[f64],
    ) -> Result<BetaPriorRunOutput, DeseqError> {
        let design = &formula.design.standard_design;
        let expanded_design = &formula.design.expanded_design;
        let plan = ExpandedPriorPlan::new(formula)?;
        let dispersions = fit.dispersion.clone().ok_or_else(|| {
            missing_prerequisite(
                ErrorKind::Dispersion,
                "beta prior fitting",
                "MAP dispersions",
            )
        })?;
        let disp_fit = fit.disp_fit.clone().ok_or_else(|| {
            missing_prerequisite(
                ErrorKind::Dispersion,
                "beta prior fitting",
                "fitted dispersion trend values",
            )
        })?;
        let input = FixedDispersionGlmInput {
            counts,
            design,
            size_factors: &fit.size_factors,
            normalization_factors: fit.normalization_factors.as_ref(),
            observation_weights: fit.observation_weights.as_ref(),
            all_zero: &fit.all_zero,
            dispersions: &dispersions,
        };
        let FixedDispersionGlmOutput {
            glm_fit: mle,
            expanded_dispersions,
        } = self.fixed_dispersion_glm_components(input)?;

        let rows = nonzero_gene_indices(&fit.all_zero);
        let expanded_names = expanded_design
            .coefficient_names()
            .ok_or_else(|| invalid_options("beta prior: the expanded design has no column names"))?
            .to_vec();
        let prior =
            plan.prior_variance(&mle.beta, &rows, &fit.base_mean, &disp_fit, &expanded_names)?;

        // The MAP refit on the expanded design, on the rows that are not all zero.
        let expanded_fit = if rows.is_empty() {
            all_zero_glm_fit(counts, expanded_design)?
        } else {
            let compact = compact_counts(counts, &rows)?;
            let compact_dispersions: Vec<f64> = rows.iter().map(|&row| dispersions[row]).collect();
            let normalization_factors = fit
                .normalization_factors
                .as_ref()
                .map(|factors| compact_matrix_rows(factors, &rows))
                .transpose()?;
            let weights = fit
                .observation_weights
                .as_ref()
                .map(|weights| compact_matrix_rows(weights, &rows))
                .transpose()?;
            let offsets = match &normalization_factors {
                Some(factors) => Offsets::PerCell(factors),
                None => Offsets::PerSample(&fit.size_factors),
            };
            let compact_fit = fit_glms_with_beta_prior(
                GlmInput {
                    counts: &compact,
                    offsets,
                    dispersions: &compact_dispersions,
                    weights: weights.as_ref(),
                },
                expanded_design,
                &prior,
                self.irls_options.clone(),
            )?;
            expand_glm_fit(compact_fit, &fit.all_zero)?
        };

        let contrast = plan.expanded_contrast(standard_contrast, expanded_design.n_coefficients());
        let mut wald =
            wald_test_contrast_with_options(&expanded_fit, &contrast, &self.wald_test_options)?;
        mask_wald_degrees_of_freedom_for_all_zero_rows(&mut wald.wald, &fit.all_zero)?;
        let contrast_all_zero = contrast_all_zero_numeric(counts, design, standard_contrast)?;
        apply_contrast_all_zero_to_wald_contrast(&mut wald, &contrast_all_zero, &fit.all_zero)?;

        let normalized = match fit.normalization_factors.as_ref() {
            Some(factors) => normalized_counts_with_factors(counts, factors)?,
            None => normalized_counts(counts, &fit.size_factors)?,
        };
        // DESeq2 computes Cook's distances from the MLE fit's means and hat
        // diagonals on the standard design.
        let cooks =
            calculate_cooks_distance(counts, &normalized, &mle.mu, &mle.hat_diagonal, design)?;
        let mut results = build_wald_contrast_results_unadjusted(
            &fit.base_mean,
            &expanded_fit,
            &wald,
            counts.gene_names(),
            Some(&expanded_dispersions),
        )?;
        results.apply_wald_test_options(&self.wald_test_options);
        for (gene, all_zero) in fit.all_zero.iter().copied().enumerate() {
            results.rows[gene].max_cooks = cooks.max_cooks[gene];
            if all_zero {
                results.rows[gene].converged = None;
                results.rows[gene].max_cooks = None;
            }
        }
        let cooks_cutoff = resolve_cooks_cutoff(
            self.cooks_cutoff,
            design.n_samples(),
            design.n_coefficients(),
        )?;
        apply_cooks_cutoff_for_design(
            &mut results,
            cooks_cutoff,
            counts,
            Some(&cooks.cooks),
            design,
            self.model_frame.as_ref(),
            None,
        )?;
        apply_independent_filtering(&mut results, &self.independent_filtering_options)?;
        results.set_resolved_contrast_metadata(
            "contrast",
            "beta-prior expanded-model contrast",
            standard_contrast,
        );

        fit.dispersion = Some(expanded_dispersions);
        fit.cooks = Some(cooks.cooks);
        fit.max_cooks = Some(cooks.max_cooks);
        attach_glm_fit(&mut fit, mle);
        fit.wald = Some(wald.wald.clone());
        Ok(BetaPriorRunOutput {
            fit,
            results,
            expanded_design: expanded_design.clone(),
            beta_prior_variance: prior,
            expanded_fit,
            expanded_contrast: contrast,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::pairwise_level_order;

    #[test]
    fn pairwise_order_matches_add_all_contrasts() {
        // R: is <- do.call(c, sapply(seq_len(n-1)+1, function(k) seq(k, n)))
        //    js <- rep(seq_len(n-1), rev(seq_len(n-1)))
        assert_eq!(pairwise_level_order(2), vec![(1, 0)]);
        assert_eq!(
            pairwise_level_order(4),
            vec![(1, 0), (2, 0), (3, 0), (2, 1), (3, 1), (3, 2)]
        );
    }
}
