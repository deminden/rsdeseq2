use super::*;
use crate::errors::{
    ErrorKind, Requirement, invalid_formula, invalid_parameter, validate_matrix_shape,
};
use crate::matrix::try_alloc_matrix;
use crate::numeric::math::ProfileMath;
pub(super) fn formula_size_factor_offsets(
    size_factors: &[f64],
    offsets: &[f64],
    n_samples: usize,
) -> Result<Option<Vec<f64>>, DeseqError> {
    if !formula_offsets_are_active(offsets) {
        return Ok(None);
    }
    if size_factors.len() != n_samples {
        return Err(invalid_dimensions(
            "formula offset size factors",
            n_samples,
            size_factors.len(),
        ));
    }
    let mut factors = formula_offset_scales(offsets, n_samples)?;
    for (sample, (factor, &size_factor)) in factors.iter_mut().zip(size_factors).enumerate() {
        if !size_factor.is_finite() || size_factor <= 0.0 {
            return Err(invalid_parameter(
                ErrorKind::Options,
                format!("size factor at sample {sample}"),
                Requirement::FiniteAndPositive,
            ));
        }
        let product = size_factor * *factor;
        if !product.is_finite() || product <= 0.0 {
            return Err(invalid_formula(format!(
                "formula offset normalization factor at sample {sample} must be finite and positive"
            )));
        }
        *factor = product;
    }
    Ok(Some(factors))
}

pub(super) fn formula_normalization_factor_offsets(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    offsets: &[f64],
) -> Result<Option<RowMajorMatrix<f64>>, DeseqError> {
    if !formula_offsets_are_active(offsets) {
        return Ok(None);
    }
    validate_matrix_shape(
        "formula offset normalization factors",
        (counts.n_genes(), counts.n_samples()),
        (
            normalization_factors.n_rows(),
            normalization_factors.n_cols(),
        ),
    )?;
    let offset_scales = formula_offset_scales(offsets, counts.n_samples())?;
    let mut values = try_alloc_matrix(
        counts.n_genes(),
        counts.n_samples(),
        "formula offset normalization factors",
    )?;
    for gene in 0..normalization_factors.n_rows() {
        for (sample, value) in normalization_factors.row(gene)?.iter().copied().enumerate() {
            if !value.is_finite() || value <= 0.0 {
                return Err(invalid_parameter(
                    ErrorKind::Options,
                    format!("normalization factor at gene {gene}, sample {sample}"),
                    Requirement::FiniteAndPositive,
                ));
            }
            let factor = value * offset_scales[sample];
            if !factor.is_finite() || factor <= 0.0 {
                return Err(invalid_formula(format!(
                    "formula offset normalization factor at gene {gene}, sample {sample} must be finite and positive"
                )));
            }
            values.push(factor);
        }
    }
    RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values).map(Some)
}

fn formula_offsets_are_active(offsets: &[f64]) -> bool {
    offsets.iter().any(|value| *value != 0.0)
}

pub(super) fn attach_cooks_to_results(
    results: &mut DeseqResults,
    max_cooks: &[Option<f64>],
) -> Result<(), DeseqError> {
    if max_cooks.len() != results.rows.len() {
        return Err(invalid_dimensions(
            "Cook's result rows",
            results.rows.len(),
            max_cooks.len(),
        ));
    }
    for (row, max_cook) in results.rows.iter_mut().zip(max_cooks.iter().copied()) {
        row.max_cooks = max_cook;
        row.cooks_outlier = None;
    }
    Ok(())
}

pub(super) fn merge_beta_prior_replacement_results(
    original_results: &DeseqResults,
    refit_results: Option<&DeseqResults>,
    refit_plan: &CooksRefitPlan,
) -> Result<DeseqResults, DeseqError> {
    if original_results.rows.len() != refit_plan.replacement.replace.len() {
        return Err(invalid_dimensions(
            "beta-prior replacement result rows",
            refit_plan.replacement.replace.len(),
            original_results.rows.len(),
        ));
    }
    if let Some(refit_results) = refit_results
        && refit_results.rows.len() != original_results.rows.len()
    {
        return Err(invalid_dimensions(
            "beta-prior replacement refit result rows",
            original_results.rows.len(),
            refit_results.rows.len(),
        ));
    }
    if refit_plan.replaced_base_mean.len() != original_results.rows.len() {
        return Err(invalid_dimensions(
            "beta-prior replacement baseMean rows",
            original_results.rows.len(),
            refit_plan.replaced_base_mean.len(),
        ));
    }
    if refit_plan.post_refit_max_cooks.len() != original_results.rows.len() {
        return Err(invalid_dimensions(
            "beta-prior replacement maxCooks rows",
            original_results.rows.len(),
            refit_plan.post_refit_max_cooks.len(),
        ));
    }

    if !refit_plan.has_replaceable_samples() {
        return Ok(original_results.clone());
    }

    // Keep the original result-table cardinality and patch in refit rows only;
    // this mirrors DESeq2's replacement-count refit rather than rebuilding the
    // entire result table from replacement counts.
    let has_refit = refit_plan.n_refit > 0 && refit_plan.requires_refit();
    let mut merged = original_results.clone();
    for (gene, row) in merged.rows.iter_mut().enumerate() {
        row.base_mean = refit_plan.replaced_base_mean[gene];
        if has_refit {
            row.max_cooks = refit_plan.post_refit_max_cooks[gene];
            row.cooks_outlier = None;
            row.filtered = None;
        }
    }

    if let Some(refit_results) = refit_results {
        for gene in refit_plan.refit_rows.iter().copied() {
            merged.rows[gene] = refit_results.rows[gene].clone();
            merged.rows[gene].base_mean = refit_plan.replaced_base_mean[gene];
            merged.rows[gene].max_cooks = refit_plan.post_refit_max_cooks[gene];
            merged.rows[gene].cooks_outlier = None;
            merged.rows[gene].filtered = None;
        }
    }

    for gene in refit_plan.new_all_zero_rows.iter().copied() {
        clear_replacement_all_zero_result(&mut merged.rows[gene]);
        merged.rows[gene].base_mean = refit_plan.replaced_base_mean[gene];
        if has_refit {
            merged.rows[gene].max_cooks = refit_plan.post_refit_max_cooks[gene];
        }
    }

    merged.independent_filtering = None;
    Ok(merged)
}

pub(super) fn clear_replacement_all_zero_result(row: &mut DeseqResultRow) {
    row.log2_fold_change = None;
    row.lfc_se = None;
    row.stat = None;
    row.pvalue = None;
    row.padj = None;
    row.dispersion = None;
    row.converged = None;
    row.cooks_outlier = None;
    row.filtered = None;
}

fn formula_offset_scales(offsets: &[f64], n_samples: usize) -> Result<Vec<f64>, DeseqError> {
    if offsets.len() != n_samples {
        return Err(invalid_dimensions(
            "formula offsets",
            n_samples,
            offsets.len(),
        ));
    }
    let mut scales = try_alloc_matrix(n_samples, 1, "formula offset scales")?;
    for (sample, &offset) in offsets.iter().enumerate() {
        let scale = offset.exp_p();
        if !scale.is_finite() || scale <= 0.0 {
            return Err(invalid_formula(format!(
                "formula offset scale at sample {sample} must be finite and positive"
            )));
        }
        scales.push(scale);
    }
    Ok(scales)
}
