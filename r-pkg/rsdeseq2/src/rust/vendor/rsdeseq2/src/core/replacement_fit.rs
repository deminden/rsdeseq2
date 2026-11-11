//! Consume replacement outputs into one reusable, row-aligned fitted state.
use super::{CooksReplacementTestOutput, DeseqFit};
use crate::{
    cooks::CooksRefitPlan, errors::DeseqError, matrix::RowMajorMatrix, results::DeseqResults,
};
fn invalid(reason: impl Into<String>) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: reason.into(),
    }
}

fn merge_vec<T: Clone>(
    original: &mut [T],
    refit: &[T],
    rows: &[usize],
    field: &str,
) -> Result<(), DeseqError> {
    if rows.is_empty() {
        return Ok(());
    }
    if original.len() != refit.len() || rows.iter().any(|&row| row >= original.len()) {
        return Err(invalid(format!(
            "replacement refit {field} dimensions do not match"
        )));
    }
    for &row in rows {
        original[row] = refit[row].clone();
    }
    Ok(())
}

fn merge_optional_vec<T: Clone>(
    original: &mut Option<Vec<T>>,
    refit: &Option<Vec<T>>,
    rows: &[usize],
    field: &str,
) -> Result<(), DeseqError> {
    if rows.is_empty() {
        return Ok(());
    }
    if original.is_none() && refit.is_none() {
        return Ok(());
    }
    let target = original
        .as_mut()
        .ok_or_else(|| invalid(format!("original fit is missing {field}")))?;
    let source = refit
        .as_ref()
        .ok_or_else(|| invalid(format!("replacement refit is missing {field}")))?;
    merge_vec(target, source, rows, field)
}

fn merge_option_matrix(
    original: &mut Option<RowMajorMatrix<f64>>,
    refit: &Option<RowMajorMatrix<f64>>,
    rows: &[usize],
    field: &str,
) -> Result<(), DeseqError> {
    if rows.is_empty() {
        return Ok(());
    }
    if original.is_none() && refit.is_none() {
        return Ok(());
    }
    let target = original
        .as_ref()
        .ok_or_else(|| invalid(format!("original fit is missing {field}")))?;
    let source = refit
        .as_ref()
        .ok_or_else(|| invalid(format!("replacement refit is missing {field}")))?;
    if target.n_rows() != source.n_rows()
        || target.n_cols() != source.n_cols()
        || rows.iter().any(|&row| row >= target.n_rows())
    {
        return Err(invalid(format!(
            "replacement refit {field} dimensions do not match"
        )));
    }
    let n_rows = target.n_rows();
    let n_cols = target.n_cols();
    let owned = original
        .take()
        .ok_or_else(|| invalid(format!("original fit is missing {field}")))?;
    let mut values = owned.into_values();
    for &row in rows {
        let start = row * n_cols;
        values[start..start + n_cols].copy_from_slice(&source.as_slice()[start..start + n_cols]);
    }
    *original = Some(RowMajorMatrix::from_row_major(n_rows, n_cols, values)?);
    Ok(())
}

fn merge_refit_fit_rows(
    original: &mut DeseqFit,
    refit: &DeseqFit,
    rows: &[usize],
) -> Result<(), DeseqError> {
    macro_rules! merge_required_vec {
        ($field:ident) => {
            merge_vec(
                &mut original.$field,
                &refit.$field,
                rows,
                stringify!($field),
            )?;
        };
    }
    macro_rules! merge_optional_vec {
        ($field:ident) => {
            merge_optional_vec(
                &mut original.$field,
                &refit.$field,
                rows,
                stringify!($field),
            )?;
        };
    }
    macro_rules! merge_matrix {
        ($field:ident) => {
            merge_option_matrix(
                &mut original.$field,
                &refit.$field,
                rows,
                stringify!($field),
            )?;
        };
    }
    merge_required_vec!(base_mean);
    merge_required_vec!(base_var);
    merge_required_vec!(all_zero);
    merge_optional_vec!(weights_fail);
    merge_optional_vec!(disp_gene_est);
    merge_optional_vec!(disp_gene_iter);
    merge_optional_vec!(disp_fit);
    merge_optional_vec!(disp_map);
    merge_optional_vec!(dispersion);
    merge_optional_vec!(disp_iter);
    merge_optional_vec!(disp_outlier);
    merge_optional_vec!(dispersion_converged);
    merge_matrix!(beta);
    merge_matrix!(beta_se);
    merge_matrix!(beta_optim_start);
    merge_matrix!(beta_covariance);
    merge_optional_vec!(beta_converged);
    merge_optional_vec!(beta_iter);
    merge_optional_vec!(beta_optim_iter);
    merge_optional_vec!(beta_optim_start_objective);
    merge_optional_vec!(beta_optim_objective);
    merge_optional_vec!(beta_optim_gradient_norm);
    merge_optional_vec!(log_like);
    merge_optional_vec!(full_deviance);
    merge_matrix!(mu);
    merge_matrix!(hat_diagonal);
    merge_optional_vec!(reduced_log_like);
    merge_optional_vec!(reduced_beta_converged);
    merge_optional_vec!(reduced_beta_iter);
    merge_matrix!(reduced_mu);
    merge_matrix!(reduced_hat_diagonal);
    if let (Some(target), Some(source)) = (&mut original.wald, &refit.wald) {
        merge_vec(&mut target.stat, &source.stat, rows, "Wald statistics")?;
        merge_vec(&mut target.pvalue, &source.pvalue, rows, "Wald p-values")?;
        if target.degrees_of_freedom.is_some() || source.degrees_of_freedom.is_some() {
            merge_optional_vec(
                &mut target.degrees_of_freedom,
                &source.degrees_of_freedom,
                rows,
                "Wald df",
            )?;
        }
    }
    if let (Some(target), Some(source)) = (&mut original.lrt, &refit.lrt) {
        merge_vec(
            &mut target.deviance,
            &source.deviance,
            rows,
            "LRT statistics",
        )?;
        merge_vec(&mut target.pvalue, &source.pvalue, rows, "LRT p-values")?;
        merge_vec(
            &mut target.reduced_converged,
            &source.reduced_converged,
            rows,
            "LRT reduced convergence",
        )?;
    }
    // Keep the original Cook's distances: DESeq2 records those before replacement.
    Ok(())
}

impl CooksReplacementTestOutput {
    /// Consume an output into a row-aligned fit, final results and replacement plan.
    ///
    /// Only selected refit rows are updated. Global trend/normalization and the
    /// original Cook's distances/maxima stay attached; other per-gene fitting
    /// diagnostics follow their refitted rows. Newly all-zero rows receive the
    /// established zero base metadata without fabricating a replacement fit.
    /// Existing separated original/refit APIs remain available before consuming.
    pub fn into_merged_fit(self) -> Result<(DeseqFit, DeseqResults, CooksRefitPlan), DeseqError> {
        let (mut fit, refit, plan, results) = match self {
            Self::Wald(o) => (o.original_fit, o.refit_fit, o.refit_plan, o.results),
            Self::Lrt(o) => (o.original_fit, o.refit_fit, o.refit_plan, o.results),
        };
        // Original refitWithoutOutliers recomputes base metadata for every row
        // after actual replacement, resetting earlier weight-failure skip flags
        // to raw count zeros. Selected refit rows may acquire skip flags again.
        if plan.has_replaceable_samples() && plan.n_refit > 0 {
            if plan.replaced_all_zero.len() != fit.all_zero.len() {
                return Err(invalid("replacement all-zero dimensions do not match"));
            }
            fit.all_zero.copy_from_slice(&plan.replaced_all_zero);
        }
        if let Some(refit) = &refit {
            if fit.design != refit.design || fit.reduced_design != refit.reduced_design {
                return Err(invalid("replacement refit designs do not match"));
            }
            merge_refit_fit_rows(&mut fit, refit, &plan.refit_rows)?;
        } else if plan.has_replaceable_samples() && !plan.refit_rows.is_empty() {
            return Err(invalid("replacement rows require a refitted model"));
        }
        if plan.new_all_zero_rows.iter().any(|&i| {
            i >= fit.base_mean.len() || i >= fit.base_var.len() || i >= fit.all_zero.len()
        }) {
            return Err(invalid("replacement new all-zero row is out of range"));
        }
        for &i in &plan.new_all_zero_rows {
            fit.base_mean[i] = 0.0;
            fit.base_var[i] = 0.0;
            fit.all_zero[i] = true;
        }
        Ok((fit, results, plan))
    }
}
