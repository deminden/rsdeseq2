// Merge fitted rows and diagnostics after Cook's replacement.

fn merge_vec<T: Clone>(
    original: &mut [T],
    refit: &[T],
    rows: &[usize],
    field: &str,
) -> Result<(), String> {
    if rows.is_empty() {
        return Ok(());
    }
    if original.len() != refit.len() || rows.iter().any(|&row| row >= original.len()) {
        return Err(format!("replacement refit {field} dimensions do not match"));
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
) -> Result<(), String> {
    if rows.is_empty() {
        return Ok(());
    }
    let target = original
        .as_mut()
        .ok_or_else(|| format!("original fit is missing {field}"))?;
    let source = refit
        .as_ref()
        .ok_or_else(|| format!("replacement refit is missing {field}"))?;
    merge_vec(target, source, rows, field)
}

fn merge_option_matrix(
    original: &mut Option<RowMajorMatrix<f64>>,
    refit: &Option<RowMajorMatrix<f64>>,
    rows: &[usize],
    field: &str,
) -> Result<(), String> {
    if rows.is_empty() {
        return Ok(());
    }
    let target = original
        .as_ref()
        .ok_or_else(|| format!("original fit is missing {field}"))?;
    let source = refit
        .as_ref()
        .ok_or_else(|| format!("replacement refit is missing {field}"))?;
    if target.n_rows() != source.n_rows()
        || target.n_cols() != source.n_cols()
        || rows.iter().any(|&row| row >= target.n_rows())
    {
        return Err(format!("replacement refit {field} dimensions do not match"));
    }
    let mut values = target.as_slice().to_vec();
    for &row in rows {
        let start = row * target.n_cols();
        values[start..start + target.n_cols()]
            .clone_from_slice(&source.as_slice()[start..start + target.n_cols()]);
    }
    *original = Some(
        RowMajorMatrix::from_row_major(target.n_rows(), target.n_cols(), values)
            .map_err(|error| error.to_string())?,
    );
    Ok(())
}

fn merge_refit_fit_rows(
    original: &mut rsdeseq2::prelude::DeseqFit,
    refit: &rsdeseq2::prelude::DeseqFit,
    rows: &[usize],
) -> Result<(), String> {
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
    merge_matrix!(beta_covariance);
    merge_optional_vec!(beta_converged);
    merge_optional_vec!(beta_iter);
    merge_optional_vec!(log_like);
    merge_optional_vec!(full_deviance);
    merge_matrix!(mu);
    // Keep the original Cook's distances: DESeq2 records those before replacement.
    Ok(())
}
