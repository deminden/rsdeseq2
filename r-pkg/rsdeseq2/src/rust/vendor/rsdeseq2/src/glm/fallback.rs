use crate::errors::{DeseqError, invalid_dimensions};
use crate::matrix::RowMajorMatrix;

/// DESeq2-shaped row routing for bounded optim fallback fitting.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct OptimFallbackRows {
    /// Rows that should be sent through fallback optimization.
    pub rows: Vec<usize>,
    /// Rows whose IRLS beta estimates are not finite.
    pub unstable_rows: Vec<usize>,
    /// Rows whose IRLS coefficient variances are not positive finite values.
    pub non_positive_variance_rows: Vec<usize>,
    /// Rows whose IRLS convergence flag is false.
    pub non_converged_rows: Vec<usize>,
}

/// Identify rows that DESeq2 would route to its `optim` backup path.
///
/// When `use_optim` is false, only unstable rows or rows with invalid
/// coefficient variances are returned. When `force_optim` is true, every row is
/// returned. This helper isolates the routing rule before the fixed-dispersion
/// IRLS path applies bounded fallback refits.
pub fn optim_fallback_rows(
    beta_converged: &[bool],
    beta: &RowMajorMatrix<f64>,
    beta_covariance: &RowMajorMatrix<f64>,
    use_optim: bool,
    force_optim: bool,
) -> Result<OptimFallbackRows, DeseqError> {
    validate_fallback_inputs(beta_converged, beta, beta_covariance)?;
    let mut output = OptimFallbackRows::default();
    let p = beta.n_cols();

    for (gene, converged) in beta_converged.iter().copied().enumerate() {
        let stable = row_values_finite(beta.row(gene)?);
        let variance_positive = covariance_diagonal_positive(beta_covariance.row(gene)?, p);
        if !stable {
            output.unstable_rows.push(gene);
        }
        if !variance_positive {
            output.non_positive_variance_rows.push(gene);
        }
        if !converged {
            output.non_converged_rows.push(gene);
        }
        if force_optim || !stable || !variance_positive || (use_optim && !converged) {
            output.rows.push(gene);
        }
    }

    Ok(output)
}

fn validate_fallback_inputs(
    beta_converged: &[bool],
    beta: &RowMajorMatrix<f64>,
    beta_covariance: &RowMajorMatrix<f64>,
) -> Result<(), DeseqError> {
    if beta_converged.len() != beta.n_rows() {
        return Err(invalid_dimensions(
            "optim fallback convergence rows",
            beta.n_rows(),
            beta_converged.len(),
        ));
    }
    if beta_covariance.n_rows() != beta.n_rows() {
        return Err(invalid_dimensions(
            "optim fallback covariance rows",
            beta.n_rows(),
            beta_covariance.n_rows(),
        ));
    }
    let expected_covariance_cols = beta.n_cols() * beta.n_cols();
    if beta_covariance.n_cols() != expected_covariance_cols {
        return Err(invalid_dimensions(
            "optim fallback covariance columns",
            expected_covariance_cols,
            beta_covariance.n_cols(),
        ));
    }
    Ok(())
}

fn row_values_finite(values: &[f64]) -> bool {
    values.iter().copied().all(f64::is_finite)
}

fn covariance_diagonal_positive(values: &[f64], p: usize) -> bool {
    (0..p).all(|idx| {
        let value = values[idx * p + idx];
        value.is_finite() && value > 0.0
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fallback_rows_follow_deseq2_use_optim_switch() {
        let beta = RowMajorMatrix::from_row_major(
            4,
            2,
            vec![
                1.0,
                2.0, //
                1.0,
                f64::NAN, //
                1.0,
                2.0, //
                1.0,
                2.0,
            ],
        )
        .unwrap();
        let covariance = RowMajorMatrix::from_row_major(
            4,
            4,
            vec![
                1.0, 0.0, 0.0, 2.0, //
                1.0, 0.0, 0.0, 2.0, //
                1.0, 0.0, 0.0, 0.0, //
                1.0, 0.0, 0.0, 2.0,
            ],
        )
        .unwrap();
        let beta_converged = vec![true, true, true, false];

        let no_optim =
            optim_fallback_rows(&beta_converged, &beta, &covariance, false, false).unwrap();
        assert_eq!(no_optim.rows, vec![1, 2]);
        assert_eq!(no_optim.unstable_rows, vec![1]);
        assert_eq!(no_optim.non_positive_variance_rows, vec![2]);
        assert_eq!(no_optim.non_converged_rows, vec![3]);

        let use_optim =
            optim_fallback_rows(&beta_converged, &beta, &covariance, true, false).unwrap();
        assert_eq!(use_optim.rows, vec![1, 2, 3]);
    }

    #[test]
    fn force_optim_routes_every_row() {
        let beta = RowMajorMatrix::from_row_major(2, 1, vec![1.0, 2.0]).unwrap();
        let covariance = RowMajorMatrix::from_row_major(2, 1, vec![1.0, 1.0]).unwrap();
        let rows = optim_fallback_rows(&[true, true], &beta, &covariance, false, true).unwrap();
        assert_eq!(rows.rows, vec![0, 1]);
    }

    #[test]
    fn fallback_rows_validate_dimensions() {
        let beta = RowMajorMatrix::from_row_major(1, 2, vec![1.0, 2.0]).unwrap();
        let covariance = RowMajorMatrix::from_row_major(1, 3, vec![1.0, 0.0, 1.0]).unwrap();
        assert!(optim_fallback_rows(&[true], &beta, &covariance, false, false).is_err());
        let good_covariance =
            RowMajorMatrix::from_row_major(1, 4, vec![1.0, 0.0, 0.0, 1.0]).unwrap();
        assert!(optim_fallback_rows(&[], &beta, &good_covariance, false, false).is_err());
    }
}
