use rsdeseq2::prelude::*;

fn matrix(rows: usize, columns: usize, values: &[f64]) -> MatrixView<'_, f64> {
    MatrixView::new(rows, columns, values, MatrixOrder::RowMajor).unwrap()
}

fn column_major(rows: usize, columns: usize, values: &[f64]) -> Vec<f64> {
    (0..columns)
        .flat_map(|column| (0..rows).map(move |row| values[row * columns + column]))
        .collect()
}

fn bits(values: &[Option<f64>]) -> Vec<Option<u64>> {
    values.iter().map(|v| v.map(f64::to_bits)).collect()
}

#[test]
fn layouts_share_logical_covariance_order_and_every_inference_bit() {
    // Deliberately asymmetric entries expose a swapped covariance index;
    // unequal gene/coefficient dimensions expose a transposed gene axis.
    let beta = [1., 3., -2., 4., 5., 2.];
    let covariance = [4., 0.5, 0.25, 9., 1., -0.5, 0.25, 4., 9., 1., 0.5, 16.];
    let beta_cm = column_major(3, 2, &beta);
    let mut cov_cm = vec![0.; covariance.len()];
    for gene in 0..3 {
        for row in 0..2 {
            for column in 0..2 {
                cov_cm[gene + 3 * (row + column * 2)] = covariance[gene * 4 + row * 2 + column];
            }
        }
    }
    let row = InferenceEstimates::new(
        matrix(3, 2, &beta),
        None,
        Some(CovarianceView::new(2, matrix(3, 4, &covariance)).unwrap()),
    )
    .unwrap();
    let column = InferenceEstimates::new(
        MatrixView::new(3, 2, &beta_cm, MatrixOrder::ColumnMajor).unwrap(),
        None,
        Some(
            CovarianceView::new(
                2,
                MatrixView::new(3, 4, &cov_cm, MatrixOrder::ColumnMajor).unwrap(),
            )
            .unwrap(),
        ),
    )
    .unwrap();
    for contrast in [[-1., 1.], [0.25, -2.], [0., -1.]] {
        let a = row.project_contrast(&contrast).unwrap();
        let b = column.project_contrast(&contrast).unwrap();
        assert_eq!(bits(&a.log2_fold_change), bits(&b.log2_fold_change));
        assert_eq!(bits(&a.lfc_se), bits(&b.lfc_se));
        for options in [
            WaldTestOptions::normal(),
            WaldTestOptions::t_per_gene_degrees_of_freedom(vec![2., f64::NAN, 10.])
                .with_lfc_threshold(0.75, WaldAlternative::Greater),
        ] {
            let a = wald_test_contrast_from_estimates(&row, &contrast, &options).unwrap();
            let b = wald_test_contrast_from_estimates(&column, &contrast, &options).unwrap();
            assert_eq!(bits(&a.log2_fold_change), bits(&b.log2_fold_change));
            assert_eq!(bits(&a.lfc_se), bits(&b.lfc_se));
            assert_eq!(bits(&a.wald.stat), bits(&b.wald.stat));
            assert_eq!(bits(&a.wald.pvalue), bits(&b.wald.pvalue));
        }
    }
    let projected = row.project_contrast(&[-1., 1.]).unwrap();
    assert_eq!(
        projected.log2_fold_change,
        vec![Some(2.), Some(6.), Some(-3.)]
    );
    for (se, variance) in projected.lfc_se.iter().zip([12.25_f64, 5.25, 23.5]) {
        assert!((se.unwrap() - variance.sqrt()).abs() < 2e-15);
    }
}

#[test]
fn marginal_se_projection_and_wald_invalid_row_policy_are_shared() {
    let beta = [1., 2., 1., 3., 1., 4., 1., 5.];
    let se = [1., 0.5, 1., 0., 1., f64::NAN, 1., 2.];
    // Covariance is deliberately unusable: a single coefficient must use the
    // same supplied marginal SE as a named coefficient in either interface.
    let cov = [f64::NAN; 16];
    let estimates = InferenceEstimates::new(
        matrix(4, 2, &beta),
        Some(matrix(4, 2, &se)),
        Some(CovarianceView::new(2, matrix(4, 4, &cov)).unwrap()),
    )
    .unwrap();
    let projected = estimates.project_contrast(&[0., -2.]).unwrap();
    assert_eq!(
        projected.log2_fold_change,
        vec![Some(-4.), Some(-6.), Some(-8.), Some(-10.)]
    );
    assert_eq!(projected.lfc_se, vec![Some(1.), Some(0.), None, Some(4.)]);
    let options = WaldTestOptions::t_per_gene_degrees_of_freedom(vec![10., 10., 10., f64::NAN]);
    let wald = wald_test_contrast_from_estimates(&estimates, &[0., -2.], &options).unwrap();
    assert_eq!(
        wald.log2_fold_change,
        vec![Some(-4.), None, None, Some(-10.)]
    );
    assert_eq!(wald.lfc_se, vec![Some(1.), None, None, Some(4.)]);
    assert_eq!(wald.wald.stat, vec![Some(-4.), None, None, Some(-2.5)]);
    assert!(wald.wald.pvalue[0].is_some());
    assert_eq!(&wald.wald.pvalue[1..], &[None, None, None]);
    // Reuse the borrowed view for an actual covariance request: missing
    // covariance yields missing rows, rather than substituting diagonal SEs.
    assert_eq!(
        estimates.project_contrast(&[1., -1.]).unwrap().lfc_se,
        vec![None; 4]
    );
}

#[test]
fn borrowed_estimates_reject_shape_aliases_and_missing_context() {
    assert!(MatrixView::new(usize::MAX, 2, &[] as &[f64], MatrixOrder::RowMajor).is_err());
    assert!(MatrixView::new(1, 0, &[] as &[f64], MatrixOrder::ColumnMajor).is_err());
    assert!(MatrixView::new(2, 2, &[1., 2., 3.], MatrixOrder::RowMajor).is_err());
    let beta = [1.; 6];
    let wrong_se =
        InferenceEstimates::new(matrix(2, 3, &beta), Some(matrix(3, 2, &beta)), None).unwrap_err();
    assert!(wrong_se.to_string().contains("betaSE rows"));
    let covariance = [1.; 8];
    let wrong_cov = InferenceEstimates::new(
        matrix(2, 3, &beta),
        None,
        Some(CovarianceView::new(2, matrix(2, 4, &covariance)).unwrap()),
    )
    .unwrap_err();
    assert!(wrong_cov.to_string().contains("covariance coefficients"));
    let view =
        InferenceEstimates::new(matrix(2, 3, &beta), Some(matrix(2, 3, &beta)), None).unwrap();
    assert!(view.project_contrast(&[1., 0.]).is_err());
    assert!(view.project_contrast(&[0.; 3]).is_err());
    assert!(view.project_contrast(&[1., f64::INFINITY, 0.]).is_err());
    assert!(
        view.project_contrast(&[1., -1., 0.])
            .unwrap_err()
            .to_string()
            .contains("betaCovariance")
    );
    assert!(
        wald_test_contrast_from_estimates(
            &view,
            &[0., 1., 0.],
            &WaldTestOptions::t_residual_degrees_of_freedom()
        )
        .is_err()
    );
    assert!(matrix(2, 3, &beta).get(2, 0).is_err());
    assert!(matrix(2, 3, &beta).get(0, 3).is_err());
    // Borrowing does not materialize another matrix.
    assert_eq!(matrix(2, 3, &beta).as_slice().as_ptr(), beta.as_ptr());
}
