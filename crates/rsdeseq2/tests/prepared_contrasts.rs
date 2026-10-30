use rsdeseq2::prelude::*;

fn matrix<T>(rows: usize, columns: usize, values: &[T]) -> MatrixView<'_, T> {
    MatrixView::new(rows, columns, values, MatrixOrder::RowMajor).unwrap()
}

#[test]
fn prepared_contrasts_stream_independent_wald_raw_and_zero_requests() {
    let beta = [1., 2., 2., 4., 0., -1.];
    let se = [0.1, 0.5, 0.1, 0.5, 0.1, 0.5];
    let counts = [0, 0, 0, 0, 5, 6, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0];
    let design = [1., 0., 1., 0., 1., 1., 1., 1., 1., 0., 1., 0.];
    let estimates =
        InferenceEstimates::new(matrix(3, 2, &beta), Some(matrix(3, 2, &se)), None).unwrap();
    let context = PreparedContrastView::new(
        estimates,
        Some(ContrastCounts::from_unsigned(matrix(3, 6, &counts))),
        Some(matrix(6, 2, &design)),
    )
    .unwrap();
    let levels = ["A", "A", "B", "B", "C", "C"];
    let zero = ContrastZeroRule::Factor {
        levels: &levels,
        numerator: "B",
        denominator: "A",
    };
    let options = WaldTestOptions::normal();
    let factor = context.evaluate(&[0., 1.], Some(&options), zero).unwrap();
    assert_eq!(
        factor.projection.log2_fold_change,
        vec![Some(0.), Some(4.), Some(-1.)]
    );
    assert_eq!(factor.projection.lfc_se, vec![Some(0.5); 3]);
    assert_eq!(
        factor.wald.unwrap().stat,
        vec![Some(0.), Some(8.), Some(-2.)]
    );
    let numeric = context
        .evaluate(&[0., 1.], Some(&options), ContrastZeroRule::Numeric)
        .unwrap();
    assert_eq!(
        numeric.projection.log2_fold_change,
        vec![Some(2.), Some(4.), Some(-1.)]
    );
    let raw = context.evaluate(&[0., 1.], None, zero).unwrap();
    assert_eq!(
        raw.projection.log2_fold_change,
        vec![Some(0.), Some(4.), Some(-1.)]
    );
    assert!(raw.wald.is_none());
    let residual = context
        .evaluate(
            &[0., 1.],
            Some(&WaldTestOptions::t_residual_degrees_of_freedom()),
            ContrastZeroRule::None,
        )
        .unwrap();
    let explicit = context
        .evaluate(
            &[0., 1.],
            Some(&WaldTestOptions::t_degrees_of_freedom(4.)),
            ContrastZeroRule::None,
        )
        .unwrap();
    assert_eq!(residual, explicit);
    // Reusing the view must preserve each request's sign and scale.
    for weight in [1., -1., 2.] {
        let contrast = [0., weight];
        let actual = context
            .evaluate(&contrast, Some(&options), ContrastZeroRule::None)
            .unwrap();
        let expected = wald_test_contrast_from_estimates(&estimates, &contrast, &options).unwrap();
        assert_eq!(
            actual.projection.log2_fold_change,
            expected.log2_fold_change
        );
        assert_eq!(actual.wald.unwrap(), expected.wald);
    }
}

#[test]
fn prepared_shapes_and_primitive_context_are_validated_without_fabricated_fit_state() {
    let values = [1.; 6];
    let estimates =
        InferenceEstimates::new(matrix(3, 2, &values), Some(matrix(3, 2, &values)), None).unwrap();
    let wrong_genes = ContrastCounts::from_unsigned(matrix(2, 3, &[0; 6]));
    assert!(PreparedContrastView::new(estimates, Some(wrong_genes), None).is_err());
    assert!(PreparedContrastView::new(estimates, None, Some(matrix(2, 3, &values))).is_err());
    let context = PreparedContrastView::new(estimates, None, None).unwrap();
    assert!(
        context
            .evaluate(
                &[0., 1.],
                Some(&WaldTestOptions::t_residual_degrees_of_freedom()),
                ContrastZeroRule::None
            )
            .is_err()
    );
    let counts = [0; 12];
    let covariance_values = [1., 0., 0., 1., 1., 0., 0., 1., 1., 0., 0., 1.];
    let estimates = InferenceEstimates::new(
        matrix(3, 2, &values),
        Some(matrix(3, 2, &values)),
        Some(CovarianceView::new(2, matrix(3, 4, &covariance_values)).unwrap()),
    )
    .unwrap();
    assert!(
        PreparedContrastView::new(
            estimates,
            Some(ContrastCounts::from_unsigned(matrix(3, 4, &counts))),
            Some(matrix(3, 2, &values)),
        )
        .is_err()
    );
    let context = PreparedContrastView::new(
        estimates,
        Some(ContrastCounts::from_unsigned(matrix(3, 4, &counts))),
        None,
    )
    .unwrap();
    assert!(
        context
            .evaluate(
                &[-1., 1.],
                Some(&WaldTestOptions::normal()),
                ContrastZeroRule::Numeric
            )
            .unwrap_err()
            .to_string()
            .contains("modelMatrix")
    );
    assert!(
        context
            .evaluate(
                &[0., 1.],
                Some(&WaldTestOptions::normal()),
                ContrastZeroRule::Numeric
            )
            .is_ok()
    );
}
