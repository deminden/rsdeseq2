use rsdeseq2::prelude::*;

fn matrix<T>(rows: usize, columns: usize, values: &[T], order: MatrixOrder) -> MatrixView<'_, T> {
    MatrixView::new(rows, columns, values, order).unwrap()
}

#[test]
fn borrowed_mask_retains_low_count_threshold_missing_and_first_tie_rules_in_both_layouts() {
    let counts = [
        1, 5, 6, 7, 0, 1, 5, 6, 0, 0, 1, 9, 5, 6, 7, 1, 5, 6, 7, 0, 9, 5, 6, 7, 0, 1, 2, 3, 4, 5,
    ];
    let cooks = [
        10.,
        0.1,
        0.2,
        0.3,
        0.4,
        10.,
        0.1,
        0.2,
        0.3,
        0.4,
        10.,
        10.,
        0.2,
        0.3,
        0.4,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
        f64::NAN,
        5.,
        0.1,
        0.2,
        0.3,
        0.4,
        10.,
        0.1,
        0.2,
        0.3,
        0.4,
    ];
    let pvalues = [
        Some(0.01),
        Some(0.02),
        Some(0.03),
        Some(0.04),
        Some(0.05),
        Some(0.06),
    ];
    let maximum = [Some(10.), Some(10.), Some(10.), Some(10.), Some(5.), None];
    let native = LowCountHeuristic::new(
        CountView::from_unsigned(matrix(6, 5, &counts, MatrixOrder::RowMajor)),
        matrix(6, 5, &cooks, MatrixOrder::RowMajor),
    )
    .unwrap();
    let mut count_columns = Vec::new();
    let mut cook_columns = Vec::new();
    for sample in 0..5 {
        for gene in 0..6 {
            count_columns.push(f64::from(counts[gene * 5 + sample]));
            cook_columns.push(cooks[gene * 5 + sample]);
        }
    }
    let foreign = LowCountHeuristic::new(
        CountView::from_doubles(matrix(6, 5, &count_columns, MatrixOrder::ColumnMajor)).unwrap(),
        matrix(6, 5, &cook_columns, MatrixOrder::ColumnMajor),
    )
    .unwrap();
    for heuristic in [native, foreign] {
        let actual = CooksPolicy::new(Some(5.), Some(heuristic))
            .unwrap()
            .evaluate(pvalues, maximum)
            .unwrap();
        assert_eq!(
            actual.outlier,
            vec![
                Some(false),
                Some(true),
                Some(false),
                Some(true),
                Some(false),
                None
            ]
        );
        assert_eq!(
            actual.pvalue,
            vec![Some(0.01), None, Some(0.03), None, Some(0.05), Some(0.06)]
        );
        let refitted = [false, true, false, false, false, false];
        let zeros = [0, 0, 0, 1, 0, 0];
        let rows = ReplacementRows::new(
            6,
            Some(RowFlags::from_bools(&refitted)),
            Some(RowFlags::from_integers(&zeros).unwrap()),
        )
        .unwrap();
        let replaced = CooksPolicy::new(Some(5.), Some(heuristic))
            .unwrap()
            .with_replacement(rows)
            .evaluate(pvalues, maximum)
            .unwrap();
        assert_eq!(
            replaced.outlier,
            vec![
                Some(false),
                Some(false),
                Some(false),
                Some(false),
                Some(false),
                None
            ]
        );
        assert_eq!(
            replaced.pvalue,
            vec![
                Some(0.01),
                Some(0.02),
                Some(0.03),
                Some(1.),
                Some(0.05),
                Some(0.06)
            ]
        );
    }
}

#[test]
fn replacement_columns_validate_before_mutation_and_preserve_unselected_missing_bits() {
    let zero = [false, true, false];
    let rows = ReplacementRows::new(3, None, Some(RowFlags::from_bools(&zero))).unwrap();
    let missing = f64::from_bits(0x7ff00000000007a2);
    let mut effect = [missing, 4., 7.];
    let mut se = [missing, 0.5, 1.];
    let mut stat = [missing, 8., 7.];
    let mut p = [missing, 0.01, 0.2];
    rows.apply_columns(&mut effect, &mut se, &mut stat, &mut p)
        .unwrap();
    assert_eq!([effect[1], se[1], stat[1], p[1]], [0., 0., 0., 1.]);
    assert_eq!([effect[2], se[2], stat[2], p[2]], [7., 1., 7., 0.2]);
    assert_eq!(effect[0].to_bits(), missing.to_bits());
    assert_eq!(p[0].to_bits(), missing.to_bits());
    let before = effect.map(f64::to_bits);
    assert!(
        rows.apply_columns(&mut effect, &mut se, &mut stat, &mut p[..2])
            .is_err()
    );
    assert_eq!(effect.map(f64::to_bits), before);
    assert!(RowFlags::from_integers(&[0, i32::MIN]).is_err());
    assert!(ReplacementRows::new(3, None, Some(RowFlags::from_bools(&[false]))).is_err());
    assert!(CooksPolicy::new(Some(f64::NAN), None).is_err());
    let policy = CooksPolicy::new(Some(5.), None).unwrap();
    assert!(policy.evaluate([Some(0.1)], [Some(f64::INFINITY)]).is_err());
    assert!(policy.evaluate([Some(0.1)], []).is_err());
    assert!(policy.evaluate([Some(1.2)], [Some(0.)]).is_err());
    // As on the legacy path, a masked row has no probability left to validate.
    assert_eq!(
        policy.evaluate([Some(1.2)], [Some(10.)]).unwrap().pvalue,
        vec![None]
    );
    for invalid in [-1., 0.5, f64::NAN, f64::INFINITY, f64::from(u32::MAX) + 1.] {
        assert!(
            CountView::from_doubles(matrix(1, 2, &[1., invalid], MatrixOrder::RowMajor)).is_err()
        );
    }
}

#[test]
fn two_level_eligibility_matches_partitions_without_quadratic_outer_storage() {
    let values = [1., -2., 1., -2., 1., 3., 1., 3.];
    let design = matrix(4, 2, &values, MatrixOrder::RowMajor);
    assert!(two_level_design_matches(design, &[2, 2, 1, 1], 2).unwrap());
    assert!(!two_level_design_matches(design, &[2, 1, 2, 1], 2).unwrap());
    assert!(!two_level_design_matches(design, &[1, 1, 2, 3], 3).unwrap());
    assert!(two_level_design_matches(design, &[1, 1, 2], 2).is_err());
    let additive = [1., 0., 0., 1., 0., 1., 1., 1., 0., 1., 1., 1.];
    assert!(
        !two_level_design_matches(
            matrix(4, 3, &additive, MatrixOrder::RowMajor),
            &[1, 1, 2, 2],
            2
        )
        .unwrap()
    );
    let invalid = [1., f64::NAN, 1., 0., 1., 1., 1., 1.];
    assert!(
        !two_level_design_matches(
            matrix(4, 2, &invalid, MatrixOrder::RowMajor),
            &[1, 1, 2, 2],
            2
        )
        .unwrap()
    );
}

#[test]
fn stored_requests_infer_matching_formula_context_and_preserve_explicit_override() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![1, 5, 6, 7]).unwrap();
    let design =
        DesignMatrix::from_row_major(4, 2, vec![1., 0., 1., 0., 1., 1., 1., 1.], None).unwrap();
    let matrix_origin = design.clone();
    let design = design.with_origin(DesignOrigin::formula(vec!["condition".into()]).unwrap());
    let (mut fit, raw) = DeseqBuilder::new()
        .size_factors(vec![1.; 4])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald(&counts, &design, &[0.1], 1)
        .unwrap();
    assert!(raw.rows[0].pvalue.is_some());
    fit.max_cooks = Some(vec![Some(10.)]);
    fit.cooks = Some(RowMajorMatrix::from_row_major(1, 4, vec![10., 0.1, 0.2, 0.3]).unwrap());
    fit.model_frame = Some(FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".into(),
            sample_levels: ["A", "A", "B", "B"].map(String::from).to_vec(),
            levels: Some(vec!["A".into(), "B".into()]),
            reference: Some("A".into()),
        }],
        ..Default::default()
    });
    let request = ResultRequest::new(
        ResultsContrast::Numeric(vec![0., 1.]),
        WaldTestOptions::normal(),
        CooksCutoff::Threshold(5.),
        IndependentFilteringOptions {
            enabled: false,
            ..Default::default()
        },
    )
    .unwrap();
    let unmasked = fit.results_for_request(&counts, &request, None).unwrap();
    assert_eq!(unmasked.rows[0].pvalue, raw.rows[0].pvalue);
    assert_eq!(unmasked.rows[0].cooks_outlier, Some(false));
    let forced = fit
        .results_for_request(&counts, &request.clone().low_count_heuristic(false), None)
        .unwrap();
    assert_eq!(forced.rows[0].pvalue, None);
    assert_eq!(forced.rows[0].padj, None);
    fit.design = Some(matrix_origin);
    let matrix = fit.results_for_request(&counts, &request, None).unwrap();
    assert_eq!(matrix.rows[0].pvalue, None);
    assert_eq!(matrix.rows[0].cooks_outlier, Some(true));
    fit.model_frame = None;
    assert_eq!(
        fit.results_for_request(&counts, &request, None)
            .unwrap()
            .rows[0]
            .pvalue,
        None
    );
}

#[test]
fn signed_zero_maxima_keep_the_first_numerically_equal_sample() {
    let counts = [9, 1, 5, 6, 7];
    for cooks in [[-0., 0., -1., -2., -3.], [0., -0., -1., -2., -3.]] {
        let heuristic = LowCountHeuristic::new(
            CountView::from_unsigned(matrix(1, 5, &counts, MatrixOrder::RowMajor)),
            matrix(1, 5, &cooks, MatrixOrder::RowMajor),
        )
        .unwrap();
        // Selecting the second sample would spare this row. First-sample
        // selection must mask it, regardless of the sign of either zero.
        let mask = CooksPolicy::new(Some(-0.5), Some(heuristic))
            .unwrap()
            .evaluate([Some(0.1)], [Some(0.)])
            .unwrap();
        assert_eq!(mask.pvalue, vec![None]);
    }
}
