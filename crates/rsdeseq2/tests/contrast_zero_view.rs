use rsdeseq2::prelude::*;

fn matrix<T>(rows: usize, columns: usize, values: &[T], order: MatrixOrder) -> MatrixView<'_, T> {
    MatrixView::new(rows, columns, values, order).unwrap()
}

#[test]
fn borrowed_count_layouts_share_zero_flags_and_reusable_global_context() {
    let counts = [
        0, 0, 0, 0, 5, 6, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 5, 6, 0, 0, 1, 0, 0, 0,
    ];
    let column = (0..6)
        .flat_map(|sample| (0..4).map(move |gene| f64::from(counts[gene * 6 + sample])))
        .collect::<Vec<_>>();
    let row = ContrastCounts::from_unsigned(matrix(4, 6, &counts, MatrixOrder::RowMajor));
    let column =
        ContrastCounts::from_doubles(matrix(4, 6, &column, MatrixOrder::ColumnMajor)).unwrap();
    let design = [
        1., 0., 0., 1., 0., 0., 0., 1., 0., 0., 1., 0., 0., 0., 1., 0., 0., 1.,
    ];
    for view in [&row, &column] {
        assert_eq!(view.all_zero(), &[false, true, false, false]);
        let numeric = view
            .numeric_flags(
                Some(matrix(6, 3, &design, MatrixOrder::RowMajor)),
                &[-1., 1., 0.],
            )
            .unwrap();
        assert_eq!(numeric, vec![true, true, false, false]);
        let factor = view
            .factor_flags(
                &["A cell", "A cell", "B cell", "B cell", "C cell", "C cell"],
                "B.cell",
                "A cell",
            )
            .unwrap();
        assert_eq!(numeric, factor);
        assert_eq!(
            view.numeric_flags(None, &[0., 1., 0.]).unwrap(),
            vec![false; 4]
        );
        let next = view
            .numeric_flags(
                Some(matrix(6, 3, &design, MatrixOrder::RowMajor)),
                &[0., -1., 1.],
            )
            .unwrap();
        assert_eq!(next, vec![false, true, false, false]);
        assert_eq!(view.all_zero(), &[false, true, false, false]);
        let policy = ContrastZeroPolicy::new(&factor, view.all_zero()).unwrap();
        let mut effect = [Some(2.); 4];
        let mut stat = [Some(4.); 4];
        let mut p = [Some(0.01); 4];
        policy.apply_wald(&mut effect, &mut stat, &mut p).unwrap();
        assert_eq!(effect, [Some(0.), Some(2.), Some(2.), Some(2.)]);
        assert_eq!(stat, [Some(0.), Some(4.), Some(4.), Some(4.)]);
        assert_eq!(p, [Some(1.), Some(0.01), Some(0.01), Some(0.01)]);
        let mut display = [2.; 4];
        policy.apply_effects(display.iter_mut()).unwrap();
        assert_eq!(display, [0., 2., 2., 2.]);
    }
}

#[test]
fn primitive_zero_context_requires_only_inputs_needed_by_the_request() {
    let counts = [0, 0, 0, 2, 0, 1, 0, 0];
    let view = ContrastCounts::from_unsigned(matrix(2, 4, &counts, MatrixOrder::RowMajor));
    assert!(
        view.numeric_flags(None, &[-1., 1.])
            .unwrap_err()
            .to_string()
            .contains("modelMatrix")
    );
    assert!(view.numeric_flags(None, &[0., 0.]).is_err());
    assert!(view.numeric_flags(None, &[1., f64::NAN]).is_err());
    assert!(
        view.numeric_flags(
            Some(matrix(2, 4, &[1.; 8], MatrixOrder::RowMajor)),
            &[1., -1., 0., 0.]
        )
        .is_err()
    );
    assert!(
        view.numeric_flags(
            Some(matrix(4, 1, &[1.; 4], MatrixOrder::RowMajor)),
            &[-1., 1.]
        )
        .is_err()
    );
    for value in [f64::NAN, f64::INFINITY, -1., 0.5, f64::from(u32::MAX) + 1.] {
        assert!(
            ContrastCounts::from_doubles(matrix(1, 2, &[0., value], MatrixOrder::ColumnMajor))
                .is_err()
        );
    }
    let extreme_counts = [-0., f64::from(u32::MAX)];
    let extreme =
        ContrastCounts::from_doubles(matrix(1, 2, &extreme_counts, MatrixOrder::RowMajor)).unwrap();
    assert_eq!(extreme.all_zero(), &[false]);
    let empty = ContrastCounts::from_doubles(matrix(0, 2, &[], MatrixOrder::ColumnMajor)).unwrap();
    assert!(empty.all_zero().is_empty());
    assert!(empty.numeric_flags(None, &[1., 0.]).unwrap().is_empty());
}

#[test]
fn zero_policy_validates_every_column_before_mutation_and_excludes_global_zeros() {
    let flags = [true, true, true];
    let global = [false, true, false];
    let policy = ContrastZeroPolicy::new(&flags, &global).unwrap();
    let mut effect = [Some(2.); 3];
    let mut stat = [Some(4.); 3];
    let mut p = [Some(0.01); 2];
    assert!(policy.apply_wald(&mut effect, &mut stat, &mut p).is_err());
    assert_eq!(effect, [Some(2.); 3]);
    assert_eq!(stat, [Some(4.); 3]);
    assert_eq!(p, [Some(0.01); 2]);
    assert!(ContrastZeroPolicy::new(&flags, &global[..2]).is_err());
    let mut effect = [None; 3];
    let mut stat = [None; 3];
    let mut p = [None; 3];
    policy.apply_wald(&mut effect, &mut stat, &mut p).unwrap();
    assert_eq!(effect, [Some(0.), None, Some(0.)]);
    assert_eq!(stat, [Some(0.), None, Some(0.)]);
    assert_eq!(p, [Some(1.), None, Some(1.)]);
}
