use rsdeseq2::prelude::{CountMatrix, DeseqBuilder, DesignMatrix};

#[test]
fn default_rank_is_invariant_to_covariate_units() {
    // R 4.6.1 qr() returns rank 2 for these three independent designs.
    for scale in [1e-12, 1.0, 1e12] {
        let values: Vec<_> = (0..8)
            .flat_map(|i| [1.0, (i as f64 - 3.5) * scale])
            .collect();
        let design = DesignMatrix::from_row_major(8, 2, values.clone(), None).unwrap();
        assert_eq!(design.rank().unwrap(), 2);
        design.validate_full_rank("mixed units").unwrap();
        assert_eq!(design.matrix().as_slice(), values);
        if scale == 1e-12 {
            assert_eq!(design.rank_with_tolerance(1e-10).unwrap(), 1);
        }
    }
    // Exact and near aliases also agree with the pinned R rank decisions.
    for slope in [0.0, 1e-12] {
        let values = (0..8)
            .flat_map(|i| [1.0, 1.0 + (i as f64 - 3.5) * slope])
            .collect();
        let design = DesignMatrix::from_row_major(8, 2, values, None).unwrap();
        assert_eq!(design.rank().unwrap(), 1);
        assert!(design.validate_full_rank("alias").is_err());
    }
}

#[test]
fn non_nested_lrt_fails_before_dispersion_estimation() {
    let counts = CountMatrix::from_row_major_u32(2, 8, vec![0; 16]).unwrap();
    let full = DesignMatrix::from_row_major(
        8,
        2,
        (0..8)
            .flat_map(|i| [1.0, if i < 4 { 0.0 } else { 1.0 }])
            .collect(),
        None,
    )
    .unwrap();
    let reduced = DesignMatrix::from_row_major(
        8,
        1,
        (0..8).map(|i| if i % 2 == 0 { 1.0 } else { 2.0 }).collect(),
        None,
    )
    .unwrap();
    let err = DeseqBuilder::new()
        .fit_lrt_with_results(&counts, &full, &reduced)
        .unwrap_err();
    assert!(err.to_string().contains("nested"), "{err}");
}
