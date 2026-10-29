use rsdeseq2::prelude::{CountMatrix, DeseqBuilder, DesignMatrix};

#[test]
fn default_rank_is_invariant_to_covariate_units() {
    // Independent R 4.6.2 decisions are frozen in rank_r462.tsv below.
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
    // Exact and sufficiently near aliases agree with both pinned policies.
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
fn default_rank_matches_frozen_scaled_r462_qr_decisions() {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(include_bytes!("data/rank_r462.tsv").as_slice());
    for row in reader.records() {
        let row = row.unwrap();
        let values: Vec<f64> = row[5].split(',').map(|v| v.parse().unwrap()).collect();
        let design = DesignMatrix::from_row_major(
            row[1].parse().unwrap(),
            row[2].parse().unwrap(),
            values.clone(),
            None,
        )
        .unwrap();
        assert_eq!(
            design.rank().unwrap(),
            row[4].parse::<usize>().unwrap(),
            "{}",
            &row[0]
        );
        // Validation never changes coefficient/fitting units.
        assert_eq!(design.matrix().as_slice(), values);
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
