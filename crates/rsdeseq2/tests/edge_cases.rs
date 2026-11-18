//! Degenerate inputs must return a typed error or a well-formed result, never
//! panic. Where DESeq2 1.52 stops with a specific message, the message is
//! checked too.

use rsdeseq2::prelude::*;

const GROUPS: [usize; 6] = [0, 0, 0, 1, 1, 1];

fn design(groups: &[usize]) -> DesignMatrix {
    let values = groups
        .iter()
        .flat_map(|group| [1.0, if *group == 1 { 1.0 } else { 0.0 }])
        .collect();
    DesignMatrix::from_row_major(groups.len(), 2, values, None).unwrap()
}

fn counts(genes: usize, samples: usize, value: impl Fn(usize, usize) -> u32) -> CountMatrix {
    let values = (0..genes * samples)
        .map(|index| value(index / samples, index % samples))
        .collect();
    CountMatrix::from_row_major_u32(genes, samples, values).unwrap()
}

fn varied(genes: usize, samples: usize) -> CountMatrix {
    counts(genes, samples, |gene, sample| {
        (10 + gene * 37 + sample * 13 + (gene * sample * 7919) % 97) as u32
    })
}

fn kind<T: std::fmt::Debug>(result: Result<T, DeseqError>) -> (ErrorKind, String) {
    let error = result.unwrap_err();
    (error.kind(), error.to_string())
}

fn fit(
    counts: &CountMatrix,
    design: &DesignMatrix,
) -> Result<(DeseqFit, DeseqResults), DeseqError> {
    DeseqBuilder::new().fit_with_results(counts, design)
}

const WITHIN_TWO_ORDERS: &str = "within 2 orders of magnitude from the minimum value";

#[test]
fn empty_and_mismatched_shapes_are_rejected() {
    assert_eq!(
        kind(CountMatrix::from_row_major_u32(0, 6, vec![])).0,
        ErrorKind::Counts
    );
    assert_eq!(
        kind(CountMatrix::from_row_major_u32(5, 0, vec![])).0,
        ErrorKind::Counts
    );
    assert_eq!(
        kind(CountMatrix::from_row_major_u32(2, 3, vec![1; 5])).0,
        ErrorKind::Shape
    );
    assert_eq!(
        kind(DesignMatrix::from_row_major(0, 2, vec![], None)).0,
        ErrorKind::Shape
    );
    assert_eq!(
        kind(fit(&varied(50, 6), &design(&[0, 0, 1, 1]))).0,
        ErrorKind::Shape
    );
}

#[test]
fn designs_without_replicates_or_full_rank_are_rejected() {
    let counts = varied(50, 6);
    let constant = DesignMatrix::from_row_major(6, 2, vec![1.0; 12], None).unwrap();
    let (error_kind, message) = kind(fit(&counts, &constant));
    assert_eq!(error_kind, ErrorKind::Options);
    assert!(message.contains("not full rank"), "{message}");
    assert_eq!(kind(fit(&counts, &design(&[1; 6]))).0, ErrorKind::Options);
    // DESeq2: "the same number of samples and coefficients to fit".
    assert_eq!(
        kind(fit(&varied(50, 2), &design(&[0, 1]))).0,
        ErrorKind::Shape
    );
    assert_eq!(
        kind(DesignMatrix::from_row_major(
            2,
            1,
            vec![1.0, f64::NAN],
            None
        ))
        .0,
        ErrorKind::Numeric
    );
}

#[test]
fn size_factors_need_a_gene_without_zeros() {
    let design = design(&GROUPS);
    // DESeq2: "every gene contains at least one zero, cannot compute log
    // geometric means".
    let one_zero_each = counts(50, 6, |gene, sample| {
        if sample == gene % 6 {
            0
        } else {
            10 + gene as u32
        }
    });
    for counts in [one_zero_each, counts(50, 6, |_, _| 0)] {
        let (error_kind, message) = kind(fit(&counts, &design));
        assert_eq!(error_kind, ErrorKind::Normalization);
        assert!(message.contains("no usable genes"), "{message}");
    }
}

#[test]
fn flat_dispersions_stop_before_any_trend() {
    let design = design(&GROUPS);
    let flat = counts(50, 6, |gene, _| 5 + gene as u32);
    let near_max = counts(50, 6, |gene, sample| {
        u32::MAX - (gene * 1000 + sample * 77) as u32
    });
    for counts in [&flat, &near_max, &varied(1, 6)] {
        for fit_type in [FitType::Parametric, FitType::Local, FitType::Mean] {
            let (error_kind, message) = kind(
                DeseqBuilder::new()
                    .fit_type(fit_type)
                    .fit_with_results(counts, &design),
            );
            assert_eq!(error_kind, ErrorKind::Dispersion, "{fit_type:?}");
            assert!(
                message.contains(WITHIN_TWO_ORDERS),
                "{fit_type:?}: {message}"
            );
        }
    }
}

#[test]
fn extreme_but_valid_inputs_give_well_formed_results() {
    let design = design(&GROUPS);
    let huge_gene = counts(200, 6, |gene, sample| {
        if gene == 0 {
            u32::MAX - sample as u32
        } else {
            (10 + gene * 37 + sample * 13 + (gene * sample * 7919) % 97) as u32
        }
    });
    for counts in [varied(2, 6), varied(200, 6), huge_gene] {
        let (fit, results) = fit(&counts, &design).unwrap();
        assert_eq!(results.rows.len(), counts.n_genes());
        for row in &results.rows {
            assert!(row.base_mean.is_finite() && row.base_mean > 0.0);
            if let Some(pvalue) = row.pvalue {
                assert!((0.0..=1.0).contains(&pvalue), "{pvalue}");
            }
        }
        assert!(
            fit.dispersion
                .unwrap()
                .iter()
                .all(|value| value.is_finite() && *value > 0.0)
        );
    }
}

#[test]
fn invalid_user_supplied_factors_and_dispersions_are_rejected() {
    let counts = varied(50, 6);
    let design = design(&GROUPS);
    for factors in [
        vec![1.0, f64::NAN, 1.0, 1.0, 1.0, 1.0],
        vec![1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
    ] {
        let result = DeseqBuilder::new()
            .size_factors(factors)
            .fit_with_results(&counts, &design);
        assert_eq!(kind(result).0, ErrorKind::Normalization);
    }
    let short = DeseqBuilder::new()
        .size_factors(vec![1.0; 5])
        .fit_with_results(&counts, &design);
    assert_eq!(kind(short).0, ErrorKind::Shape);

    let small = varied(5, 6);
    let wald = |dispersions: &[f64], coefficient| {
        DeseqBuilder::new().fit_fixed_dispersion_wald(&small, &design, dispersions, coefficient)
    };
    let (negative_kind, negative) = kind(wald(&[0.1, -0.1, 0.1, 0.1, 0.1], 1));
    assert_eq!(negative_kind, ErrorKind::Numeric);
    assert!(negative.starts_with("out-of-range value"), "{negative}");
    let (_, nan) = kind(wald(&[0.1, f64::NAN, 0.1, 0.1, 0.1], 1));
    assert!(nan.starts_with("non-finite value"), "{nan}");
    assert_eq!(kind(wald(&[0.1; 4], 1)).0, ErrorKind::Shape);
    assert_eq!(kind(wald(&[0.1; 5], 2)).0, ErrorKind::Shape);
}

#[test]
fn invalid_contrasts_and_threads_are_rejected() {
    let small = varied(5, 6);
    let design = design(&GROUPS);
    let contrast = |weights: &[f64]| {
        DeseqBuilder::new().fit_fixed_dispersion_wald_contrast(&small, &design, &[0.1; 5], weights)
    };
    assert_eq!(kind(contrast(&[0.0, 1.0, 0.0])).0, ErrorKind::Shape);
    assert_eq!(kind(contrast(&[0.0, 0.0])).0, ErrorKind::Options);
    let threads = DeseqBuilder::new()
        .threads(0)
        .fit_with_results(&varied(50, 6), &design);
    assert_eq!(kind(threads).0, ErrorKind::Options);
}

#[test]
fn oversized_matrices_fail_with_an_error_instead_of_aborting() {
    // 2^41 values (16 TiB): the size fits in usize, the allocation cannot succeed.
    let error = RowMajorMatrix::from_elem(1usize << 31, 1usize << 10, 0.0f64).unwrap_err();
    assert!(
        matches!(&error, DeseqError::InvalidDimensions { context, .. } if context == "row-major matrix allocation failed"),
        "{error:?}"
    );
    let error = RowMajorMatrix::from_elem(usize::MAX, 2, 0.0f64).unwrap_err();
    assert!(matches!(error, DeseqError::InvalidDimensions { .. }));
}
