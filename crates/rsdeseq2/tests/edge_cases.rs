//! Degenerate inputs must return a typed error or a well-formed result, never
//! panic. The validation cases check the error fields exposed to callers.

use crate::support::run::run_results;
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

fn kind<T: std::fmt::Debug>(result: Result<T, DeseqError>) -> ErrorKind {
    result.unwrap_err().kind()
}

fn fit(
    counts: &CountMatrix,
    design: &DesignMatrix,
) -> Result<(DeseqFit, DeseqResults), DeseqError> {
    run_results(DeseqBuilder::new().run(&Analysis::new(counts, design)))
}

#[test]
fn empty_and_mismatched_shapes_are_rejected() {
    assert_eq!(
        kind(CountMatrix::from_row_major_u32(0, 6, vec![])),
        ErrorKind::Counts
    );
    assert_eq!(
        kind(CountMatrix::from_row_major_u32(5, 0, vec![])),
        ErrorKind::Counts
    );
    assert_eq!(
        kind(CountMatrix::from_row_major_u32(2, 3, vec![1; 5])),
        ErrorKind::Shape
    );
    assert_eq!(
        kind(DesignMatrix::from_row_major(0, 2, vec![], None)),
        ErrorKind::Shape
    );
    assert_eq!(
        kind(fit(&varied(50, 6), &design(&[0, 0, 1, 1]))),
        ErrorKind::Shape
    );
}

#[test]
fn designs_without_replicates_or_full_rank_are_rejected() {
    let counts = varied(50, 6);
    let constant = DesignMatrix::from_row_major(6, 2, vec![1.0; 12], None).unwrap();
    std::assert_matches!(
        fit(&counts, &constant).unwrap_err(),
        DeseqError::DesignRankDeficient { rank: 1, columns: 2, zero_columns, .. }
            if zero_columns.is_empty()
    );
    assert_eq!(kind(fit(&counts, &design(&[1; 6]))), ErrorKind::Options);
    // DESeq2: "the same number of samples and coefficients to fit".
    assert_eq!(
        kind(fit(&varied(50, 2), &design(&[0, 1]))),
        ErrorKind::Shape
    );
    assert_eq!(
        kind(DesignMatrix::from_row_major(
            2,
            1,
            vec![1.0, f64::NAN],
            None
        )),
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
        std::assert_matches!(
            fit(&counts, &design).unwrap_err(),
            DeseqError::NoUsableGenesForSizeFactors
        );
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
            std::assert_matches!(
                DeseqBuilder::new()
                    .fit_type(fit_type)
                    .run(&Analysis::new(counts, &design)).unwrap_err(),
                DeseqError::NoUsableDispersionTrendRows { min_disp: 1e-8, genes_checked }
                    if genes_checked == counts.n_genes(),
                "{fit_type:?}"
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
        let result = run_results(
            DeseqBuilder::new()
                .size_factors(factors)
                .run(&Analysis::new(&counts, &design)),
        );
        assert_eq!(kind(result), ErrorKind::Normalization);
    }
    let short = run_results(
        DeseqBuilder::new()
            .size_factors(vec![1.0; 5])
            .run(&Analysis::new(&counts, &design)),
    );
    assert_eq!(kind(short), ErrorKind::Shape);

    let small = varied(5, 6);
    let wald = |dispersions: &[f64], coefficient| {
        run_results(
            DeseqBuilder::new().run(
                &Analysis::new(&small, &design)
                    .fixed_dispersion(dispersions)
                    .result(ResultColumn::Coefficient(coefficient)),
            ),
        )
    };
    for rejected in [-0.1, f64::NAN] {
        std::assert_matches!(
            wald(&[0.1, rejected, 0.1, 0.1, 0.1], 1).unwrap_err(),
            DeseqError::NonFiniteValue { context, index: Some(1), value }
                if context == "dispersion" && value.to_bits() == rejected.to_bits()
        );
    }
    assert_eq!(kind(wald(&[0.1; 4], 1)), ErrorKind::Shape);
    assert_eq!(kind(wald(&[0.1; 5], 2)), ErrorKind::Shape);
}

#[test]
fn invalid_contrasts_and_threads_are_rejected() {
    let small = varied(5, 6);
    let design = design(&GROUPS);
    let contrast = |weights: &[f64]| {
        run_results(
            DeseqBuilder::new().run(
                &Analysis::new(&small, &design)
                    .fixed_dispersion(&[0.1; 5])
                    .result(ResultColumn::Numeric(weights)),
            ),
        )
    };
    assert_eq!(kind(contrast(&[0.0, 1.0, 0.0])), ErrorKind::Shape);
    assert_eq!(kind(contrast(&[0.0, 0.0])), ErrorKind::Options);
    let threads = run_results(
        DeseqBuilder::new()
            .threads(0)
            .run(&Analysis::new(&varied(50, 6), &design)),
    );
    assert_eq!(kind(threads), ErrorKind::Options);
}

#[test]
fn oversized_matrices_fail_with_an_error_instead_of_aborting() {
    // 2^58 values (2 EiB): the size passes the `isize::MAX` byte check but is
    // beyond every 64-bit virtual address space, so the reservation fails
    // on every host. (16 TiB, used before, is reserved lazily on macOS, and
    // `from_elem` then filled it until the machine ran out of memory.)
    let error = RowMajorMatrix::from_elem(1usize << 48, 1usize << 10, 0.0f64).unwrap_err();
    let bytes = (1usize << 58) * std::mem::size_of::<f64>();
    let expected_context = format!("row-major matrix allocation failed ({bytes} bytes requested)");
    assert!(
        matches!(&error, DeseqError::InvalidDimensions { context, .. } if context == &expected_context),
        "{error:?}"
    );
    let error = RowMajorMatrix::from_elem(usize::MAX, 2, 0.0f64).unwrap_err();
    assert!(matches!(error, DeseqError::InvalidDimensions { .. }));
}
