//! Independently frozen R arithmetic plus shared raw/fitted exclusion policy.
use rsdeseq2::matrix::{CountRows, CountView, MatrixOrder, MatrixView};
use rsdeseq2::prelude::*;
use serde::Deserialize;

#[derive(Deserialize)]
struct Reference {
    offsets: Vec<Offset>,
    rows: Vec<Rows>,
}
#[derive(Deserialize)]
struct Offset {
    name: String,
    genes: usize,
    samples: usize,
    size_factors: Option<Vec<String>>,
    factors: Option<Vec<String>>,
    offsets: Vec<String>,
    expected: Option<Vec<String>>,
}
#[derive(Deserialize)]
struct Rows {
    name: String,
    genes: usize,
    samples: usize,
    counts: Vec<u32>,
    dispersions: Vec<String>,
    fitted: Vec<String>,
    excluded: Option<Vec<bool>>,
    expected: Option<Vec<usize>>,
}
fn reference() -> Reference {
    serde_json::from_str(include_str!("data/row_preparation_r462.json")).unwrap()
}
fn numbers(v: &[String]) -> Vec<f64> {
    v.iter()
        .map(|v| f64::from_bits(u64::from_str_radix(v, 16).unwrap()))
        .collect()
}
fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|v| v.to_bits()).collect()
}
fn columns<T: Copy>(v: &[T], rows: usize, cols: usize) -> Vec<T> {
    (0..cols)
        .flat_map(|c| (0..rows).map(move |r| v[r * cols + c]))
        .collect()
}
#[test]
fn base_factor_offset_products_match_independent_r_bits_in_both_layouts() {
    let r = reference();
    assert_eq!(r.offsets.len(), 12);
    for case in r.offsets {
        let sf = case.size_factors.as_deref().map(numbers);
        let nf = case.factors.as_deref().map(numbers);
        let offsets = numbers(&case.offsets);
        let column = nf.as_ref().map(|v| columns(v, case.genes, case.samples));
        for order in [MatrixOrder::RowMajor, MatrixOrder::ColumnMajor] {
            let nf = match order {
                MatrixOrder::RowMajor => nf.as_ref(),
                MatrixOrder::ColumnMajor => column.as_ref(),
            }
            .map(|v| MatrixView::new(case.genes, case.samples, v, order).unwrap());
            let actual = normalization_factors_with_offsets(
                case.genes,
                case.samples,
                sf.as_deref(),
                nf,
                &offsets,
            );
            if let Some(expected) = &case.expected {
                assert_eq!(
                    bits(actual.unwrap().as_slice()),
                    bits(&numbers(expected)),
                    "{} {order:?}",
                    case.name
                );
            } else {
                assert!(actual.is_err(), "{} {order:?}", case.name);
            }
        }
    }
}
#[test]
fn shrink_rows_preserve_raw_zero_and_fitted_exclusion_semantics_across_layouts() {
    let r = reference();
    assert_eq!(r.rows.len(), 6);
    for case in r.rows {
        let unsigned = CountRows::from_unsigned(
            MatrixView::new(
                case.genes,
                case.samples,
                &case.counts,
                MatrixOrder::RowMajor,
            )
            .unwrap(),
        );
        let counts = columns(&case.counts, case.genes, case.samples)
            .into_iter()
            .map(f64::from)
            .collect::<Vec<_>>();
        let foreign = CountRows::from_doubles(
            MatrixView::new(case.genes, case.samples, &counts, MatrixOrder::ColumnMajor).unwrap(),
        )
        .unwrap();
        assert_eq!(unsigned.all_zero(), foreign.all_zero(), "{}", case.name);
        let d = numbers(&case.dispersions);
        let f = numbers(&case.fitted);
        for rows in [&unsigned, &foreign] {
            let selected = NormalShrinkageRows::new(
                rows,
                &d,
                &f,
                case.excluded.as_deref().map(RowFlags::from_bools),
            );
            if let Some(expected) = &case.expected {
                assert_eq!(selected.unwrap().indices(), expected, "{}", case.name);
            } else {
                assert!(selected.is_err(), "{}", case.name);
            }
        }
    }
}
#[test]
fn preparation_rejects_bad_shapes_flags_counts_and_unusable_bases() {
    let counts = [0, 0, 1, 2];
    let rows =
        CountRows::from_unsigned(MatrixView::new(2, 2, &counts, MatrixOrder::RowMajor).unwrap());
    assert!(NormalShrinkageRows::new(&rows, &[0.1], &[0.1, 0.1], None).is_err());
    assert!(
        NormalShrinkageRows::new(
            &rows,
            &[0.1, 0.1],
            &[0.1, 0.1],
            Some(RowFlags::from_bools(&[false]))
        )
        .is_err()
    );
    for bad in [f64::NAN, f64::INFINITY, -1., 0.5, 4294967296.] {
        let values = [0., bad];
        let matrix = MatrixView::new(1, 2, &values, MatrixOrder::RowMajor).unwrap();
        assert!(CountView::from_doubles(matrix).is_err());
        assert!(CountRows::from_doubles(matrix).is_err());
    }
    assert!(normalization_factors_with_offsets(2, 2, None, None, &[0., 0.]).is_err());
    assert!(normalization_factors_with_offsets(2, 2, Some(&[1., 1.]), None, &[0.]).is_err());
    assert!(normalization_factors_with_offsets(2, 2, Some(&[1., 0.]), None, &[0., 0.]).is_err());
    let wrong = [1.; 6];
    assert!(
        normalization_factors_with_offsets(
            2,
            2,
            None,
            Some(MatrixView::new(3, 2, &wrong, MatrixOrder::RowMajor).unwrap()),
            &[0., 0.]
        )
        .is_err()
    );
}

#[test]
fn formula_offsets_use_effective_factors_without_dividing_by_tiny_base_factors() {
    let counts = CountMatrix::from_row_major_u32(
        3,
        6,
        vec![
            1, 20, 21, 20, 20, 20, 10, 11, 12, 15, 16, 17, 30, 32, 31, 50, 52, 51,
        ],
    )
    .unwrap();
    let frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".into(),
            sample_levels: ["A", "A", "A", "B", "B", "B"].map(String::from).to_vec(),
            levels: Some(vec!["A".into(), "B".into()]),
            reference: Some("A".into()),
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "exposure".into(),
            values: vec![700.; 6],
        }],
    };
    let builder = DeseqBuilder::new()
        .threads(1)
        .fit_type(FitType::Mean)
        .model_frame(frame)
        .size_factors(vec![1e-310; 6])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..Default::default()
        });
    // Shape errors precede size-factor estimation, even when the counts would
    // make ratio estimation fail. Keep the earlier native preflight contract.
    let wrong_counts = CountMatrix::from_row_major_u32(3, 2, vec![0; 6]).unwrap();
    let error = DeseqBuilder::new()
        .model_frame(builder.current_model_frame().unwrap().clone())
        .fit_formula_with_results(&wrong_counts, "~condition+offset(exposure)")
        .unwrap_err()
        .to_string();
    assert!(error.contains("formula offsets"), "{error}");
    let formula = builder
        .expanded_formula_design_with_offsets("~condition+offset(exposure)")
        .unwrap();
    let (fit, actual) = builder
        .clone()
        .fit_formula_with_results(&counts, "~condition+offset(exposure)")
        .unwrap();
    let factors =
        normalization_factors_with_offsets(3, 6, Some(&[1e-310; 6]), None, &[700.; 6]).unwrap();
    let (matrix_fit, expected) = builder
        .normalization_factors(factors.clone())
        .fit_with_results(&counts, &formula.design.standard_design)
        .unwrap();
    assert_eq!(actual.rows, expected.rows);
    assert_eq!(fit.normalization_factors.as_ref(), Some(&factors));
    // Complete formula workflows consume sample scales into the effective
    // matrix, matching wrapped fits. Explicit matrix inputs retain their
    // supplied sample-scale metadata even though the matrix takes precedence.
    assert_eq!(fit.size_factors, vec![1.; 6]);
    assert_eq!(matrix_fit.size_factors, vec![1e-310; 6]);
}
