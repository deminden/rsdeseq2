use rsdeseq2::prelude::*;
use rsdeseq2::transform::{FrozenLocalDispersion, SavedLocalDispersionFit, vst_frozen_local};
use serde::Deserialize;

#[derive(Deserialize)]
struct References {
    cases: Vec<Case>,
    minimum_floor: Floor,
}
#[derive(Deserialize)]
struct Floor {
    means: Vec<String>,
    values: Vec<String>,
}
#[derive(Deserialize)]
struct Case {
    name: String,
    knots: Vec<String>,
    residual_values: Vec<String>,
    residual_slopes: Vec<String>,
    global_polynomial: Vec<String>,
    queries: Vec<String>,
    predictions: Vec<String>,
    genes: usize,
    samples: usize,
    counts: Vec<u32>,
    size_factors: Option<Vec<String>>,
    normalization_factors: Option<Vec<String>>,
    grid: Vec<String>,
    dispersions: Vec<String>,
    transformed: Vec<String>,
    serialization_identical: bool,
}
fn numbers(x: &[String]) -> Vec<f64> {
    x.iter().map(|x| x.parse().unwrap()).collect()
}
fn close(actual: &[f64], expected: &[String], tolerance: f64, context: &str) {
    assert_eq!(actual.len(), expected.len());
    for (i, (a, b)) in actual.iter().zip(numbers(expected)).enumerate() {
        assert!(
            a.is_finite() && (a - b).abs() <= tolerance * b.abs().max(1.0),
            "{context}[{i}]: {a:.17e} != {b:.17e}"
        );
    }
}
#[test]
fn saved_original_locfit_predictions_and_full_vst_preserve_independent_bounds() {
    let references: References =
        serde_json::from_str(include_str!("data/frozen_local_vst_r462.json")).unwrap();
    for case in references.cases {
        assert!(case.serialization_identical);
        let fit = SavedLocalDispersionFit::new(
            numbers(&case.knots),
            numbers(&case.residual_values),
            numbers(&case.residual_slopes),
            numbers(&case.global_polynomial).try_into().unwrap(),
        )
        .unwrap();
        close(
            &fit.evaluate(&numbers(&case.queries)).unwrap(),
            &case.predictions,
            2e-12,
            &case.name,
        );
        let counts =
            CountMatrix::from_row_major_u32(case.genes, case.samples, case.counts).unwrap();
        let (normalized, inverse) = if let Some(nf) = &case.normalization_factors {
            let nf = RowMajorMatrix::from_row_major(case.genes, case.samples, numbers(nf)).unwrap();
            (
                normalized_counts_with_factors(&counts, &nf).unwrap(),
                local_vst_inverse_size_factor_mean_from_normalization_factors(&nf).unwrap(),
            )
        } else {
            let sf = numbers(case.size_factors.as_ref().unwrap());
            (
                normalized_counts(&counts, &sf).unwrap(),
                local_vst_inverse_size_factor_mean(&sf).unwrap(),
            )
        };
        let curve = FrozenLocalDispersion::Saved(fit);
        let (grid, dispersions) = curve.grid(&normalized).unwrap();
        let borrowed_grid = if let Some(nf) = &case.normalization_factors {
            let nf = numbers(nf);
            frozen_local_vst_grid_from_counts(
                (&counts).into(),
                None,
                Some(
                    MatrixView::new(case.genes, case.samples, &nf, MatrixOrder::RowMajor).unwrap(),
                ),
            )
            .unwrap()
        } else {
            frozen_local_vst_grid_from_counts(
                (&counts).into(),
                Some(&numbers(case.size_factors.as_ref().unwrap())),
                None,
            )
            .unwrap()
        };
        assert_eq!(
            grid.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            borrowed_grid
                .iter()
                .map(|x| x.to_bits())
                .collect::<Vec<_>>()
        );
        close(&grid, &case.grid, 2e-12, &format!("{} grid", case.name));
        close(
            &dispersions,
            &case.dispersions,
            2e-12,
            &format!("{} dispersions", case.name),
        );
        let transformed = vst_frozen_local(&normalized, &curve, inverse).unwrap();
        close(
            transformed.as_slice(),
            &case.transformed,
            2e-6,
            &format!("{} VST", case.name),
        );
    }
}

#[test]
fn malformed_saved_fits_and_invalid_predictions_fail_before_indexing() {
    for (knots, values, slopes, polynomial) in [
        (vec![], vec![], vec![], [0.; 4]),
        (vec![1., 1.], vec![0., 0.], vec![0., 0.], [0.; 4]),
        (vec![1., 2.], vec![0.], vec![0., 0.], [0.; 4]),
        (vec![1., 2.], vec![0., 0.], vec![f64::NAN, 0.], [0.; 4]),
        (vec![1., 2.], vec![0., 0.], vec![0., 0.], [f64::INFINITY; 4]),
    ] {
        assert!(SavedLocalDispersionFit::new(knots, values, slopes, polynomial).is_err());
    }
    let fit =
        SavedLocalDispersionFit::new(vec![2., 1.], vec![0., 0.], vec![0., 0.], [0.; 4]).unwrap();
    assert_eq!(fit.evaluate(&[1., 10.]).unwrap(), vec![1., 1.]);
    for value in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(fit.evaluate(&[value]).is_err());
    }
}

#[test]
fn minimum_floor_expression_and_zero_count_missingness_are_preserved() {
    use rsdeseq2::design::expression::{Argument, Expression, Value};
    let call = |function: &str, expressions: Vec<Expression>| Expression::Call {
        function: function.into(),
        arguments: expressions
            .into_iter()
            .map(|expression| Argument {
                name: None,
                expression,
            })
            .collect(),
    };
    let curve = FrozenLocalDispersion::Expression {
        argument: "mu".into(),
        constants: vec![],
        expression: call(
            "rep",
            vec![
                Expression::Literal(Value::scalar(1e-8)),
                call("length", vec![Expression::Symbol("mu".into())]),
            ],
        ),
    };
    let reference: References =
        serde_json::from_str(include_str!("data/frozen_local_vst_r462.json")).unwrap();
    let predicted = curve
        .evaluate(&numbers(&reference.minimum_floor.means))
        .unwrap();
    assert_eq!(predicted, numbers(&reference.minimum_floor.values));
    let zero = RowMajorMatrix::from_elem(3, 2, 0.).unwrap();
    assert!(
        vst_frozen_local(&zero, &curve, 1.)
            .unwrap()
            .as_slice()
            .iter()
            .all(|v| v.is_nan())
    );
    for value in [
        Value::scalar(0.1),
        Value::numbers(vec![-1.; 3]),
        Value::numbers(vec![f64::INFINITY; 3]),
    ] {
        let invalid = FrozenLocalDispersion::Expression {
            argument: "mu".into(),
            expression: Expression::Literal(value),
            constants: vec![],
        };
        assert!(invalid.evaluate(&[1., 2., 3.]).is_err());
    }
}
