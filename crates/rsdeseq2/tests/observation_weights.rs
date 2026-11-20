use crate::support::two_group_design;
use approx::assert_relative_eq;
use rsdeseq2::errors::{CalculationLocation, ErrorKind, Requirement};
use rsdeseq2::prelude::*;

#[test]
fn observation_weights_are_normalized_by_gene_row_max() {
    let weights =
        RowMajorMatrix::from_row_major(2, 4, vec![2.0, 1.0, 0.0, 2.0, 0.5, 0.5, 0.0, 0.25])
            .unwrap();
    let output =
        preprocess_observation_weights(&weights, &two_group_design(4, Some("condition_B_vs_A")))
            .unwrap();

    assert_eq!(output.design_rank, 2);
    assert_eq!(output.weights_fail, vec![false, false]);
    assert_eq!(output.weights_ok(), vec![true, true]);
    assert_relative_eq!(output.weights.row(0).unwrap()[0], 1.0, epsilon = 1e-12);
    assert_relative_eq!(output.weights.row(0).unwrap()[1], 0.5, epsilon = 1e-12);
    assert_relative_eq!(output.weights.row(1).unwrap()[0], 1.0, epsilon = 1e-12);
    assert_relative_eq!(output.weights.row(1).unwrap()[3], 0.5, epsilon = 1e-12);
}

#[test]
fn observation_weights_flag_rows_that_drop_design_rank() {
    let weights = RowMajorMatrix::from_row_major(
        3,
        4,
        vec![
            1.0, 1.0, 1.0, 1.0, //
            1.0, 1.0, 0.0, 0.0, //
            0.0, 0.0, 0.0, 0.0,
        ],
    )
    .unwrap();
    let output =
        preprocess_observation_weights(&weights, &two_group_design(4, Some("condition_B_vs_A")))
            .unwrap();

    assert_eq!(output.weights_fail, vec![false, true, true]);
    assert_eq!(output.weights.row(2).unwrap(), &[0.0, 0.0, 0.0, 0.0]);
}

#[test]
fn observation_weights_use_thresholded_subset_for_cox_reid_check() {
    let weights = RowMajorMatrix::from_row_major(
        4,
        4,
        vec![
            1.0, 1.0, 0.001, 0.001, // The thresholded design drops the group column.
            0.001, 0.001, 1.0, 1.0, // Both retained columns are constant.
            1.0, 1.0, 1.0, 1.0, // The full thresholded design is restored.
            1.0, 1.0, 0.0, 0.0, // The weighted design itself has lost rank.
        ],
    )
    .unwrap();
    let output = preprocess_observation_weights_with_options(
        &weights,
        &two_group_design(4, Some("condition_B_vs_A")),
        ObservationWeightOptions {
            weight_threshold: 1e-2,
            ..ObservationWeightOptions::default()
        },
    )
    .unwrap();

    assert_eq!(output.weights_fail, vec![false, true, false, true]);
}

#[test]
fn observation_weights_reject_overflowed_cox_reid_column_sum() {
    let design = DesignMatrix::from_row_major(2, 1, vec![f64::MAX, f64::MAX], None).unwrap();
    let weights = RowMajorMatrix::from_row_major(2, 2, vec![1.0, 0.0, 1.0, 1.0]).unwrap();

    let err = preprocess_observation_weights(&weights, &design).unwrap_err();

    assert_eq!(err.kind(), ErrorKind::Numeric);
    assert_eq!(
        err.to_string(),
        format!(
            "out-of-range value in observation weight Cox-Reid column sum at index Some(0): {}",
            f64::MAX
        )
    );
    std::assert_matches!(
        err,
        DeseqError::CoxReidColumnSumOverflow {
            gene: 1,
            sample: 1,
            column: 0,
            term: f64::MAX,
            value: f64::INFINITY,
        }
    );
}

#[test]
fn observation_weights_for_non_full_rank_design_use_zero_column_check() {
    let design = DesignMatrix::from_row_major(
        3,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 0.0,
        ],
        None,
    )
    .unwrap();
    let weights = RowMajorMatrix::from_row_major(
        2,
        3,
        vec![
            1.0, 1.0, 1.0, //
            1.0, 0.0, 1.0,
        ],
    )
    .unwrap();
    let output = preprocess_observation_weights(&weights, &design).unwrap();

    assert_eq!(output.design_rank, 2);
    assert_eq!(output.weights_fail, vec![false, true]);

    let rank_deficient = DesignMatrix::from_row_major(
        3,
        2,
        vec![
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        None,
    )
    .unwrap();
    let output = preprocess_observation_weights(&weights, &rank_deficient).unwrap();

    assert_eq!(output.design_rank, 1);
    assert_eq!(output.weights_fail, vec![false, false]);
}

#[test]
fn observation_weights_validate_inputs() {
    let design = two_group_design(4, Some("condition_B_vs_A"));
    let bad_dims = RowMajorMatrix::from_row_major(1, 3, vec![1.0, 1.0, 1.0]).unwrap();
    std::assert_matches!(
        preprocess_observation_weights(&bad_dims, &design).unwrap_err(),
        DeseqError::InvalidDimensions { context, expected: 4, actual: 3 }
            if context == "observation weight columns"
    );

    let bad_weight = RowMajorMatrix::from_row_major(
        3,
        4,
        vec![
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            -1.0,
            1.0,
            1.0,
            f64::NAN,
            1.0,
            1.0,
            1.0,
        ],
    )
    .unwrap();
    let error = preprocess_observation_weights(&bad_weight, &design).unwrap_err();
    assert_eq!(
        error.to_string(),
        "out-of-range value in observation weight at index Some(1): -1"
    );
    std::assert_matches!(
        error,
        DeseqError::InvalidCalculation {
            domain: ErrorKind::Numeric,
            operation: "observation weight",
            location: Some(CalculationLocation::GeneSample { gene: 1, sample: 1 }),
            value: -1.0,
            ..
        }
    );

    let weights = RowMajorMatrix::from_row_major(1, 4, vec![1.0; 4]).unwrap();
    std::assert_matches!(
        preprocess_observation_weights_with_options(
            &weights,
            &design,
            ObservationWeightOptions {
                weight_threshold: f64::NAN,
                ..ObservationWeightOptions::default()
            }
        ).unwrap_err(),
        DeseqError::InvalidParameter {
            domain: ErrorKind::Options,
            parameter,
            requirement: Requirement::FiniteAndNonNegative,
        } if parameter == "weight threshold"
    );
}

#[test]
fn builder_keeps_weight_threshold_shared_with_dispersion_options() {
    let observation_options = ObservationWeightOptions {
        weight_threshold: 0.25,
        ..ObservationWeightOptions::default()
    };
    let builder = DeseqBuilder::new().observation_weight_options(observation_options);

    assert_relative_eq!(
        builder
            .current_gene_wise_dispersion_options()
            .weight_threshold,
        0.25,
        epsilon = 1e-12
    );

    let dispersion_options = GeneWiseDispersionOptions {
        weight_threshold: 0.125,
        ..GeneWiseDispersionOptions::default()
    };
    let builder = builder.gene_wise_dispersion_options(dispersion_options);

    assert_relative_eq!(
        builder
            .current_observation_weight_options()
            .weight_threshold,
        0.125,
        epsilon = 1e-12
    );
}
