use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

fn two_group_design() -> DesignMatrix {
    DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap()
}

#[test]
fn observation_weights_are_normalized_by_gene_row_max() {
    let weights =
        RowMajorMatrix::from_row_major(2, 4, vec![2.0, 1.0, 0.0, 2.0, 0.5, 0.5, 0.0, 0.25])
            .unwrap();
    let output = preprocess_observation_weights(&weights, &two_group_design()).unwrap();

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
    let output = preprocess_observation_weights(&weights, &two_group_design()).unwrap();

    assert_eq!(output.weights_fail, vec![false, true, true]);
    assert_eq!(output.weights.row(2).unwrap(), &[0.0, 0.0, 0.0, 0.0]);
}

#[test]
fn observation_weights_use_thresholded_subset_for_cox_reid_check() {
    let weights = RowMajorMatrix::from_row_major(1, 4, vec![1.0, 1.0, 0.001, 0.001]).unwrap();
    let output = preprocess_observation_weights_with_options(
        &weights,
        &two_group_design(),
        ObservationWeightOptions {
            weight_threshold: 1e-2,
            ..ObservationWeightOptions::default()
        },
    )
    .unwrap();

    assert_eq!(output.weights_fail, vec![false]);
}

#[test]
fn observation_weights_reject_overflowed_cox_reid_column_sum() {
    let design = DesignMatrix::from_row_major(2, 1, vec![f64::MAX, f64::MAX], None).unwrap();
    let weights = RowMajorMatrix::from_row_major(1, 2, vec![1.0, 1.0]).unwrap();

    let err = preprocess_observation_weights(&weights, &design).unwrap_err();

    assert!(
        err.to_string()
            .contains("observation weight Cox-Reid column sum")
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
    let design = two_group_design();
    let bad_dims = RowMajorMatrix::from_row_major(1, 3, vec![1.0, 1.0, 1.0]).unwrap();
    assert!(preprocess_observation_weights(&bad_dims, &design).is_err());

    let bad_weight = RowMajorMatrix::from_row_major(1, 4, vec![1.0, -1.0, 1.0, 1.0]).unwrap();
    assert!(preprocess_observation_weights(&bad_weight, &design).is_err());

    let weights = RowMajorMatrix::from_row_major(1, 4, vec![1.0; 4]).unwrap();
    assert!(
        preprocess_observation_weights_with_options(
            &weights,
            &design,
            ObservationWeightOptions {
                weight_threshold: f64::NAN,
                ..ObservationWeightOptions::default()
            }
        )
        .is_err()
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
