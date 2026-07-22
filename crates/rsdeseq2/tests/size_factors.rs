use approx::assert_relative_eq;
use rsdeseq2::prelude::*;
use std::f64::consts::{FRAC_1_SQRT_2, SQRT_2};

#[test]
fn ratio_size_factors_hand_computable() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![2, 4, 8, 4, 8, 16]).unwrap();
    let size_factors = estimate_size_factors_ratio(&counts).unwrap();
    assert_relative_eq!(size_factors[0], 0.5, epsilon = 1e-12);
    assert_relative_eq!(size_factors[1], 1.0, epsilon = 1e-12);
    assert_relative_eq!(size_factors[2], 2.0, epsilon = 1e-12);
}

#[test]
fn poscounts_size_factors_with_zeros() {
    let counts = CountMatrix::from_row_major_u32(3, 3, vec![0, 4, 8, 5, 5, 5, 0, 0, 0]).unwrap();
    let size_factors = estimate_size_factors_poscounts(&counts).unwrap();
    assert_relative_eq!(size_factors[0], 1.0, epsilon = 1e-12);
    assert_relative_eq!(size_factors[1], 1.122462048309373, epsilon = 1e-12);
    assert_relative_eq!(size_factors[2], 1.5874010519681994, epsilon = 1e-12);
}

#[test]
fn ratio_size_factors_use_deseq2_log_scale_location() {
    let counts = CountMatrix::from_row_major_u32(2, 2, vec![1, 4, 4, 4]).unwrap();
    let size_factors = estimate_size_factors_ratio(&counts).unwrap();
    assert_relative_eq!(size_factors[0], FRAC_1_SQRT_2, epsilon = 1e-12);
    assert_relative_eq!(size_factors[1], SQRT_2, epsilon = 1e-12);
}

#[test]
fn size_factor_geometric_means_keep_large_counts_finite() {
    let counts =
        CountMatrix::from_row_major_u32(2, 3, vec![u32::MAX, u32::MAX, u32::MAX, 1, u32::MAX, 1])
            .unwrap();

    let ratio = estimate_size_factors_ratio(&counts).unwrap();
    let poscounts = estimate_size_factors_poscounts(&counts).unwrap();

    assert!(ratio.iter().all(|value| value.is_finite() && *value > 0.0));
    assert!(
        poscounts
            .iter()
            .all(|value| value.is_finite() && *value > 0.0)
    );
}

#[test]
fn original_norm_matrix_path_preempts_size_factors() {
    let counts = CountMatrix::from_row_major_u32(
        4,
        4,
        vec![
            2, 2, 2, 2, //
            4, 4, 4, 4, //
            6, 6, 6, 6, //
            8, 8, 8, 8,
        ],
    )
    .unwrap();
    let norm_matrix = RowMajorMatrix::from_row_major(
        4,
        4,
        vec![
            1.0, 1.0, 1.0, 1.0, //
            1.0, 1.0, 1.0, 1.0, //
            1.0, 1.0, 1.0, 1.0, //
            1.0, 1.0, 1.0, 1.0,
        ],
    )
    .unwrap();
    let true_size_factors = [2.0, 1.0, 1.0, 0.5];
    let normalization_factors = RowMajorMatrix::from_row_major(
        4,
        4,
        norm_matrix
            .as_slice()
            .chunks(4)
            .flat_map(|row| {
                row.iter()
                    .copied()
                    .zip(true_size_factors)
                    .map(|(norm, factor)| norm * factor)
                    .collect::<Vec<_>>()
            })
            .collect(),
    )
    .unwrap();

    let normalized = normalized_counts_with_factors(&counts, &normalization_factors).unwrap();

    assert_relative_eq!(normalized.row(0).unwrap()[0], 1.0, epsilon = 1e-12);
    assert_relative_eq!(normalized.row(0).unwrap()[3], 4.0, epsilon = 1e-12);
    for gene in 0..counts.n_genes() {
        for (sample, true_size_factor) in true_size_factors.iter().copied().enumerate() {
            assert_relative_eq!(
                normalization_factors.row(gene).unwrap()[sample]
                    / norm_matrix.row(gene).unwrap()[sample],
                true_size_factor,
                epsilon = 1e-12
            );
        }
    }
}

#[test]
fn supplied_geo_means_are_stabilized() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![2, 4, 8, 4, 8, 16]).unwrap();
    let size_factors =
        estimate_size_factors_ratio_with_options(&counts, Some(&[4.0, 8.0]), None).unwrap();
    let geometric_mean = (size_factors.iter().map(|value| value.ln()).sum::<f64>()
        / size_factors.len() as f64)
        .exp();
    assert_relative_eq!(geometric_mean, 1.0, epsilon = 1e-12);
}

#[test]
fn supplied_geo_means_length_is_validated() {
    let counts = CountMatrix::from_row_major_u32(4, 4, (1..=16).collect()).unwrap();
    let err = estimate_size_factors_ratio_with_options(&counts, Some(&[1.0, 2.0, 3.0]), None)
        .unwrap_err();
    assert!(err.to_string().contains("geometric means"));
}

#[test]
fn supplied_geo_means_all_zero_are_rejected() {
    let counts = CountMatrix::from_row_major_u32(4, 4, (1..=16).collect()).unwrap();
    let err = estimate_size_factors_ratio_with_options(&counts, Some(&[0.0, 0.0, 0.0, 0.0]), None)
        .unwrap_err();
    assert!(matches!(err, DeseqError::NoUsableGenesForSizeFactors));
}

#[test]
fn control_genes_subset_size_factors() {
    let counts =
        CountMatrix::from_row_major_u32(3, 3, vec![2, 4, 8, 100, 100, 100, 4, 8, 16]).unwrap();
    let size_factors = estimate_size_factors_ratio_with_options(&counts, None, Some(&[1])).unwrap();
    assert_relative_eq!(size_factors[0], 1.0, epsilon = 1e-12);
    assert_relative_eq!(size_factors[1], 1.0, epsilon = 1e-12);
    assert_relative_eq!(size_factors[2], 1.0, epsilon = 1e-12);
}

#[test]
fn control_gene_indices_are_validated() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![2, 4, 8, 4, 8, 16]).unwrap();
    let err = estimate_size_factors_ratio_with_options(&counts, None, Some(&[2])).unwrap_err();
    assert!(err.to_string().contains("control gene index"));
}

#[test]
fn builder_accepts_control_gene_mask() {
    let counts =
        CountMatrix::from_row_major_u32(3, 3, vec![2, 4, 8, 100, 100, 100, 4, 8, 16]).unwrap();
    let fit = DeseqBuilder::new()
        .control_gene_mask(vec![false, true, false])
        .fit_size_factors_and_base_means(&counts)
        .unwrap();
    assert_eq!(fit.size_factors, vec![1.0, 1.0, 1.0]);
}

#[test]
fn builder_validates_control_gene_mask_length() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![2, 4, 8, 4, 8, 16]).unwrap();
    let err = DeseqBuilder::new()
        .control_gene_mask(vec![true])
        .fit_size_factors_and_base_means(&counts)
        .unwrap_err();
    assert!(err.to_string().contains("control gene mask"));
}

#[test]
fn builder_uses_supplied_geometric_means() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![2, 4, 8, 4, 8, 16]).unwrap();
    let fit = DeseqBuilder::new()
        .geometric_means(vec![4.0, 8.0])
        .fit_size_factors_and_base_means(&counts)
        .unwrap();
    let geometric_mean = (fit.size_factors.iter().map(|value| value.ln()).sum::<f64>()
        / fit.size_factors.len() as f64)
        .exp();
    assert_relative_eq!(geometric_mean, 1.0, epsilon = 1e-12);
}
