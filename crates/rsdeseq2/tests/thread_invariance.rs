use rsdeseq2::prelude::*;

#[test]
fn strict_builder_is_deterministic_for_initial_stages() {
    let counts = CountMatrix::from_row_major_u32(
        4,
        4,
        vec![10, 12, 20, 24, 0, 0, 5, 7, 100, 80, 90, 120, 3, 6, 9, 12],
    )
    .unwrap();

    let fit_a = DeseqBuilder::new()
        .execution_mode(ExecutionMode::Strict)
        .threads(1)
        .fit_size_factors_and_base_means(&counts)
        .unwrap();
    let fit_b = DeseqBuilder::new()
        .execution_mode(ExecutionMode::Strict)
        .threads(8)
        .fit_size_factors_and_base_means(&counts)
        .unwrap();

    assert_eq!(fit_a.size_factors, fit_b.size_factors);
    assert_eq!(fit_a.base_mean, fit_b.base_mean);
}
