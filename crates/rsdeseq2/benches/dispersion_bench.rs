use criterion::{Criterion, criterion_group, criterion_main};
use rsdeseq2::prelude::*;

fn dispersion_linear_mu_bench(c: &mut Criterion) {
    let counts = CountMatrix::from_row_major_u32(
        4,
        4,
        vec![
            10, 12, 20, 24, //
            0, 0, 5, 7, //
            100, 80, 90, 120, //
            3, 6, 9, 12,
        ],
    )
    .unwrap();
    let size_factors = vec![1.0, 1.0, 1.0, 1.0];
    let normalized = normalized_counts(&counts, &size_factors).unwrap();
    let base_mean = base_mean(&normalized).unwrap();
    let base_var = base_variance(&normalized).unwrap();
    let all_zero = counts.all_zero_flags();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        None,
    )
    .unwrap();

    c.bench_function("dispersion_linear_mu", |bench| {
        bench.iter(|| {
            estimate_gene_wise_dispersions_linear_mu(
                GeneWiseDispersionInput {
                    counts: &counts,
                    design: &design,
                    size_factors: &size_factors,
                    normalization_factors: None,
                    normalized_counts: &normalized,
                    base_mean: &base_mean,
                    base_var: &base_var,
                    all_zero: &all_zero,
                    observation_weights: None,
                },
                GeneWiseDispersionOptions::default(),
            )
        })
    });
}

criterion_group!(benches, dispersion_linear_mu_bench);
criterion_main!(benches);
