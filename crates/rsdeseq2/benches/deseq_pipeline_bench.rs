use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rsdeseq2::prelude::*;

fn initial_pipeline_bench(c: &mut Criterion) {
    let n_genes = 2_000;
    let n_samples = 8;
    let values = (0..n_genes * n_samples)
        .map(|idx| ((idx % 101) + 1) as u32)
        .collect::<Vec<_>>();
    let counts = CountMatrix::from_row_major_u32(n_genes, n_samples, values).unwrap();

    c.bench_function("fit_size_factors_and_base_means", |bench| {
        bench.iter(|| {
            DeseqBuilder::new()
                .fit_size_factors_and_base_means(black_box(&counts))
                .unwrap()
        })
    });
}

criterion_group!(benches, initial_pipeline_bench);
criterion_main!(benches);
