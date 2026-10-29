use criterion::{Criterion, criterion_group, criterion_main};
use rsdeseq2::prelude::*;
use std::hint::black_box;
mod support;

fn pipelines(c: &mut Criterion) {
    let data = support::fixture(2048, 16);
    let builder = DeseqBuilder::new()
        .threads(1)
        .fit_type(FitType::Mean)
        .disable_cooks_cutoff()
        .disable_independent_filtering();
    c.bench_function("normalization_2048_genes_16_samples", |b| {
        b.iter(|| {
            builder
                .fit_size_factors_and_base_means(black_box(&data.counts))
                .unwrap()
        })
    });
    c.bench_function("native_wald_mean_2048_genes_16_samples", |b| {
        b.iter(|| {
            builder
                .fit_with_results(black_box(&data.counts), &data.design)
                .unwrap()
        })
    });
    let transform = support::fixture(256, 24);
    c.bench_function("rlog_mean_256_genes_24_samples", |b| {
        b.iter(|| {
            builder
                .rlog_glm_mu(black_box(&transform.counts), &transform.design)
                .unwrap()
        })
    });
}
criterion_group!(benches, pipelines);
criterion_main!(benches);
