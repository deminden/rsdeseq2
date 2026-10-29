use criterion::{Criterion, criterion_group, criterion_main};
use rsdeseq2::prelude::*;
use std::hint::black_box;
mod support;

fn glm_stages(c: &mut Criterion) {
    let data = support::fixture(512, 16);
    let builder = DeseqBuilder::new()
        .threads(1)
        .size_factors(data.size_factors.clone())
        .fit_type(FitType::Mean)
        .disable_cooks_cutoff()
        .disable_independent_filtering();
    c.bench_function("irls_512_genes_16_samples", |b| {
        b.iter(|| {
            builder
                .with_execution_context(|| {
                    fit_fixed_dispersion_irls(
                        black_box(&data.counts),
                        &data.design,
                        &data.size_factors,
                        &data.dispersions,
                        IrlsOptions::default(),
                    )
                    .unwrap()
                })
                .unwrap()
        })
    });
    let (fit, _) = builder
        .fit_with_results(&data.counts, &data.design)
        .unwrap();
    let request = ResultsContrast::Numeric(vec![0.0, 1.0]);
    c.bench_function("fitted_contrast_512_genes_16_samples", |b| {
        b.iter(|| {
            fit.results_contrast_request(black_box(&data.counts), &request, None)
                .unwrap()
        })
    });
}
criterion_group!(benches, glm_stages);
criterion_main!(benches);
