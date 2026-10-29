use criterion::{Criterion, criterion_group, criterion_main};
use rsdeseq2::prelude::*;
use std::hint::black_box;
mod support;

fn dispersion_stages(c: &mut Criterion) {
    let data = support::fixture(512, 16);
    let builder = DeseqBuilder::new().threads(1);
    let normalized = normalized_counts(&data.counts, &data.size_factors).unwrap();
    let base_mean = base_mean(&normalized).unwrap();
    let base_var = base_variance(&normalized).unwrap();
    let all_zero = data.counts.all_zero_flags();
    c.bench_function("dispersion_linear_mu_512_genes_16_samples", |b| {
        b.iter(|| {
            builder
                .with_execution_context(|| {
                    estimate_gene_wise_dispersions_linear_mu(
                        GeneWiseDispersionInput {
                            counts: black_box(&data.counts),
                            design: &data.design,
                            size_factors: &data.size_factors,
                            normalization_factors: None,
                            normalized_counts: &normalized,
                            base_mean: &base_mean,
                            base_var: &base_var,
                            all_zero: &all_zero,
                            observation_weights: None,
                        },
                        GeneWiseDispersionOptions::default(),
                    )
                    .unwrap()
                })
                .unwrap()
        })
    });
}
criterion_group!(benches, dispersion_stages);
criterion_main!(benches);
