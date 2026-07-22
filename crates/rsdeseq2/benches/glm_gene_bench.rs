use criterion::{Criterion, criterion_group, criterion_main};
use rsdeseq2::glm::nb::nbinom_log_likelihood;

fn glm_placeholder_bench(c: &mut Criterion) {
    c.bench_function("glm_placeholder", |bench| {
        bench.iter(|| nbinom_log_likelihood(&[10, 12, 20, 24], &[11.0, 13.0, 19.0, 23.0], 0.1))
    });
}

criterion_group!(benches, glm_placeholder_bench);
criterion_main!(benches);
