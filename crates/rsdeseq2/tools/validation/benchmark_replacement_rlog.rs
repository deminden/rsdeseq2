//! Measure Cook's replacement and rlog prior estimation at fixed matrix sizes.
use rsdeseq2::prelude::*;
use std::{hint::black_box, time::Instant};
fn main() {
    let mut args = std::env::args().skip(1);
    let shape = args.next().unwrap_or_else(|| "20k".into());
    let (genes, samples) = match shape.as_str() {
        "20k" => (20_000, 100),
        "60k" => (60_000, 800),
        _ => panic!("expected matrix shape 20k or 60k"),
    };
    let threads = args
        .next()
        .map_or(1, |n| n.parse::<usize>().expect("positive worker count"));
    assert!(
        threads > 0 && args.next().is_none(),
        "expected a matrix shape and optional positive worker count"
    );
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let mut state = 20261119_u64;
    let values = (0..genes * samples)
        .map(|cell| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            12 + ((state >> 32) as u32 % 60) + ((cell / samples) as u32 % 200)
        })
        .collect::<Vec<_>>();
    let counts = CountMatrix::from_row_major_u32(genes, samples, values).unwrap();
    let design = DesignMatrix::from_row_major(
        samples,
        2,
        (0..samples)
            .flat_map(|i| [1., (i >= samples / 2) as u8 as f64])
            .collect(),
        None,
    )
    .unwrap();
    let factors = vec![1.; samples];
    let dispersions = vec![0.1; genes];
    let normalized = rsdeseq2::normalization::normalized_counts(&counts, &factors).unwrap();
    let means = rsdeseq2::normalization::base_mean(&normalized).unwrap();
    let cooks = RowMajorMatrix::from_elem(genes, samples, 0.5).unwrap();
    pool.install(|| {
        let replacement = CooksReplacementOptions::new(1.);
        let start = Instant::now();
        let fit = rsdeseq2::cooks::replace_outlier_counts(
            &counts,
            &normalized,
            &factors,
            None,
            &cooks,
            &design,
            &replacement,
        )
        .unwrap();
        let replacement_s = start.elapsed();
        let first = black_box(fit.candidate_replacement_counts.as_slice()[0]);
        drop(fit);
        let start = Instant::now();
        let prior = rsdeseq2::transform::estimate_rlog_sample_prior_variance(
            &normalized,
            &means,
            &dispersions,
        )
        .unwrap();
        println!(
            "{genes},{samples},{threads},{:.6},{:.6},{first},{:016x}",
            replacement_s.as_secs_f64(),
            start.elapsed().as_secs_f64(),
            prior.to_bits()
        );
    });
}
