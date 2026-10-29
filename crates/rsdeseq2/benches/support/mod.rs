#![allow(dead_code)]
use rsdeseq2::prelude::{CountMatrix, DesignMatrix};

pub struct Fixture {
    pub counts: CountMatrix,
    pub design: DesignMatrix,
    pub size_factors: Vec<f64>,
    pub dispersions: Vec<f64>,
}

pub fn fixture(genes: usize, samples: usize) -> Fixture {
    let mut state = 20261029_u64;
    let mut counts = Vec::with_capacity(genes * samples);
    for gene in 0..genes {
        let base = (15 + gene % 500) as f64;
        for sample in 0..samples {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let noise = 0.5 + ((state >> 32) % 1000) as f64 / 1000.0;
            let fold = if sample >= samples / 2 && gene % 3 == 0 {
                1.5
            } else {
                1.0
            };
            counts.push((base * noise * fold) as u32);
        }
    }
    Fixture {
        counts: CountMatrix::from_row_major_u32(genes, samples, counts).unwrap(),
        design: DesignMatrix::from_row_major(
            samples,
            2,
            (0..samples)
                .flat_map(|s| [1.0, f64::from(s >= samples / 2)])
                .collect(),
            Some(vec!["Intercept".into(), "conditionB".into()]),
        )
        .unwrap(),
        size_factors: vec![1.0; samples],
        dispersions: (0..genes).map(|g| 0.05 + (g % 20) as f64 * 0.01).collect(),
    }
}
