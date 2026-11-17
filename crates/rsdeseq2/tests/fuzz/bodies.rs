//! Bodies of the cargo-fuzz targets in `fuzz/`, shared with the stable smoke
//! test in `fuzz_smoke.rs`. Each body takes arbitrary bytes and must return
//! without panicking; errors are expected and ignored.

use rsdeseq2::design::formula_has_offset_terms;
use rsdeseq2::io::{
    read_count_matrix_tsv_from_reader, read_labeled_design_matrix_tsv_from_reader,
    read_labeled_sample_numeric_tsv_from_reader, read_sample_levels_tsv_from_reader,
};
use rsdeseq2::prelude::*;

/// TSV readers for counts, designs, sample values and sample levels.
pub fn tsv_readers(data: &[u8]) {
    let selector = data.first().copied().unwrap_or(0) % 4;
    let body = data.get(1..).unwrap_or(&[]);
    match selector {
        0 => drop(read_count_matrix_tsv_from_reader(body)),
        1 => drop(read_labeled_design_matrix_tsv_from_reader(body)),
        2 => drop(read_labeled_sample_numeric_tsv_from_reader(body, "fuzz")),
        _ => drop(read_sample_levels_tsv_from_reader(body)),
    }
}

/// Six samples: factors `condition` and `batch`, numeric `age` and `dose`.
pub fn model_frame() -> FormulaModelFrame {
    let factor = |name: &str, levels: [&str; 6]| FormulaFactorColumn {
        name: name.to_string(),
        sample_levels: levels.iter().map(|level| level.to_string()).collect(),
        levels: None,
        reference: None,
    };
    FormulaModelFrame {
        factors: vec![
            factor("condition", ["A", "A", "A", "B", "B", "B"]),
            factor("batch", ["x", "y", "z", "x", "y", "z"]),
        ],
        numeric_covariates: vec![
            FormulaNumericColumn {
                name: "age".to_string(),
                values: vec![30.0, 41.5, 52.0, 33.0, 47.0, 60.0],
            },
            FormulaNumericColumn {
                name: "dose".to_string(),
                values: vec![0.0, 0.0, 1.0, 1.0, 2.0, 2.0],
            },
        ],
    }
}

/// Formula parsing and design expansion against a fixed six-sample frame.
pub fn formula(data: &[u8]) {
    let Ok(text) = std::str::from_utf8(data) else {
        return;
    };
    let _ = formula_has_offset_terms(text);
    let _ = expanded_formula_design_with_offsets_from_model_frame(text, &model_frame());
}

/// A fixed-dispersion Wald fit on counts taken from the bytes: four samples
/// in two groups, any number of genes, any u32 count.
pub fn fixed_dispersion_wald(data: &[u8]) {
    const SAMPLES: usize = 4;
    let Some((&dispersion_byte, body)) = data.split_first() else {
        return;
    };
    let counts = body
        .chunks_exact(4)
        .map(|chunk| u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
        .collect::<Vec<_>>();
    let genes = counts.len() / SAMPLES;
    if genes == 0 || genes > 64 {
        return;
    }
    let Ok(matrix) =
        CountMatrix::from_row_major_u32(genes, SAMPLES, counts[..genes * SAMPLES].to_vec())
    else {
        return;
    };
    let design = DesignMatrix::from_row_major(
        SAMPLES,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        None,
    )
    .expect("fixed design");
    let dispersion = 1e-8 + f64::from(dispersion_byte) / 16.0;
    let _ = DeseqBuilder::new().fit_fixed_dispersion_wald(
        &matrix,
        &design,
        &vec![dispersion; genes],
        1,
    );
}
