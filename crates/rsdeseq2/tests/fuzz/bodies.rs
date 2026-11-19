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
        .as_chunks::<4>()
        .0
        .iter()
        .copied()
        .map(u32::from_le_bytes)
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
    let _ = DeseqBuilder::new().run(
        &Analysis::new(&matrix, &design)
            .fixed_dispersion(&vec![dispersion; genes])
            .result(ResultColumn::Coefficient(1)),
    );
}

/// Decode arbitrary formula transport bytes, then validate a canonical round trip.
pub fn formula_wire(data: &[u8]) {
    use rsdeseq2::design::expression::wire;
    let mut bytes = b"RFX2".to_vec();
    bytes.extend_from_slice(data);
    if let Ok(value) = wire::decode(&bytes) {
        let encoded = wire::encode(&value).expect("validated value can be encoded");
        let decoded = wire::decode(&encoded).expect("encoded value can be decoded");
        assert_eq!(wire::encode(&decoded).unwrap(), encoded);
    }
}

/// Parse byte-derived R seed words using an optional supported RNG header.
pub fn r_seed(data: &[u8]) {
    use rsdeseq2::shrinkage::PosteriorRng;
    let Some((&header, state)) = data.split_first() else {
        return;
    };
    let mut seed = match header % 3 {
        0 => vec![10403],
        1 => vec![10407],
        _ => Vec::new(),
    };
    seed.extend(
        state
            .as_chunks::<4>()
            .0
            .iter()
            .copied()
            .map(i32::from_le_bytes),
    );
    if let Ok(mut original) = PosteriorRng::from_r_seed(&seed) {
        let mut restored =
            PosteriorRng::from_r_seed(&original.r_seed()).expect("generated state is restorable");
        assert_eq!(original.uniform().to_bits(), restored.uniform().to_bits());
    }
}

/// Parse argument sets without dispatching a command or touching the filesystem.
pub fn cli_args(data: &[u8]) {
    use clap::Parser;
    let Ok(text) = std::str::from_utf8(data) else {
        return;
    };
    let mut fields = text.split('\0');
    let Some(selector) = fields.next() else {
        return;
    };
    let mut argv = vec!["rsdeseq2"];
    let command = match selector.bytes().next().unwrap_or(0) % 6 {
        0 => Some("wald"),
        1 => Some("lrt"),
        2 => Some("rlog"),
        3 => Some("size-factors"),
        4 => Some("normalized-counts"),
        _ => None,
    };
    argv.extend(command);
    argv.extend(fields);
    drop(rsdeseq2::cli::Cli::try_parse_from(argv));
}
