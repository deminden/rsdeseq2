# rsdeseq2

`rsdeseq2` is an independent Rust implementation of DESeq2-compatible
normalization, dispersion estimation, negative-binomial Wald/LRT workflows,
transformations, and result assembly. It is available as a Rust crate, CLI,
and R access layer; it does not contain or call DESeq2 code at runtime.

## Evidence at a Glance

All precision figures compare with saved DESeq2/R 4.6.0 outputs. Speed figures
include process startup, parsing, computation, and output writing.

| evidence | data | result |
| --- | ---: | --- |
| Normalization | 17 tissues, 8,731 samples, 612.7M normalized cells | zero finite/NA mismatches; max relative error `9.89e-15` |
| End-to-end Wald | kidney full-blocked contrast, 65,580 genes | LFC median / p99 / max absolute error: `2.46e-14` / `3.03e-12` / `7.70e-04`; p-value: `2.37e-12` / `4.38e-11` / `6.50e-05` |
| Optimizer isolation | 512 bounded NB objectives | `rcompat-lbfgsb` 0.2.0: 512/512 exact endpoints, values, and counts; 0.1.6: 0/512 exact endpoints+values |
| Process benchmark | 10k/50k genes × 16 samples, DESeq2 1.52.0 | 33x–406x faster and 38x–100x lower peak RSS for checked primitive stages |

The dependency-only optimizer upgrade is a large isolated numerical
improvement but barely changes full real-data parity. Supplying the existing
closed-form beta gradient through the new 0.2.0 API does improve the production
path: versus finite differences, LFC median/p99 error fell by 22%/20%, SE by
22%/20%, and statistic by 22%/20% on the same 65,580-gene kidney replay. Only
26 genes (0.040%) entered L-BFGS-B; across eight real contrasts the rate was
305/535,178 (0.057%). For the non-replaced optimizer-tail rows, upstream
dispersion drift was amplified by the sensitive optimizer objective, so exact
L-BFGS-B cannot by itself reproduce a DESeq2 endpoint built from a slightly
different dispersion. Three whole-workflow runs showed overlapping timing
ranges and only a descriptive 0.8% median reduction, so no material speedup is
claimed for this rare route.

Tracked evidence: [controlled before/after errors][real-data-precision] and
[route, input-drift, and runtime summary][real-data-routes].
Commands, larger tables, and interpretation are in [benchmarks][benchmarks];
implemented and missing surfaces are in [compatibility][compatibility] and the
[gap analysis][gap-analysis].

## Scope

Implemented: ratio/poscounts normalization, normalized counts and row metadata,
native and supplied-dispersion NB GLMs, Wald/LRT results and contrasts, Cook's
handling, independent filtering, beta-prior refits, supported formula/model
frames, `normTransform`, VST, and rlog building blocks.

Still incomplete: full Bioconductor object mutation/metadata behavior,
glmGamPoi parity, high-level rlog semantics, `lfcShrink`, plotting, arbitrary R
formula expressions, and broader convenience APIs.

## Rust Usage

```rust
use rsdeseq2::prelude::*;

fn main() -> Result<(), DeseqError> {
    let counts = CountMatrix::from_row_major_u32(
        3,
        4,
        vec![
            10, 12, 20, 24,
            0,  0,  5,  7,
            100, 80, 90, 120,
        ],
    )?;

    let fit = DeseqBuilder::new()
        .size_factor_method(SizeFactorMethod::Ratio)
        .execution_mode(ExecutionMode::Strict)
        .fit_size_factors_and_base_means(&counts)?;

    println!("{:?}", fit.size_factors);
    println!("{:?}", fit.base_mean);
    Ok(())
}
```

## CLI

```bash
cargo run -p rsdeseq2 -- wald \
  --counts counts.tsv \
  --design design.tsv \
  --fit-type parametric \
  --coefficient 1 \
  --output results.tsv

cargo run -p rsdeseq2 -- --help
```

## Development

Requires current stable Rust, tracked in `Cargo.toml` via `rust-version`.

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```

Generate DESeq2 reference fixtures:

```bash
Rscript scripts/generate_deseq2_references.R
cargo test -p rsdeseq2 --test dispersion_reference
cargo test -p rsdeseq2 --test wald_reference
cargo test -p rsdeseq2 --test lrt_reference
```

Run speed/RAM benchmarks for current apples-to-apples primitives:

```bash
scripts/benchmark_rsdeseq2.sh
```

[gap-analysis]: https://github.com/deminden/rsdeseq2/blob/main/docs/deseq2-gap-analysis.md
[compatibility]: https://github.com/deminden/rsdeseq2/blob/main/docs/compatibility.md
[benchmarks]: https://github.com/deminden/rsdeseq2/blob/main/docs/benchmarks.md
[real-data-precision]: https://github.com/deminden/rsdeseq2/blob/main/docs/data/lbfgsb_real_data_precision.tsv
[real-data-routes]: https://github.com/deminden/rsdeseq2/blob/main/docs/data/lbfgsb_real_data_route_summary.tsv
