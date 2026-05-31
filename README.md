# rsdeseq2

`rsdeseq2` is an experimental Rust toolkit for DESeq2-compatible workflow
primitives. It focuses on deterministic, inspectable building blocks for
normalization, dispersion experiments, GLM tests, and result-table assembly.

Use it today as a Rust crate or narrow CLI for validated primitives. It is not
yet a drop-in replacement for full DESeq2 differential expression analysis.

## Current Scope

Implemented areas include size-factor estimation, normalized counts and base
row metadata, fixed-dispersion NB GLM Wald/LRT primitives, native dispersion
foundations, Cook's and independent-filtering helpers, beta-prior refit
primitives, and `normTransform`/VST/rlog building blocks.

Still in progress: full `DESeq()` workflow parity, formula-aware high-level
result handling, expanded-model beta-prior workflows, full glmGamPoi behavior,
high-level rlog object semantics, lfcShrink, plotting, and mature high-level
interfaces.

Detailed status lives in
[docs/deseq2-gap-analysis.md](docs/deseq2-gap-analysis.md) and
[docs/compatibility.md](docs/compatibility.md).

## Real-Data Parity

Current README benchmarks are shown only for outputs with matching reference
checks. A fresh real-data run compared `rsdeseq2` with DESeq2 1.46.0 on a
73,321 gene x 818 sample count matrix:

| primitive | parity check | rsdeseq2 | DESeq2 reference | speedup | rsdeseq2 RSS | reference RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `size-factors` | max diff `3.15e-14` | 3.48 s | 24.67 s | 7.1x | 237 MiB | 2.03 GiB |
| `base-mean` | max diff `5.47e-09` | 4.07 s | 25.88 s | 6.4x | 695 MiB | 2.47 GiB |

A five-tissue saved-reference sweep also matches offline DESeq2 outputs for
implemented primitive outputs:

| output | real-data coverage | harshest max diff | max RSS |
| --- | ---: | ---: | ---: |
| `size-factors` | 5 tissues, 1,998 samples | `2.62e-14` | 237 MiB |
| `normalized-counts` | 5 tissues, 138,321,118 cells | `1.19e-07` | 693 MiB |
| `base-mean` | 5 tissues, 341,286 genes | `4.66e-09` | 694 MiB |
| `wald_results` | 65,580 genes, 78 samples | median LFC diff `1.04e-13`; max lfcSE diff `3.27e-04`; max p-value diff `4.79e-05` | 610 MiB |

These are validated primitive CLI paths, not full-workflow `DESeq()` timings.
The latest sweep completed with zero swaps. The local dispersion trend now uses
a pure-Rust locfit-compatible backend; on the same real-data fixture its 64,344
finite fitted values match DESeq2 with median relative error `4.04e-10`, p99
`2.80e-09`, and max `3.19e-09`. Existing committed GLM-mu local
MAP/Wald/LRT fixture metrics were already at machine precision and remain
unchanged.
Methodology and synthetic
benchmark results are in
[docs/benchmarks.md](docs/benchmarks.md).

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
cargo run -p rsdeseq2 -- size-factors \
  --counts counts.tsv \
  --method ratio \
  --output size_factors.tsv

cargo run -p rsdeseq2 -- base-mean \
  --counts counts.tsv \
  --size-factors size_factors.tsv \
  --output base_mean.tsv

cargo run -p rsdeseq2 -- normalized-counts \
  --counts counts.tsv \
  --size-factors size_factors.tsv \
  --output normalized_counts.tsv

cargo run -p rsdeseq2 -- vst \
  --counts counts.tsv \
  --design design.tsv \
  --blind=false \
  --fit-type mean \
  --output vst.tsv

cargo run -p rsdeseq2 -- rlog \
  --counts counts.tsv \
  --design design.tsv \
  --blind=false \
  --fit-type mean \
  --output rlog.tsv

cargo run -p rsdeseq2 -- wald \
  --counts counts.tsv \
  --design design.tsv \
  --normalization-factors normalization_factors.tsv \
  --observation-weights observation_weights.tsv \
  --fit-type parametric \
  --coefficient 1 \
  --output results.tsv

cargo run -p rsdeseq2 -- lrt \
  --counts counts.tsv \
  --design design.tsv \
  --reduced-design reduced_design.tsv \
  --fit-type parametric \
  --coefficient 1 \
  --output lrt_results.tsv
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
