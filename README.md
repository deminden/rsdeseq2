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
primitives, and `normTransform`/VST building blocks.

Still in progress: full `DESeq()` workflow parity, formula-aware high-level
result handling, expanded-model beta-prior workflows, full glmGamPoi behavior,
rlog, lfcShrink, plotting, and mature high-level interfaces.

Detailed status lives in
[docs/deseq2-gap-analysis.md](docs/deseq2-gap-analysis.md) and
[docs/compatibility.md](docs/compatibility.md).

## Real-Data Benchmark

Current README benchmarks are shown only for primitives with matching reference
outputs. On a real muscle raw-count matrix with 56,937 genes and 881 samples,
five process-level CLI runs gave these medians:

| primitive | parity check | rsdeseq2 | DESeq2 reference | speedup | rsdeseq2 RSS | reference RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `size-factors` | max diff `3.86e-14` | 1.15 s | 26.71 s | 23.2x | 199 MiB | 1.90 GiB |
| `base-mean` | max diff `4.47e-07` | 1.38 s | 27.55 s | 20.0x | 581 MiB | 2.28 GiB |

These are validated primitive CLI paths, not full-workflow `DESeq()` timings.
Methodology and synthetic benchmark results are in
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
  --method poscounts \
  --output base_mean.tsv
```

## Development

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
