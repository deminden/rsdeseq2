# rsdeseq2

`rsdeseq2` is a Rust toolkit for DESeq2-compatible workflow
primitives. It focuses on deterministic, inspectable building blocks for
normalization, dispersion experiments, GLM tests, and result-table assembly.

`rsdeseq2` does not contain or reuse DESeq2 implementation code. It implements
the same documented and observed behavior independently, with parity checked
against reference outputs.

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
[docs/deseq2-gap-analysis.md][gap-analysis] and
[docs/compatibility.md][compatibility].

## Real-Data Parity

Current README benchmarks are shown only for outputs with matching reference
checks. A fresh real-data run compared `rsdeseq2` with DESeq2 1.46.0 on a
73,321 gene x 818 sample count matrix:

| primitive | parity check | rsdeseq2 | DESeq2 reference | speedup | rsdeseq2 RSS | reference RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `size-factors` | max diff `3.15e-14` | 3.48 s | 24.67 s | 7.1x | 237 MiB | 2.03 GiB |
| `base-mean` | max diff `5.47e-09` | 4.07 s | 25.88 s | 6.4x | 695 MiB | 2.47 GiB |

A five-tissue saved-reference sweep also matches offline DESeq2 outputs for
implemented primitive outputs. These rows are at floating-point parity for the
reported primitive:

| output | real-data coverage | max abs diff | max rel diff | mismatches | max RSS |
| --- | ---: | ---: | ---: | ---: | ---: |
| `size-factors` | 5 tissues, 1,998 samples | `2.62e-14` | `1.99e-14` | 0 | 237 MiB |
| `normalized-counts` | 5 tissues, 138,321,118 cells | `1.19e-07` | `9.74e-15` | 0 | 693 MiB |
| `base-mean` | 5 tissues, 341,286 genes | `4.66e-09` | `6.73e-15` | 0 | 694 MiB |

The current hard real-data Wald contrast is much tighter after the MAP
dispersion start fix, but still has visible tail differences in beta-fallback
and statistic rows:

| metric | mean abs diff | median abs diff | p99 abs diff | max abs diff | mismatches |
| --- | ---: | ---: | ---: | ---: | ---: |
| `baseMean` | `1.13e-12` | `8.88e-16` | `6.82e-12` | `6.52e-09` | 0 |
| `log2FoldChange` | `2.17e-08` | `3.77e-14` | `3.33e-12` | `7.70e-04` | 0 |
| `lfcSE` | `1.57e-10` | `2.33e-12` | `1.66e-10` | `8.26e-07` | 0 |
| `stat` | `3.19e-08` | `6.07e-12` | `3.44e-11` | `1.25e-03` | 0 |
| `pvalue` | `3.64e-09` | `3.03e-12` | `4.20e-11` | `6.50e-05` | 0 |
| `padj` | `2.12e-08` | `0` | `7.87e-11` | `4.50e-05` | 0 |

That focused contrast covers 65,580 genes and 78 retained samples; the latest
run took 151.0 s with 610 MiB peak RSS and zero swaps. These are validated
primitive CLI paths, not full-workflow `DESeq()` timings.
The local dispersion trend now uses
a pure-Rust locfit-compatible backend; on the same real-data fixture its 64,344
finite fitted values match DESeq2 with median relative error `3.74e-13`, p99
`5.85e-12`, and max `1.47e-11`. Existing committed GLM-mu local
MAP/Wald/LRT fixture metrics were already at machine precision and remain
unchanged.
Methodology and synthetic
benchmark results are in
[docs/benchmarks.md][benchmarks].

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

[gap-analysis]: https://github.com/deminden/rsdeseq2/blob/main/docs/deseq2-gap-analysis.md
[compatibility]: https://github.com/deminden/rsdeseq2/blob/main/docs/compatibility.md
[benchmarks]: https://github.com/deminden/rsdeseq2/blob/main/docs/benchmarks.md
