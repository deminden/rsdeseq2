# rsdeseq2

`rsdeseq2` is a Rust implementation of DESeq2-compatible differential
expression workflows. It focuses on deterministic, inspectable execution for
normalization, dispersion estimation, GLM tests, transformations, and
result-table assembly.

`rsdeseq2` does not contain or reuse DESeq2 implementation code. It implements
the same documented and observed behavior independently, with parity checked
against reference outputs.

Use it today as a Rust crate, CLI, or R access layer for validated
DESeq2-compatible workflows.

## Current Scope

Implemented areas include size-factor estimation, normalized counts and base
row metadata, fixed-dispersion and native-dispersion NB GLM Wald/LRT workflows,
DESeq2-style result contrasts, Cook's and independent-filtering helpers,
beta-prior refit workflows, and `normTransform`/VST/rlog building blocks.

Interface work still in progress: complete Bioconductor `DESeqDataSet`
mutation/metadata plumbing, full glmGamPoi behavior, high-level rlog object
semantics, lfcShrink, plotting, and broader convenience APIs.

Detailed status lives in
[docs/deseq2-gap-analysis.md][gap-analysis] and
[docs/compatibility.md][compatibility].

## Real-Data Parity

README benchmarks are shown only for outputs with matching DESeq2 1.46.0
reference checks. The real-data parity sweep uses five publication-data tissues
for normalization outputs and one kidney Wald contrast for the full result-table
path.

| workflow | reference case | coverage | runtime / peak RSS | precision |
| --- | --- | ---: | ---: | --- |
| `size-factors` | five tissues | 1,998 samples | 1.55 s / 237 MiB | max abs `2.62e-14`, max rel `1.99e-14`, 0 mismatches |
| `normalized-counts` | five tissues | 138,321,118 count cells | 7.03 s / 693 MiB | max abs `1.19e-07`, max rel `9.74e-15`, 0 mismatches |
| `base-mean` | five tissues | 341,286 genes | 1.64 s / 694 MiB | max abs `4.66e-09`, max rel `6.73e-15`, 0 mismatches |
| `wald-results` | kidney Wald contrast `condition_B_vs_A` with design `~ perm_block + condition` | 65,580 genes, 78 samples | 151.0 s / 610 MiB | largest max abs `1.25e-03` (`stat`), 0 mismatches |
| `local-dispersion-trend` | publication-data local trend fixture | 64,344 finite fitted values | fixture check | median rel `3.74e-13`, p99 rel `5.85e-12`, max rel `1.47e-11` |

The Wald result row includes Cook's outlier replacement/refit, final Cook's
masking, and independent filtering; the full per-column Wald precision table is
in [docs/benchmarks.md][benchmarks].

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
