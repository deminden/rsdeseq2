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

## Numeric Evidence

| evidence | scope | result |
| --- | ---: | --- |
| Normalization | 17 tissues, 8,731 samples, 612.7M cells | zero finite/NA mismatches; max relative error `9.89e-15` |
| End-to-end Wald | 65,580-gene kidney contrast | LFC median / p99 / max abs error `2.46e-14` / `3.03e-12` / `7.70e-04`; p-value `2.37e-12` / `4.38e-11` / `6.50e-05` |
| L-BFGS-B isolation | 512 R 4.6.0 stress objectives | 0.2.0 matches 512/512 endpoints, values, and counts exactly |
| Process benchmark | 10k/50k genes × 16 samples | 33x–406x faster, 38x–100x lower peak RSS than DESeq2 1.52.0 for checked primitives |

The 0.2.0 optimizer is dramatically more precise in isolation. A
dependency-only replay moved end-to-end errors by less than 2%, but using its
analytic-gradient API reduced median/p99 LFC, SE, and statistic errors by about
20–22% versus finite differences. Only 26/65,580 kidney
genes (0.040%) and 305/535,178 fitted rows across eight real contrasts (0.057%)
used the fallback; tiny upstream dispersion differences still change these
sensitive optimizer targets. Three-run whole-workflow timing ranges overlapped,
so this is a precision improvement rather than a claimed speedup. See the repository [benchmark
documentation][benchmarks], tracked [before/after data][real-data-precision],
and [route/input-drift summary][real-data-routes].

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
[real-data-precision]: https://github.com/deminden/rsdeseq2/blob/main/docs/data/lbfgsb_real_data_precision.tsv
[real-data-routes]: https://github.com/deminden/rsdeseq2/blob/main/docs/data/lbfgsb_real_data_route_summary.tsv
