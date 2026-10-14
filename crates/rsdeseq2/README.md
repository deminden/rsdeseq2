# rsdeseq2

`rsdeseq2` is a Rust implementation of DESeq2-compatible differential
expression workflows. It focuses on deterministic, inspectable execution for
normalization, dispersion estimation, GLM tests, transformations, and
result-table assembly.

`rsdeseq2` does not contain or reuse DESeq2 implementation code. It implements
the same documented and observed behavior independently, with parity checked
against reference outputs.

The crate includes both a Rust library and the `rsdeseq2` command-line tool.
It is not a drop-in replacement for the complete Bioconductor API.

## Install

Add the library to a Rust project:

```bash
cargo add rsdeseq2
```

Install the CLI:

```bash
cargo install rsdeseq2 --locked
```

## Supported Scope

Implemented areas include size-factor estimation, normalized counts and base
row metadata, fixed-dispersion and selected native-dispersion NB GLM Wald/LRT
workflows, DESeq2-style result contrasts, Cook's and independent-filtering
helpers, beta-prior refit workflows, `normTransform`/VST/rlog transformations,
and native CLI VST/rlog commands. The R package also provides dataset-based
VST/rlog objects and selected `lfcShrink()` methods; see its package README for
the supported Bioconductor integrations.

Unsupported interfaces include complete Bioconductor `DESeqDataSet` mutation
and metadata support, full glmGamPoi behavior, plotting, and other convenience
APIs.

Detailed status lives in
[docs/deseq2-gap-analysis.md][gap-analysis] and
[docs/compatibility.md][compatibility].

## Validation and Performance

Current validation covers 11 full native Wald contrasts from GTEx v8 Kidney
Cortex, including biological comparisons and reproducible blocked/unblocked
random splits. The [current parity, speed and memory comparison][gtex-benchmark]
records measurements, settings, hardware, and reproduction commands. The
[dispersion correction report][gtex-parity] documents the automatic mean-fitting
correction and remaining error bounds.

The test suite also covers documented normalization, GLM, Wald/LRT,
transformation, shrinkage, and R-wrapper paths. These measurements establish
behavior for the tested inputs and settings, rather than complete DESeq2
workflow parity.

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

Requires Rust 1.97.1 or newer and uses the Rust 2024 edition, as declared by
`rust-version` and `edition` in the workspace `Cargo.toml`.

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
python3 -m unittest discover -s scripts/tests -p 'test_*.py'
```

Generate DESeq2 reference fixtures:

```bash
Rscript scripts/generate_deseq2_references.R
cargo test -p rsdeseq2 --test dispersion_reference
cargo test -p rsdeseq2 --test wald_reference
cargo test -p rsdeseq2 --test lrt_reference
```

Run matched speed/RAM benchmarks for the supported primitives:

```bash
scripts/benchmark_rsdeseq2.sh
```

[gap-analysis]: https://github.com/deminden/rsdeseq2/blob/v0.3.0/docs/deseq2-gap-analysis.md
[compatibility]: https://github.com/deminden/rsdeseq2/blob/v0.3.0/docs/compatibility.md
[benchmarks]: https://github.com/deminden/rsdeseq2/blob/v0.3.0/docs/benchmarks.md
[gtex-parity]: https://github.com/deminden/rsdeseq2/blob/v0.3.0/docs/gtex-dispersion-improvement.md

[gtex-benchmark]: https://github.com/deminden/rsdeseq2/blob/v0.3.0/docs/gtex-native-wald-benchmark.md
