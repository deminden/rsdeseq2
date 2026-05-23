# rsdeseq2

`rsdeseq2` is an experimental Rust implementation of DESeq2 workflow
primitives. It is not a drop-in replacement for DESeq2 and is not
production-ready for differential expression analysis.

The current product surface is the Rust crate, the Rust CLI, and an
experimental R package scaffold for wrapper development. Mature R wrapper
paths must call the Rust implementation and must not fall back to
R/Bioconductor DESeq2 for runtime computation.
R/DESeq2 is used only to generate offline parity fixtures and benchmarks.

## What Works

- Size factors: `ratio`, `poscounts`, supplied geometric means, control genes,
  and supplied size factors.
- Normalized counts, gene/sample normalization factors, `baseMean`, `baseVar`,
  `allZero`, and weighted base metadata.
- Fixed-dispersion NB GLM primitives, selected Wald/LRT paths, Wald
  thresholds, t p-values, primitive contrasts, Cook's distances, Cook's
  masking, outlier-replacement planning, and independent filtering.
- Linear-mu and GLM-mu native dispersion foundations with parametric/mean
  trends, prior variance, MAP shrinkage, selected native Wald/LRT paths, and
  observation-weight handling for the implemented GLM-mu branch.
- DESeq2 reference fixtures for stage-by-stage Rust tests, including weighted
  fixed Wald/LRT and weighted GLM-mu dispersion/MAP/Wald/LRT anchors.

## Still Missing

- Full DESeq2 end-to-end `DESeq()` parity.
- Local/glmGamPoi dispersion trends and glmGamPoi MAP behavior.
- Full beta-prior variance, expanded-model handling, optim fallback, and all
  contrast/refit edge cases.
- Full formula/metadata-aware result handling, automatic Cook's heuristics, VST,
  rlog, lfcShrink, plotting, mature CLI, and a mature Rust-backed R wrapper.

See [docs/deseq2-gap-analysis.md](docs/deseq2-gap-analysis.md) for the detailed
comparison with original DESeq2.

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

Benchmark details are in [docs/benchmarks.md](docs/benchmarks.md).
