# rsdeseq2

> Experimental work in progress. `rsdeseq2` is not a drop-in replacement for
> DESeq2 and is not production-ready for differential expression analysis.
> Current development focuses on stage-by-stage parity for implemented
> normalization, dispersion, GLM, Cook's, and result-table primitives.

`rsdeseq2` is an early Rust implementation of the core DESeq2 workflow. The
goal is DESeq2-like intermediate and final results with a fast, inspectable
Rust numerical core and R-first compatibility.

This repository follows the same broad structure as `rsfgsea`: a Rust
workspace, a core crate, R package scaffolding, scripts, docs, benchmarks, and
reference/parity outputs. The statistical reference is the official
Bioconductor DESeq2 package. This project is a clean-room reimplementation:
DESeq2 behavior is documented and validated against, not copied line by line.

## Current Status

Currently implemented:

- Row-major `CountMatrix` for genes x samples count data.
- Generic row-major matrix storage.
- Basic `DesignMatrix` wrapper for R-generated model matrices, including
  deterministic full-rank checks for GLM-facing paths.
- `ratio` and DESeq2-style `poscounts` size-factor estimation.
- Optional caller-supplied size factors, supplied geometric means, and
  control-gene subsets for size factors.
- Size-factor normalized counts.
- Gene/sample normalization factors for normalized counts, base row metadata,
  supplied-dispersion fixed Wald/LRT GLM offsets, and the current native
  linear-mu dispersion/Wald subset.
- `baseMean`, `baseVar`, and `allZero` early row metadata, including
  DESeq2-style weighted base metadata helpers.
- Builder-owned observation weights for weighted base metadata, design-aware
  `weights_fail` flags, supplied-dispersion fixed Wald/LRT GLM paths, and the
  current GLM-mu native dispersion/Wald branch.
- Benjamini-Hochberg adjusted p-values with missing-value support.
- Negative-binomial log PMF and row/matrix log-likelihood helpers using
  DESeq2's `mu`/dispersion parameterization.
- Intercept-only fixed-dispersion NB GLM shortcut.
- Initial fixed-dispersion IRLS for supplied design matrices, with optional
  observation weights, per-coefficient natural-log-scale ridge values, and
  selectable normal-equation or DESeq2-style augmented QR solvers. GLM fit
  state exposes log likelihoods, DESeq2-style full deviance, beta convergence,
  and iteration counts.
- DESeq2-style observation-weight preprocessing helper: row-max normalization
  and design-rank failure flags.
- Default coefficient-level Wald statistics and p-values.
- Optional DESeq2-style Wald t p-values with residual, scalar, or per-gene
  degrees of freedom.
- Log2-scale beta covariance storage exposed in `DeseqFit` for implemented GLM
  fits and primitive linear-contrast Wald statistics.
- Result-row assembly for precomputed primitive numeric Wald contrasts.
- Primitive coefficient-name, positive/negative coefficient-list, and common
  factor-level contrast resolution for already-built design matrices.
- Selected-coefficient Wald LFC-threshold alternatives for `greaterAbs`,
  `greaterAbs2014`, `greaterAbsUPSHOT` without t p-values, `lessAbs`,
  `greater`, and `less`.
- Selected-coefficient Wald result rows with BH-adjusted p-values.
- Supplied-dispersion fixed-dispersion Wald pipeline for one coefficient and
  primitive numeric contrasts.
- DESeq2-style numeric/expanded `contrastAllZero` handling for primitive Wald
  contrasts where selected contrast samples can be identified from the design
  matrix.
- DESeq2-style character/factor-level `contrastAllZero` handling for primitive
  factor-level Wald contrasts when caller-supplied sample levels are available.
- Supplied-dispersion fixed-dispersion LRT pipeline for full vs reduced
  designs, including full deviance and full/reduced log-likelihood/convergence
  diagnostics in `DeseqFit`.
- Limited native-dispersion LRT pipeline for the current linear-mu and GLM-mu
  MAP dispersion branches.
- Initial linear-mu gene-wise dispersion estimator with DESeq2-style
  rough/moments starts, unweighted Cox-Reid objective scoring, Armijo line
  search, and grid fallback.
- Initial GLM-mu gene-wise dispersion estimator that alternates
  fixed-dispersion NB GLM mean fitting with fixed-mean dispersion optimization
  using DESeq2's `niter`/`fitidx` shape, with optional preprocessed
  observation weights.
- DESeq2-style parametric dispersion trend foundation:
  `dispersion = asymptDisp + extraPois / mean`, robust residual trimming, and
  Gamma identity-link IRLS.
- DESeq2-style mean dispersion trend: the same initial `100 * minDisp`
  viability gate, constant fitted dispersion from the `10 * minDisp`
  filtered trimmed mean, and `FitType::Mean` builder dispatch.
- DESeq2-style log-dispersion prior objective plus first and second derivative
  support for MAP dispersion fitting, including prior-aware line-search and
  grid optimizer entry points.
- Low-level observation-weighted dispersion objectives, weighted Cox-Reid
  adjustment first/second derivatives, and weighted prior-aware MAP optimizer
  entry points.
- Deterministic DESeq2-style dispersion prior variance estimation:
  MAD-squared log-residual variance, `trigamma((m - p) / 2)` sampling-variance
  subtraction, low-residual-df histogram/KL matching, and `0.25` floor.
- Initial DESeq2-style MAP dispersion fitting for the linear-mu and GLM-mu
  branches with `dispInit`, `log(dispFit)` prior means, prior-aware line
  search, grid fallback, `dispMAP`, `dispOutlier`, and final dispersion
  outputs. GLM-mu fit states can carry normalized observation weights into MAP.
- Limited native Wald pipeline for the current linear-mu no-weight and GLM-mu
  optionally weighted MAP dispersion paths, selected by `FitType::Parametric`
  or `FitType::Mean`.
- DESeq2-style all-zero row expansion for those fixed-dispersion pipelines.
- DESeq2-style Cook's distance matrix and `maxCooks` for that pipeline.
- Cook's cutoff p-value masking with BH recomputation for result rows.
- Explicit primitive helper for DESeq2's two-group low-count Cook's heuristic,
  for callers that know the design is a one-factor two-level formula case.
- Primitive Cook's outlier replacement-count transform with trimmed normalized
  means and size-factor or normalization-factor rescaling.
- Cook's replacement-refit planning metadata for replacement-count base
  statistics, `refitReplace`, `newAllZero`, and post-refit `maxCooks` masking.
- Limited Cook's replacement-refit execution for the implemented GLM-mu native
  Wald and LRT branches, preserving original size factors and merging refit rows.
- Base-mean independent filtering with filtered BH adjustment metadata and an
  R `stats::lowess`-shaped rejection-curve smoother for the default DESeq2
  threshold grid.
- DESeq2 `mcols(dds)`-style diagnostic alias view for implemented Wald/LRT
  fit-state fields.
- Inspectable `DeseqFit` skeleton and `DeseqBuilder`.
- Minimal CLI for size factors and base means.
- R package scaffold and R reference-generation script.
- R primitive helper for Cook's cutoff result masking and the explicit
  two-group low-count Cook's heuristic.
- Skip-safe Rust golden tests for generated DESeq2 normalization,
  observation-weight, fixed-dispersion GLM, and weighted GLM-mu native
  dispersion/Wald references.

Numerically reproduced against generated DESeq2 1.46.0 fixtures so far:
size factors (`ratio`, `poscounts`), normalized counts, `baseMean`/`baseVar`,
normalization-factor metadata, weighted base metadata, supplied-dispersion
unweighted Wald/LRT GLM beta/SE/log-likelihood/statistic/p-value fields,
fitted means, hat diagonals, Cook's distances, parametric dispersion trend,
dispersion prior variance, one MAP-dispersion fixture, and Cook's replacement
bookkeeping.

Not yet implemented:

- Full DESeq2 dispersion estimation, including complete weighted dispersion
  parity, local/glmGamPoi trend types, and production-ready end-to-end
  dispersion parity.
- Full negative-binomial GLM fitting parity, including DESeq2 beta-prior
  variance estimation and expanded-model handling, formula/colData-aware
  `results(contrast=...)` semantics, high-level integration of weight-failure
  flags, automatic formula-aware Cook's heuristics, contrast-aware Cook's/refit
  edge cases, and optim fallback.
- Full DESeq2 Wald and LRT parity beyond the current linear-mu/GLM-mu native
  Wald/LRT subsets.
- Full Cook's outlier replacement-triggered refit support for contrasts,
  beta priors, Bioconductor assay preservation, and all edge cases.
- DESeqDataSet integration.
- VST, rlog, lfcShrink, plotting, mature CLI, and Python bindings.

`rsdeseq2` is not production-ready for differential expression analysis until
stage-by-stage DESeq2 parity tests mature.

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

The CLI is intentionally narrow.

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

Input count files are tab-delimited with a header row, gene IDs in the first
column, and samples in remaining columns.

## R Usage

The current R package exposes primitive matrix helpers for the implemented
early normalization stages:

```r
library(rsdeseq2)

sf <- estimateSizeFactorsRust(counts, method = "ratio")
norm <- normalizedCountsRust(counts, sf)
baseMean <- baseMeanRust(counts, sf)
baseMetadata <- baseMetadataRust(counts, sf)

# Gene/sample normalization factors preempt size factors, matching DESeq2.
norm_nf <- normalizedCountsRust(counts, normalizationFactors = nf)
baseMetadata_nf <- baseMetadataRust(counts, normalizationFactors = nf, weights = weights)
```

These helpers currently use an R fallback that mirrors the Rust-supported
algorithms while the native Rust bridge is being wired. `baseMetadataRust()`
returns primitive `baseMean`, `baseVar`, and `allZero` row metadata and can
apply raw observation weights before row summaries. `estimateSizeFactorsRust()`,
`normalizedCountsRust()`, `baseMeanRust()`, and `baseMetadataRust()` also have
opt-in `native = TRUE` bridge attempts for these primitives, with R fallback
when the shared library is unavailable. `applyCooksCutoffRust(native = TRUE)`
can use the same bridge pattern for primitive Cook's masking while R keeps BH
adjustment and output assembly. Full DESeqDataSet integration is still future
work. The eventual high-level API should accept DESeqDataSet-like inputs and
return DESeq2-like outputs:

```r
library(rsdeseq2)

dds <- rsdeseq2::DESeq(
    dds,
    fitType = "parametric",
    test = "Wald",
    engine = "rust",
    parityMode = "strict",
    nproc = 8
)

res <- rsdeseq2::results(dds)
```

Unsupported integration still returns clear errors rather than pretending to be
a drop-in DESeq2 replacement.

## Development

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```

DESeq2 references can be generated with:

```bash
Rscript scripts/generate_deseq2_references.R
cargo test -p rsdeseq2 --test results_reference
cargo test -p rsdeseq2 --test wald_reference
cargo test -p rsdeseq2 --test lrt_reference
```

The generated references live under
`crates/rsdeseq2/tests/data/deseq2_reference/`. Without those files, the
golden-reference tests skip and the hand-computable Rust tests still run.

See `docs/` for algorithm notes, compatibility status, reproducibility, and the
development roadmap. The detailed implementation TODO list is maintained in
`docs/implementation-plan.md`.
