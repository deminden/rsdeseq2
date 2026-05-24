# DESeq2 Gap Analysis

This page compares `rsdeseq2` with original Bioconductor DESeq2. The goal is
not to vendor or translate DESeq2 line by line; the goal is a Rust
implementation that matches DESeq2 behavior stage by stage.

R/DESeq2 is allowed only as an offline reference generator for tests and
benchmarks.

## Current Scope

`rsdeseq2` currently works best as a Rust library for primitive matrices and
explicit model matrices. It has substantial coverage for normalization,
fixed-dispersion GLM primitives, parts of dispersion estimation, and selected
Wald/LRT result paths.

The current CLI is deliberately narrow:

- `size-factors`
- `base-mean`

The richer pipeline exists in the Rust API and tests, not yet as a mature user
CLI.

## Matched Or Partly Matched

See [compatibility.md](compatibility.md) for a more detailed numerical parity
snapshot with fixture evidence. In short, the strongest current matches are
the normalization/base-metadata primitives, fixed-dispersion GLM/Wald/LRT
internals with supplied dispersions, implemented dispersion trend/MAP pieces,
and selected result/filtering/Cook's primitives.

| DESeq2 area | rsdeseq2 status |
| --- | --- |
| Count matrix shape | Implemented as genes x samples, row-major. |
| Size factors | `ratio`, `poscounts`, supplied geometric means, control genes, and supplied factors. |
| Normalized counts | Size-factor and gene/sample normalization-factor paths. |
| Base metadata | `baseMean`, `baseVar`, `allZero`, weighted base metadata. |
| Observation weights | Row-max normalization, design checks, `weights_fail`, and implemented weighted fixed/native paths. |
| Fixed-dispersion GLM | Intercept shortcut, IRLS, QR option, log likelihood, deviance, SEs, hats, beta-prior variance/refit primitives, Wald/LRT. |
| Wald tests | Selected coefficients, fixed-dispersion primitive contrasts, native linear-mu/GLM-mu primitive contrasts, threshold alternatives, normal and t p-values. |
| LRT | Fixed-dispersion full vs reduced and limited native-dispersion branches. |
| Cook's diagnostics | Cook's matrix, `maxCooks`, p-value masking, low-count heuristic primitive, replacement planning, limited Wald/LRT/contrast refit. |
| Independent filtering | BaseMean-driven filtered BH and DESeq2-shaped lowess threshold selection. |
| Dispersion gene estimates | Linear-mu and GLM-mu foundations, rough/moments starts, Cox-Reid, Armijo, grid fallback. |
| Dispersion trends | Parametric, mean, and initial pure-Rust local trends. |
| Transformations | `normTransform` plus mean-fit, parametric, and local VST primitives for normalized counts; VST dispatch from fitted dispersion trends. |
| Dispersion prior/MAP | Prior variance, MAP shrinkage, outlier replacement, weighted low-level objective pieces. |
| Diagnostics | `DeseqFit` fields plus DESeq2-style alias view for implemented row metadata. |
| Reference validation | Generated DESeq2 1.46.0 fixtures for implemented stages, including unweighted and weighted GLM-mu mean/local MAP/Wald/LRT anchors plus local Cox-Reid MAP/Wald/LRT anchors and compact result-row checks. |

## Major Missing Pieces

### End-To-End DESeq Parity

Original DESeq2 exposes high-level `DESeq()`, `results()`, `nbinomWaldTest()`,
and `nbinomLRT()` workflows over `DESeqDataSet`. `rsdeseq2` does not yet expose
a full equivalent end-to-end workflow. Current high-level Rust builder paths
cover only selected fixed-dispersion and limited native-dispersion branches.

### Full Dispersion Estimation

Still missing:

- exact local `locfit` edge-case parity,
- glmGamPoi trend and MAP behavior,
- exact seeded Monte Carlo/loess numerical identity for low-df dispersion prior
  variance,
- complete weighted dispersion parity beyond the current deterministic
  weighted GLM-mu mean/local fixtures,
- broader stage-by-stage native LRT references,
- more edge cases around convergence, skipped rows, and replacement refits.

### Full GLM Fitting

Still missing:

- expanded model-matrix beta-prior averaging and high-level workflow plumbing,
- broader validation of the new bounded optim fallback against DESeq2 rows that
  actually require backup fitting,
- complete weighted low-level `fitNbinomGLMs` behavior for rows DESeq2 marks
  `weightsFail` but still fits when called directly.

### Results And Contrasts

Still missing:

- full `results(contrast=...)` semantics,
- formula/metadata-aware factor-level handling,
- complete coefficient-name cleanup,
- complete formula-aware contrast Cook's/refit behavior,
- full Bioconductor result-object metadata and formula-aware result metadata
  beyond the typed primitive table view.

### Outliers And Refits

Still missing:

- full Cook's replacement-triggered refit for beta-prior and wrapper-object paths,
- all beta-prior interactions,
- full metadata preservation around replacement counts and final result tables,
- automatic formula-aware low-count Cook's heuristic selection.

### Transformations And Secondary APIs

Still missing:

- high-level VST object workflow and exact local `splinefun` parity,
- rlog,
- lfcShrink-compatible hooks,
- plotting helpers,
- mature CLI commands for full differential expression.

### R Wrapper

A mature R wrapper is planned after the Rust core reaches stronger parity. R
should prepare familiar inputs and present outputs, while Rust performs the
statistical computation.

Rules for that wrapper:

- no fallback to R/Bioconductor DESeq2 for runtime computation,
- unsupported paths fail explicitly,
- DESeq2 remains a test/reference dependency only,
- wrapper tests should compare wrapper output to Rust and generated DESeq2
  fixtures, not use DESeq2 as a hidden execution path.

## Reference Anchors

Primary DESeq2 source anchors used for parity work:

- `external/DESeq2/R/core.R`
- `external/DESeq2/R/results.R`
- `external/DESeq2/R/methods.R`
- `external/DESeq2/src/DESeq2.cpp`
- `external/DESeq2/tests/testthat/`

Important generated references live under:

- `crates/rsdeseq2/tests/data/deseq2_reference/`

The generator is:

- `scripts/generate_deseq2_references.R`

## Benchmark Interpretation

Current benchmarks compare only implemented apples-to-apples primitives, mainly
size-factor and base-mean CLI paths against equivalent DESeq2/R reference
operations. They are not claims about full `DESeq()` speed, because full DESeq2
workflow parity is not implemented yet.

The latest local runs use `/usr/bin/time -v` with repeated process-level runs.
On the 56,937 gene x 881 sample real count matrix, five-repeat medians were:

- `size-factors`: `rsdeseq2` 1.15 s and 199 MiB RSS versus DESeq2/R 26.71 s
  and 1.90 GiB RSS, max absolute difference `3.86e-14`.
- `base-mean`: `rsdeseq2` 1.38 s and 581 MiB RSS versus DESeq2/R 27.55 s and
  2.28 GiB RSS, max absolute difference `4.47e-07`.

See [benchmarks.md](benchmarks.md) for running time/RAM benchmarks.
