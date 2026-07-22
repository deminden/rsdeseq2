# DESeq2 Gap Analysis

This page compares `rsdeseq2` with original Bioconductor DESeq2. The goal is
not to vendor or translate DESeq2 line by line; the goal is a Rust
implementation that matches DESeq2 behavior stage by stage.

R/DESeq2 is allowed only as an offline reference generator for tests and
benchmarks.

## Implemented Scope

`rsdeseq2` provides a Rust library for primitive matrices and explicit model
matrices. Its implemented scope includes normalization, fixed-dispersion GLM
primitives, parts of dispersion estimation, and selected Wald/LRT result paths.

The CLI is deliberately narrower than the Rust API. It covers the main
implemented primitive workflows and can write result,
independent-filtering, Cook's, and replacement/refit sidecar tables for the
implemented Wald/LRT paths:

- `size-factors`
- `normalized-counts`
- `base-mean`
- `vst`
- `rlog`
- `wald`
- `lrt`

Broader object-style workflows are not exposed by the CLI. The Rust API
provides partial formula/model-frame support for implemented expanded-design
and contrast routes; it does not implement a full R formula engine.

## Matched Or Partly Matched

See [compatibility.md](compatibility.md) for a more detailed numerical parity
snapshot with fixture evidence. In short, the strongest matches are the
normalization/base-metadata primitives, fixed-dispersion GLM/Wald/LRT
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
| Wald tests | Selected coefficients, fixed-dispersion primitive contrasts, native linear-mu/GLM-mu primitive contrasts, threshold alternatives, normal and t p-values, including passing t-tail cases from DESeq2's `test_nbinomWald.R`. |
| LRT | Fixed-dispersion full vs reduced and limited native-dispersion branches. |
| Cook's diagnostics | Cook's matrix with diagnostic exports and CLI sidecars, `maxCooks`, p-value masking, low-count heuristic primitive, automatic two-level factor metadata condition for supplied-dispersion fixed and limited GLM-mu factor-level result routes plus limited GLM-mu replacement refits, replacement planning with scalar and row metadata plus assay exports, supplied-dispersion fixed Wald/LRT replacement refit, primitive expanded beta-prior Wald replacement refit, and limited native Wald/LRT/contrast refit. |
| Independent filtering | BaseMean-driven filtered BH and DESeq2-shaped lowess threshold selection. |
| Dispersion gene estimates | Linear-mu and GLM-mu foundations, rough/moments starts, Cox-Reid, Armijo, grid fallback. |
| Dispersion trends | Parametric, mean, and pure-Rust locfit-compatible local trends. |
| Transformations | `normTransform`, mean-fit, parametric, and local VST primitives, fit-state and builder-level VST dispatch, low-level rlog sample-effect fitting with retained GLM intermediates, frozen-intercept rlog primitives, rlog prior estimation, fit-state rlog and frozen-rlog reuse, builder-level design-aware/blind GLM-mu rlog and frozen-rlog helpers, and native CLI `vst`/`rlog` commands including frozen-intercept rlog input. |
| Dispersion prior/MAP | Prior variance, MAP shrinkage, outlier replacement, weighted low-level objective pieces. |
| Diagnostics | `DeseqFit` fields plus DESeq2-style alias view for implemented row metadata. |
| Formula/model-frame routing | Expanded formula designs can be built from explicit or builder-stored owned model-frame columns with factor reference inference/override, declared R factor level order, formula-local `relevel(factor, ref=...)`, `factor(...)`/`ordered(...)` identity, level, and label transforms, `as.factor(...)`/`as.ordered(...)` identity transforms, and `droplevels(...)` transforms, numeric covariates, formula offsets, and validation. Supported factor-transform string arguments are quote-aware for transform discovery, parenthesis matching, named-argument parsing, coefficient naming, and model-frame contrast metadata, including escaped quotes and backslashes. Formula design construction carries a formula-local model frame for supported derived factor/numeric columns; supported formula and model-frame fit objects retain it; and implemented Wald/LRT result contrast and supported Cook's replacement routes can resolve DESeq2-style character contrasts from explicit, builder-stored, or formula-local model-frame factor metadata, including exact-first R-cleaned aliases for factor names and levels plus reference inference from declared factor levels before observed sample order. Top-level Wald helpers can build supported formulas directly, and top-level LRT helpers can build supported full/reduced formulas directly; supported full-formula offsets are applied in fit/result, contrast-request, and Cook's replacement-refit helpers. |
| Reference validation | Generated DESeq2 1.46.0 fixtures for implemented stages, including unweighted and weighted GLM-mu mean/local MAP/Wald/LRT references plus local Cox-Reid MAP/Wald/LRT references and compact result-row checks. |

## Major Missing Pieces

### End-To-End DESeq Parity

Original DESeq2 exposes high-level `DESeq()`, `results()`, `nbinomWaldTest()`,
and `nbinomLRT()` workflows over `DESeqDataSet`. `rsdeseq2` does not expose a
full equivalent end-to-end workflow. High-level Rust builder paths cover only
selected fixed-dispersion and limited native-dispersion branches.

### Full Dispersion Estimation

Missing scope:

- broader synthetic `locfit` edge-case parity,
- glmGamPoi trend and MAP behavior,
- exact seeded Monte Carlo/loess numerical identity for low-df dispersion prior
  variance,
- complete weighted dispersion parity beyond the documented deterministic
  weighted GLM-mu mean/local fixtures,
- broader stage-by-stage native LRT references,
- more edge cases around convergence, skipped rows, and replacement refits.

The locfit-compatible local smoother mainly improves local `dispFit` parity.
On the versioned real-data local-trend fixture, median relative error improved
from `7.99e-03` to `3.74e-13`, p99 from `2.00e-01` to `5.85e-12`, and max from
`4.28e-01` to `1.47e-11`. The versioned downstream GLM-mu local MAP/Wald/LRT
fixtures show machine-precision parity, so local-smoother changes do not affect
those metrics. For the versioned high-error real-data Wald contrast, preserving
MAP dispersion line-search starts above `maxDisp` reduced the non-optimizer
lfcSE tail:
`lfcSE_max_abs` moved from `3.27e-04` to `8.26e-07`, while `lfcSE_mean_abs`
moved from `3.06e-08` to `1.57e-10`. Final MAP dispersions are bounded
before storage, matching the reference workflow shape.

### Full GLM Fitting

Missing scope:

- higher-level beta-prior integration around primitive expanded-model
  averaging and replacement refits,
- broader validation of the bounded limited-memory BFGS-style optim fallback
  against DESeq2 rows that actually require backup fitting,
- complete weighted low-level `fitNbinomGLMs` behavior for rows that DESeq2
  marks `weightsFail` but fits when called directly.

### Results And Contrasts

Missing scope:

- full `results(contrast=...)` semantics across every object-style route,
- full formula/metadata-aware factor-level handling beyond the implemented
  owned-model-frame and formula-local character contrast resolvers in supported
  fixed/native Wald/LRT routes,
- remaining coefficient-name cleanup beyond the implemented R-style aliases
  for primitive coefficient-name, list, and factor-level candidates,
- complete formula-aware contrast Cook's/refit behavior outside the
  implemented fixed and limited native Wald/LRT raw/model-frame/formula-local
  routes and formula-built top-level native helpers,
- full Bioconductor result-object metadata and formula-aware result metadata
  beyond the typed primitive table view.

### Outliers And Refits

Missing scope:

- full Cook's replacement-triggered refit for wrapper-object paths,
- complete Bioconductor assay and object metadata around beta-prior
  replacement refits,
- high-level metadata preservation around replacement counts and final result
  tables beyond the implemented primitive replacement/refit metadata and assay
  exports,
- broader formula-aware low-count Cook's heuristic selection outside the
  implemented supplied-dispersion fixed, limited GLM-mu two-level factor,
  and model-frame character contrast result/replacement-refit routes.

### Transformations And Secondary APIs

Missing scope:

- high-level VST object workflow and exact local `splinefun` parity,
- full high-level rlog object dispatch using frozen intercepts with complete
  object metadata,
- lfcShrink-compatible hooks,
- plotting helpers,
- mature CLI commands for full differential expression.

### R Wrapper

The R wrapper contract assigns familiar input preparation and output
presentation to R while Rust performs the statistical computation.

Rules for that wrapper:

- no fallback to R/Bioconductor DESeq2 for runtime computation,
- unsupported paths fail explicitly,
- DESeq2 remains a test/reference dependency only,
- wrapper tests should compare wrapper output to Rust and generated DESeq2
  fixtures, not use DESeq2 as a hidden execution path.

## DESeq2 Reference Sources

Primary DESeq2 sources used for validation:

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

Benchmarks compare matched primitives against
equivalent DESeq2/R reference operations, plus offline DESeq2 result fixtures
for broader real-data parity. They are not claims about full `DESeq()` speed,
because full DESeq2 workflow parity is outside the implemented scope.

The publication-style normalization sweep completed all 17 tissue matrices
with zero mismatches. Size factors matched to `2.132e-14`, normalized counts to
`1.94e-07`, and base means to `6.52e-09` max absolute difference. Median
per-task elapsed times in the six-worker sweep were 3.401 s for size factors,
25.231 s for normalized counts, and 3.571 s for base means; those absolute
times include concurrent CPU and filesystem contention.

IRLS fitted-mean and optimizer-fallback hat handling is validated on
high-error and diagnostic contrasts with full models and blocking factors for
heart, testis, and pancreas. These checks cover 278,257 result rows
with zero missing-row or finite/NA-pattern mismatches across `baseMean`,
`log2FoldChange`, `lfcSE`, `stat`, `pvalue`, and `padj`. The remaining
differences are finite numeric drift.

The L-BFGS-B fallback uses `rcompat-lbfgsb` 0.2.1. Against 512 bounded
negative-binomial stress objectives saved from R 4.6.1, it reproduced 512/512
endpoints, objective values, and function/gradient counts exactly. Version
0.1.6 matched endpoint plus objective in 493/512 cases and the objective alone
in 507/512 cases at the practical scan thresholds, 0/512 exactly, and 311/512
evaluation counts. These exact results apply to the recorded
x86_64 Linux, OpenBLAS 0.3.32, one-thread environment.

That isolated dependency-only gain did not translate into a large end-to-end
change on the controlled 65,580-gene kidney replay. For a separately measured
analytic-gradient configuration, median/p99 errors changed from
`3.15e-14` / `3.79e-12` to `2.46e-14` / `3.03e-12` for LFC, from
`1.77e-12` / `1.88e-10` to `1.38e-12` / `1.50e-10` for SE, and from
`5.33e-12` / `3.70e-11` to `4.18e-12` / `2.96e-11` for the Wald statistic
(about 22%/20% lower), while maximum errors remained unchanged. The
comparative measurements are recorded in
[`data/lbfgsb_real_data_precision.tsv`](data/lbfgsb_real_data_precision.tsv).

The final/refit implementation matches the R 4.6.1 callback arithmetic,
finite-difference trajectory, fused predictor accumulation, and
DESeq2-compatible natural-log QR starting values for log2-optimizer
initialization. It uses a dual-solver stability check that retains the analytic
endpoint only when the compatible solution is materially worse under
their common objective and the coefficient vectors disagree beyond numerical
resolution. Gene-wise dispersion remains analytic.

On the fixed 100-row high-error v0.2.4 set, the measured v0.2.5 median absolute
error is `6.064e-10`, mean error is `1.261e-4`, and maximum error is
`1.526e-3`. The v0.2.4 measurements are `1.464e-4`, `3.794e-4`, and
`3.094e-3`; using unrounded values, v0.2.5 is 241,402x lower at the median,
3.009x lower at the mean, and 50.68% lower at the maximum. Of the fixed rows,
89/100 improved and 78/100 improved by at least 10x. The dominant median error
was removed by compensated log accumulation in size-factor estimation, which
keeps near-boundary dispersion routing aligned with R's extended-precision
row-mean calculation.

Three additional full-blocked contrasts measured maximum Wald-statistic errors
of `1.768e-3` (heart), `5.615e-5` (testis), and `2.803e-4` (pancreas), with
zero missing-row or finite/NA-pattern mismatches.

L-BFGS-B was used by 26/65,580 kidney genes (`0.040%`) and 305/535,178 rows
across eight high-error real-data contrasts (`0.057%`). For eight kidney rows
that used L-BFGS-B and did not undergo replacement, median
dispersion drift of `4.77e-08` was amplified into median beta-target drift of
`9.82e-05`.
The stability check addresses this general sensitivity without reference data
or gene-specific exceptions. Optimizer-use and timing measurements are stored in
[`data/lbfgsb_real_data_route_summary.tsv`](data/lbfgsb_real_data_route_summary.tsv).

The 2026-07-22 moderate process benchmark used R 4.6.1 and DESeq2 1.52.0 on
10k/50k genes x 16 samples. Across the checked primitive CLI stages,
`rsdeseq2` used `0.010–0.140 s` and `6.01–16.54 MiB`, while the DESeq2/R process
used `3.865–4.625 s` and `601.64–638.58 MiB`. These absolute measurements
correspond to 33.04x–386.5x faster execution and 38.62x–100.05x lower peak RSS.

See [benchmarks.md](benchmarks.md) for running time/RAM benchmarks.
