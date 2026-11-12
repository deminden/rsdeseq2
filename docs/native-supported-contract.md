# Standalone Rust and R wrapper contract

Version 0.5.0 uses one numerical engine. This inventory maps supported public
entry points and option families to their Rust owners and independent tests.
[Compatibility](compatibility.md) lists explicit limits; [testing](testing-parity.md)
records current passing CI.

## Numerical ownership

R validates and translates arguments, aligns identifiers, inspects stored
primitive metadata, and assembles R/Bioconductor objects. C validates and
marshals buffers. Statistical calculations belong to the standalone Rust core;
neither interface may select an R arithmetic implementation. Independent R
calculations belong in reference generators and tests only.

All 25 exports in the current R `NAMESPACE` fall into these groups. Paths in
the numerical-owner column are relative to `crates/rsdeseq2/src/`.

| Public R entry points | Numerical owner | Independent accuracy and interface evidence |
| --- | --- | --- |
| `DESeq`, `fitRust`, `estimateDispersionsRust` | `core/implementation/builder/`, `normalization.rs`, `dispersion/`, `glm/`; `core/replacement_fit.rs` merges replacement fits | Original DESeq2 fixtures and `check_fresh_fit_identity.R`: 31 matrix requests and 113 later queries, including weighted replacement state and execution decisions. |
| `results`, `resultsList`, `resultsNames`, `resultsNamesRust`, `resolveResultsContrastRust` | `contrasts/`, `core/fitted_results.rs`, `results/request.rs`, `glm/inference.rs`, `glm/inference_primitives.rs`, `results/cooks_policy.rs`, `results/replacement_policy.rs`, `independent_filtering.rs` | Original contrast/tail/filter references, primitive interface checks and complete later-query panels. Names and selector translation are adapter metadata; statistics use the shared kernels. |
| `waldFitRust`, `lrtFitRust`, `nbinomWaldTestRust`, `nbinomLRTRust` | Primitive fit packaging forwards inference to the same `glm/`, `contrasts/` and `results/` kernels | `test-primitive-routes.R`, `test-results-object-inputs.R`, and `check_rust_r_identity.R` exercise packaging, overrides, stored Wald/LRT dispatch and missing context. These constructors do not perform a new count fit. |
| `estimateSizeFactorsRust`, `normalizedCountsRust`, `baseMeanRust`, `baseMetadataRust` | `normalization.rs`, `normalization/offsets.rs`, `matrix/` | Original normalization/weighted metadata fixtures, `row_preparation.rs`, `results_reference.rs`, and primitive interface checks. Matrix normalization factors take precedence over sample factors. |
| `applyCooksCutoffRust` | `results/cooks_policy.rs` and shared BH adjustment | Original policy fixtures and `test-result-policy.R`; caller-supplied low-count context is validated, not inferred from a primitive table. |
| `resultsTableRust`, `applyIndependentFilteringRust` | `independent_filtering.rs`, `independent_filtering/lowess.rs`, shared result formatting | Original quantile/lowess/BH fixtures; `test-results-filtering.R` and primitive interface checks. R assembles already-computed columns and metadata. |
| `deseq2McolsDiagnosticsRust`, `rsdeseq2DiagnosticSchemaRust` | Shared diagnostics/schema and supplied primitive fit fields | Schema and primitive-route tests. These helpers expose available diagnostics; they do not reconstruct absent fits or implement every Bioconductor metadata convention. |
| `vst`, `varianceStabilizingTransformation`, `rlog` | `core/implementation/builder/stages/`, `transform/` | Original transform fixtures and integration tests; `check_frozen_local_identity.R` covers eight exact standalone/wrapped cases with R numerical calls blocked. |
| `lfcShrink` | `glm/beta/`, `shrinkage/`, shared formula/model-matrix construction | Original normal, ashr and apeglm references plus actual workflow/interface checks listed below. No production call to the original shrinkage packages computes the result. |

The R Rust modules under `src/rust/src/lib/` and `formula_request/` are
transport adapters. They construct the same core options and borrow or import
primitive storage. Shared kernels, rather than identical adapter source text,
are the authority. A lexical guard is useful evidence against accidentally
reintroducing R calls; it does not replace numerical or complete-workflow tests.

## Supported options and boundaries

- **Fitting and results:** Wald and full/reduced LRT, supported coefficient,
  numeric/list and factor-level contrasts, supplied normalization factors and
  observation weights, Cook's masking/replacement, independent filtering, and
  BH adjustment. Wald threshold alternatives and supplied t degrees of freedom
  use shared probability kernels. Later replacement queries require retained
  replacement context in the standalone API. Primitive fits reject requests
  whose required model/count context is absent.
- **Stored state and batching:** supported gene subsets map to the original
  fitted rows. Changed fitted sample inputs require a refit. Explicit `NULL`
  overrides are preserved. A batch validates shared inputs once and evaluates
  chunks of at most 64 requests; each request keeps its own masking/filtering
  decisions. Legacy objects without required offset or fit provenance may need
  the original input or a refit. See [wrapper behavior](r-wrapper.md).
- **Formula preparation:** declared native operators and built-in transforms,
  factor coding, polynomial groups, and offsets run in Rust. Complete formula
  fits use R term/variable order, independently checked without column remapping.
  The `expanded_formula_design*` compatibility helpers expose their
  explicit expanded layout; use `formula_design_for_fitting_from_model_frame`
  when preparing the complete-fit layout. Matching reduced offsets are accepted;
  conflicting offsets fail. Effective matrix normalization canonicalizes unused
  sample factors to ones for formula workflows. Explicit matrix workflows retain
  their supplied sample-factor metadata. See [formulas](r-wrapper.md#formulas).
- **Formula limits:** arbitrary custom functions and the full R formula grammar,
  including splines, are not supported. The standalone primitive frame has
  separate factor/numeric columns; ordered inputs use explicit `ordered()`.
  Arbitrary R classes and relative frame-column order are not claimed by this
  representation. No custom formula function is executed to discover values.
- **Transforms:** blind/design-aware VST and rlog, strict fast-VST subset
  selection, normalization/weights/offsets, frozen rlog inputs, saved
  parametric/mean curves, standard one-dimensional quadratic Gaussian saved
  locfit profiles, and restricted dispersion expressions. Unsupported saved
  profiles, custom calls, active bindings and unevaluated captured promises
  fail before execution. Existing finite original-reference bounds apply;
  interface identity does not imply bitwise identity to original DESeq2.
- **Normal shrinkage:** named/indexed coefficients and supported expanded factor,
  numeric/list contrasts. Numeric contrasts need supplied results; list-contrast
  requirements and ordered/transformed factor limitations are described in
  [the wrapper guide](r-wrapper.md). Formula interactions are rejected.
- **ashr:** all six declared normal-likelihood mixture grids; normal/t,
  canonical normal-mixture/log-F and Poisson/binomial likelihoods; and restricted
  likelihood expressions. Shared Rust handles mixture fitting, selected outputs,
  estimated modes, posterior moments and s-values. Numeric/named output selection
  and source warning sequences have independent references. See
  [ashr options and likelihoods](r-wrapper.md#ashr). Unsupported prior/likelihood
  combinations retain explicit errors, including combinations rejected by the
  original package deeper in its dispatch.
- **ashr controls:** SQP, `mixEM`, `w_mixEM`, `mixVBEM`, `cxxMixSquarem` and native
  `mixIP` use the [shared core control policy](r-wrapper.md#ashr).
  Active controls have the documented typed contract. A supplied fixed prior
  skips unused controls but still validates the optimizer name and required
  prior. Higher-order EM cycling and SQP spectral controls are supported as
  documented. Native `mixIP` is independently KKT checked; it is not a claim to
  reproduce an unavailable original Mosek execution. Arbitrary malformed R
  values accepted incidentally by an unused source branch do not establish an
  additional supported numerical method.
- **apeglm:** `general`, `nbinomR`, `nbinomCR`, `nbinomC`, `nbinomC*`, `betabinR`,
  `betabinCR`, `betabinC`, and `betabinC*`. Restricted general likelihoods,
  posterior grids/intervals, scalar parameter uncertainty and internal posterior
  contrasts follow the [analysis contract](r-wrapper.md#apeglm).
  Public `lfcShrink` keeps its selected-coefficient table shape; its `contrast=`
  selection argument remains unsupported for apeglm. BFGS/L-BFGS-B are supported
  when refinement requires an optimizer. Other refinement families are rejected;
  unused controls follow the independently checked
  [method-dispatch policy](r-wrapper.md#apeglm).
- **Random starts:** native Mersenne-Twister and L'Ecuyer-CMRG with Inversion
  normals import/write back state. Other R RNG configurations are rejected;
  there is no R random-generation fallback. See [state and failure behavior](r-wrapper.md#random-starts-and-interruption).
- **Execution:** reusable local worker pools serve fitting, refits and transforms.
  `strict`, `fast` and `diagnostic` currently select the same numerical algorithm.
  Public forwarding methods remain available. CLI publication uses staged
  outputs and a completion manifest; its failure/alias tests are separate from
  statistical correctness.

Full Bioconductor object semantics, arbitrary original fitted objects,
glmGamPoi fitting, every unstable low-level weighted ridge row, full R syntax,
and every R RNG/optimizer family remain outside the documented supported scope.
Unsupported inputs fail explicitly; there is no R numerical fallback.

## Formula and posterior resource contracts

Formula expansion removes repeated variable identities while preserving additive
groups, dot cross-products and term order. Input is limited to 1 MiB, parenthesis
nesting to 128 and generated interaction expansion to 16,384 terms. Excessive
requests fail before unbounded matrix/term construction.

A posterior retains the coefficient identity owning its probability summaries.
`summary_coefficient()` exposes it; `coefficient_summary()` rejects mismatched
coefficient requests and inconsistent matrix/vector shapes. Construct posterior
objects through fitting APIs; their identity field is private.

`PosteriorAnalysisOptions::estimated_grid_evaluations()` reports batch-wide
likelihood work, including parameter uncertainty. `PosteriorWorkControl` defaults
to 10,000,000 grid evaluations; `max_grid_evaluations=None` removes the budget.
`PosteriorCancellation` supports cooperative cancellation and an originating-thread
poll callback. R exposes the same budget as `max.grid.evaluations`, with explicit
`NULL` meaning unlimited. Work limits do not select lower precision.

## Validation

[Current CI](testing-parity.md#ci-status) passes both Linux R lanes, five target
panels, the full Rust suite, full GTEx replay and native-boundary instrumentation.
Independent original references and standalone/wrapped identity remain separate
requirements. [Reproducibility](reproducibility.md) identifies generators, frozen
inputs and acceptance gates. [Benchmarks](benchmarks.md) defines measured scope;
version changes do not turn saved timings into newly measured results.
