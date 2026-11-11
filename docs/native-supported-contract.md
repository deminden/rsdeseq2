# Current standalone and R wrapper contract

This inventory describes local development after `6fea9b6`, following the
published 0.4.0 source. It maps the supported public surface to its numerical
owner and evidence. It is not a claim of complete DESeq2 compatibility or final
release acceptance. The original [C01–C23 plan](codebase-improvement-plan-2026-10-22.md)
and [one-engine requirement](rust-only-wrapper-migration.md) remain authoritative.

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
| `DESeq`, `fitRust`, `estimateDispersionsRust` | `core/implementation/builder/`, `normalization.rs`, `dispersion/`, `glm/`; `core/replacement_fit.rs` merges replacement fits | Original DESeq2 dispersion/Wald/LRT fixtures; [fresh fits](native-fresh-fit-development.md), [weighted state](weighted-replacement-state-correction.md), and [replacement execution](replacement-execution-correction.md). Actual interface checks cover 31 matrix requests and 113 later queries. |
| `results`, `resultsList`, `resultsNames`, `resultsNamesRust`, `resolveResultsContrastRust` | `contrasts/`, `core/fitted_results.rs`, `results/request.rs`, `glm/inference.rs`, `glm/inference_primitives.rs`, `results/cooks_policy.rs`, `results/replacement_policy.rs`, `independent_filtering.rs` | Original contrast/tail/filter references; [borrowed batches](rust-only-result-batches-2026-10-30.md), [Cook/refit policy](rust-only-cooks-policy-2026-10-30.md), and [later queries](native-replacement-query-development.md). Names and selector translation are adapter metadata; statistics use the shared kernels. |
| `waldFitRust`, `lrtFitRust`, `nbinomWaldTestRust`, `nbinomLRTRust` | Primitive fit packaging forwards inference to the same `glm/`, `contrasts/` and `results/` kernels | `test-primitive-routes.R`, `test-results-object-inputs.R`, and `check_rust_r_identity.R` exercise packaging, overrides, stored Wald/LRT dispatch and missing context. These constructors do not perform a new count fit. |
| `estimateSizeFactorsRust`, `normalizedCountsRust`, `baseMeanRust`, `baseMetadataRust` | `normalization.rs`, `normalization/offsets.rs`, `matrix/` | Original normalization/weighted metadata fixtures, `row_preparation.rs`, `results_reference.rs`, and primitive interface checks. Matrix normalization factors take precedence over sample factors. |
| `applyCooksCutoffRust` | `results/cooks_policy.rs` and shared BH adjustment | Original policy fixtures and `test-result-policy.R`; caller-supplied low-count context is validated, not inferred from a primitive table. |
| `resultsTableRust`, `applyIndependentFilteringRust` | `independent_filtering.rs`, `independent_filtering/lowess.rs`, shared result formatting | Original quantile/lowess/BH fixtures; `test-results-filtering.R` and primitive interface checks. R assembles already-computed columns and metadata. |
| `deseq2McolsDiagnosticsRust`, `rsdeseq2DiagnosticSchemaRust` | Shared diagnostics/schema and supplied primitive fit fields | Schema and primitive-route tests. These helpers expose available diagnostics; they do not reconstruct absent fits or implement every Bioconductor metadata convention. |
| `vst`, `varianceStabilizingTransformation`, `rlog` | `core/implementation/builder/dispersion_transforms.rs`, `transform/` | Original transform fixtures and integration tests; [frozen local fits and expressions](frozen-local-vst-audit-2026-10-31.md) include eight exact standalone/wrapped cases with R numerical calls blocked. |
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
  The older `expanded_formula_design*` compatibility helpers retain their
  historical explicit layout; use `formula_design_for_fitting_from_model_frame`
  when preparing the complete-fit layout. Matching reduced offsets are accepted;
  conflicting offsets fail. Effective matrix normalization canonicalizes unused
  sample factors to ones for formula workflows. Explicit matrix workflows retain
  their supplied sample-factor metadata. See [formula-fit scope](native-formula-fit-development.md).
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
  [continuous likelihoods](native-continuous-likelihood-development.md),
  [count likelihoods](native-count-likelihood-development.md),
  [mode](native-ashr-mode-development.md), and
  [warnings](native-likelihood-warning-development.md). Unsupported prior/likelihood
  combinations retain explicit errors, including combinations rejected by the
  original package deeper in its dispatch.
- **ashr controls:** SQP, `mixEM`, `w_mixEM`, `mixVBEM`, `cxxMixSquarem` and native
  `mixIP` use the [shared core control parser](native-fixed-prior-controls-development.md).
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
  contrasts follow the [analysis contract](native-posterior-analysis-development.md).
  Public `lfcShrink` keeps its selected-coefficient table shape; its `contrast=`
  selection argument remains unsupported for apeglm. BFGS/L-BFGS-B are supported
  when refinement requires an optimizer. Other refinement families are rejected;
  unused controls follow the independently checked
  [method-dispatch policy](native-regression-method-policy-development.md).
- **Random starts:** native Mersenne-Twister and L'Ecuyer-CMRG with Inversion
  normals import/write back state. Other R RNG configurations are rejected;
  there is no R random-generation fallback. See [state and failure behavior](native-random-start-development.md).
- **Execution:** reusable local worker pools serve fitting, refits and transforms.
  `strict`, `fast` and `diagnostic` currently select the same numerical algorithm.
  Public forwarding methods remain available. CLI publication uses staged
  outputs and a completion manifest; its failure/alias tests are separate from
  statistical correctness.

Full Bioconductor object semantics, arbitrary original fitted objects,
glmGamPoi fitting, every unstable low-level weighted ridge row, full R syntax,
and every R RNG/optimizer family remain outside the documented supported scope.
This inventory does not turn them into newly completed features or remove a
previously supported method. Source-bound reports retain their original scope.

## Remaining completion evidence

The [combined-suite audit](native-plan-audit-development.md) identified a default
ashr sparse SQP discrepancy in the existing 18-gene subset test. The subsequent
[shared Rust correction](native-sparse-sqp-development.md) closes that finding
and passes the unchanged original-comparison bound. Supported here describes
the intended API contract, not proof that every implementation route has
already passed final combined acceptance.

This inventory makes the remaining audit finite. C19/C21 still require a
requirement-by-requirement reconciliation of these entry routes, supported
options, warnings, state changes and forwarding adapters against current tests
and actual interface outputs. A collection of successful primitive tests alone
does not close the complete-workflow requirement.

C23 additionally requires the final combined source to pass frozen GTEx and
held-out precision gates, actual standalone/wrapped checks, current speed and
whole-process RSS measurements, Rust/R/Python checks, MSRV, boundary sanitizers,
source/package checks and required platform gates. The new formula-fitting
oracle has not yet been independently regenerated on Linux. Earlier timing and
package artifacts cannot certify later source changes. Rebuild and hash the
local release bundle only after reconciling these results. GitHub push, tagging
and publication require a new user instruction.
