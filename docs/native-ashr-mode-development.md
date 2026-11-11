# Native estimated prior modes

This unreleased checkpoint adds `mode="estimate"` to the shared native ashr
workflow. Standalone Rust exposes `estimate_location_mode` and `PriorModeSearch`,
returning the mode, profile log likelihood, bounds, actual fit count and ordered
warnings. The real R bridge calls the same function before final fitting.
Neither interface calls an R optimizer or likelihood callback.

The objective refits mixture weights at each candidate mode and maximizes the
unpenalized marginal log likelihood. Observation weights and Dirichlet parameters
affect the inner fit; the outer likelihood remains unweighted. Bounds use the
original effect estimates before alpha scaling, while the grid and likelihood
use the existing transformed observations. All supported component families,
Student-t/uniform likelihoods, restricted expressions and mixture optimizers
share the implementation.

The independent Rust search uses golden-section steps with inverse parabolic
interpolation. Its absolute tolerance is interval width times `EPS^(1/4)`, with
`sqrt(EPS)*abs(mode) + tolerance/3` point spacing, consistent with the
[documented scalar optimizer](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/optimize.html).
This is a bracketed local search; multimodal profiles have no global optimum
guarantee. The implementation does not copy or link R optimizer code.

The shared `location_fit` module separates validation/preprocessing and weight
fitting from posterior summaries. Search trials therefore do not calculate
posterior moments, CDF queries, pruning or sign-rate summaries. Existing complete
batch arithmetic retains its order. Original `optimize` performs an extra fit at
the selected point to populate its objective; Rust reuses the cached likelihood
and warnings. Reported native fit counts exclude that avoided duplicate.

## Bounds and rejected inputs

R's two-element numeric `mode` behavior is preserved explicitly: original ashr's
scalar `ifelse` keeps only the first coordinate, which expands the observation
range; the second coordinate does not restrict the search. A missing first
coordinate leaves the observation range unchanged. Frozen original traces cover
both orders and an interior pair. Two text elements also select the original
data-range search. For genuinely explicit bounds, standalone callers use
`PriorModeSearch::Interval(lower, upper)`.

Constant/nonfinite data-derived ranges, invalid explicit bounds, incompatible
supplied/fixed priors and nonfinite objective values reject clearly. Negative
standard errors retain the existing native rejection. Original ashr instead
emits NaN/replacement warnings and returns an endpoint for the negative-SE
profile; its unchanged observation is recorded separately from
`native_expected_error`. The existing safe omission of the original compiled
optimizer's empty-vector diagnostic warning also remains explicit.

## Local evidence

The independent generator runs only original ashr 2.2-63 under R 4.6.2. It records
26 full profiles and original bounds, search trials and likelihoods. Twenty-two
valid profiles cover automatic/supplied scales, all six component families,
Student-t and restricted logistic likelihoods, alpha scaling, ordinary/weighted/
variational/compiled-style optimizers, point masses, shifts and range controls.
Four invalid profiles retain original observations and explicit native policy.

Mode, posterior and likelihood comparisons retain `2e-6` bounds. Bounds and
search counts match the original traces, accounting for its duplicate final
fit. Warnings match for valid cases apart from the documented compiled warning.
Independent analytical tests recover the precision-weighted normal MLE and a
constrained optimum. An invalid posterior-moment expression leaves mode and
likelihood bits unchanged, proving it is not evaluated by the likelihood-only
search. Fifteen connected Rust tests and one scalar-search unit test pass; the
R posterior file passes 65 assertions.

All 26 complete R/native requests and independent R-free requests agree on every
returned bit and rejection decision. R statistical functions, including
`optimize`/`optimise`, are blocked with positive controls. Six public 128-gene
profiles cover normal/expanded searches, the wrapper's filtering-alpha argument,
EM, Student-t and restricted logistic fitting. Effect-scale alpha is a separate
primitive argument, covered above. Maximum mean/SD/s-value differences are
`3.61e-5`, `2.35e-6`, `6.31e-7`, within the unchanged public `1e-3`/`0.02` gates.

Twelve accepted complete payload sets replay exactly against immutable
`8093a48`, including prior output-selection, optimizer, regression and posterior
requests. Previously committed numerical oracles remain untouched. Core/backend
Clippy and vendor integrity pass. [Local evidence](data/native-ashr-mode-local.json)
binds the tested source, references and artifacts.

```sh
Rscript scripts/generate_ashr_mode_reference.R "$FRESH_REFERENCE_JSON"
cargo test -p rsdeseq2 --test ashr_mode --test location_shrinkage \
  --test adaptive_result --test ashr_output --test advanced_ashr_optimizers --test ashr_optimizers
cargo test -p rsdeseq2 --lib bounded_search_has_independent
cargo build --release --example evaluate_r_interface
Rscript scripts/check_native_ashr_mode.R "$DRIVER" "$R_LIBRARY" \
  "$FRESH_OUTPUT_DIR" "$PAIRED_FITS_RDS"
```

The cached paired-fit snapshot comes from the native random-start workflow; this
check uses its first 128 genes without repeating count fits. No runtime dependency,
quantification timing/RSS claim or counted optimization round is added. Remaining
built-in likelihood families, SQP controls, complete warning/error/authority and
final combined precision/speed/memory/sanitizer/platform/source-package gates
remain open. All development and commits stay local.
