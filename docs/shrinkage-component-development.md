# Local shrinkage component development

This unreleased work extends the standalone Rust engine with normal-likelihood
uniform and truncated-normal mixtures and preparation of all six mixture-grid
families. R transports the same inputs to the same core. The local public
`lfcShrink(type="ashr")` route now uses the native engine for its default SQP
workflow, fixed/supplied priors and restricted likelihood expressions. Advanced
ashr options and the apeglm route remain unfinished; this is not a release-ready
replacement of the full optional-method interface.

Following the user's latest instruction, further development and commits stay
local. No new GitHub publication or release is authorized. Focused checks cover
changed numerical contracts during development; the full release, performance
and memory checks are deferred until the combined implementation is ready.
Previously recorded checks retain their original source scope.

## Implemented components

`UniformMixture` supports finite intervals, shifted and one-sided components,
zero weights and zero-width point masses. It exposes component log likelihoods,
posterior probabilities, means, SDs and a right-continuous CDF. Undefined
finite-arithmetic rows remain NaN instead of silently dropping components.
`TruncatedNormalMixture` supports finite/infinite bounds, half normals, shifted
components and zero SDs. `FixedMixture` shares family dispatch and batch
conditioning with the previously implemented normal mixture.

The compatibility rules are deliberate:

- Uniform likelihoods retain the original CDF orientation and subtraction.
  This preserves the original underflow/undefined rows and its moment rounding;
  matching those values is a parity claim, not an improvement in mathematical
  accuracy.
- Truncated-normal likelihoods retain the original supplied-component
  zero-centered convolution, including when supplied means differ from zero.
- Original truncated-normal SD evaluation recycles bound vectors across a
  transposed observation matrix. `posteriors()` retains the complete batch shape
  to reproduce it. Calling the scalar `posterior()` repeatedly represents
  separate one-row calls and can produce different SDs, as the original does.
- The independent truncated-normal CDF panel appends a zero-weight untruncated
  sentinel for component CDF evaluation. It avoids the original helper's
  uninitialized-variable error when no interval contains the threshold. The
  sentinel contributes no posterior mass. This panel certifies component values;
  it does not close the public method's error-decision contract.

`prepare_mixture_grid` builds normal, uniform, half-uniform, positive/negative
uniform and half-normal grids. It handles automatic geometric scales, supplied
scales, point masses, mode/range rules, duplicate constrained intervals,
Dirichlet penalties and preparation warnings. `MixtureGrid::fit` composes these
components with the existing standalone EM/SQP weight fitter. The grid oracle
executes original ashr preparation and intercepts only its downstream optimizer;
it does not claim new complete-method fitting parity.

## Numerical refinements

The initial normal log-tail continued fraction met its scalar ULP bound but
lost useful information when nearby endpoints entered posterior moments.
The current log-tail implementation uses Cody's published rational
approximations with a split square. It leaves the existing Wald tail kernels
unchanged. The original continued-fraction checkpoint remains documented at
[its recorded source](normal-log-tails-2026-11-10.md).

The normal log-density uses explicit fused multiplication/addition, matching the
recorded native ARM reference density rounding. R's vector expression-level
`exp` uses `1+x` close to zero; the existing formula rule now lives in a shared
Rust helper used by the component expressions. Distribution/GLM kernels keep
their own mathematical exponential. These details repair rejected far-tail and
narrow-interval trials without widening a numerical gate.

A separate original-R vector panel covers `exp`, `log1p` and `expm1` at signed
zero and both sides of their small-argument boundaries. Vector `log1p`/`expm1`
return zero at exact zero and the argument for nonzero magnitude at most machine
epsilon. The shared Rust expression helpers now preserve these rules. R's
byte-code scalar intrinsics can round differently; the reference evaluates
interpreted vectors, matching formula-frame expressions. No R calculation or
runtime is introduced into the Rust engine.

## Focused validation

The independent original-package panels use `deseq2-r462-native`: R 4.6.2 and
ashr 2.2-63. Existing frozen references and public-method bounds are preserved.

| Panel | Scope | Bound |
| --- | --- | --- |
| Uniform posterior | 10 priors, 650 observations, including 10 undefined rows; far tails, point masses, narrow/ultranarrow intervals | Mean/SD absolute `2e-6`; weights/CDF absolute `2e-8`; log likelihood `1e-8 + 1e-12 * abs(reference)`; exact NaN/infinity decisions |
| Truncated-normal posterior | 7 priors, 196 observations; batch SDs, half normals, shifted/equal bounds, zero SD/weights | Same component bounds and missingness decisions |
| Mixture preparation | 24 profiles across all six families; automatic/null, supplied/shifted and clipped grids | Numeric parameters `2e-14 + 2e-14 * abs(reference)`; exact dimensions, component order and warnings |
| Elementary expressions | 3 vector expressions with signed-zero and rounding-boundary inputs | Exact values and metadata bits |

Nine focused Rust tests pass, and the combined native source installs locally.
All 17 bounded-prior profiles (846 observations), 24 grids and 191 expression
cases pass direct native/R-free bit comparisons. Grid warnings also match and
R calculations are blocked during bounded/grid bridge calls. Root Clippy and
vendor integrity pass. The [local validation record](data/shrinkage_components_local_validation.json)
binds fixture, implementation and executable hashes.

Focused Rust checks cover these tables and existing formula-expression cases.
Direct R/native versus R-free comparisons use exact binary input transport,
require identical output bits, and block R calculations during native requests.
They remain distinct from comparison with original-package values. Detailed
local results are recorded in the private development journal and ignored output
directories; the final combined release audit will bind the accepted source and
artifacts after production adoption.

The earlier full uniform checkpoint passed 1,086 Rust tests, three backend
tests, 2,678 R assertions with R source check status OK, Clippy/MSRV checks in
both dependency graphs, and 47 Python tests (one platform skip). Those totals
precede the subsequent truncated-normal/grid/expression changes and do not
certify the current combined source. The combined development checks are focused;
there is no new total-workflow speed, memory or counted optimization-round claim.

## Reproduction

Use fresh output paths and the pinned original environment for the independent
reference generators:

```sh
Rscript scripts/generate_uniform_posterior_reference.R /tmp/uniform.json
Rscript scripts/generate_truncated_posterior_reference.R /tmp/truncated.json
Rscript scripts/generate_mixture_grid_reference.R /tmp/grids.json
Rscript scripts/generate_elementary_expression_reference.R /tmp/expressions.json
cargo test --locked -p rsdeseq2 --test bounded_posterior --test mixture_grid --test formula_expression
```

After building the standalone `evaluate_r_interface` example and installing the
candidate R package locally, use `check_uniform_posterior_identity.R` for both
bounded-family fixtures (the fourth argument selects the reference),
`check_mixture_grid_identity.R`, and `check_formula_expression_identity.R`.
Pass the same standalone executable, candidate library and a fresh output
directory. The public methods still require additional implementation and the
final combined checks described in the [migration plan](rust-only-wrapper-migration.md).


## Subsequent batch implementation in progress

Further **local work in progress** adds a complete location-scale batch engine:
normal/Student-t likelihood selection, default native weight fitting, supplied
and fixed priors, pruning, alpha scaling, excluded-row prior summaries, raw
unpruned CDF requests, sign/discovery probabilities and ordered s/q-values.
Twenty independently generated original-package batches cover automatic and
supplied grids, all six normal-likelihood families, observation weights,
missing/zero/infinite SEs, unnormalized priors and actual upstream rejection.
Their numerical values, NA/NaN distinctions, warnings and acceptance decisions
pass the focused Rust table. The subsequent public-adapter checks below cover
adoption of the default workflow, separately from the component table.

The user also approved restricted custom likelihood expressions translated into
Rust. The implementation imports one-expression function syntax and scalar
numeric/logical captures without executing R functions, active bindings or
promises. Arithmetic, comparisons, elementary math, `ifelse`, `pmin` and `pmax`
use the existing standalone expression evaluator. First moments without an
explicit expression use native adaptive Gauss integration; absent second moments
retain missing SDs, matching the original helper. Captured objects and array
constants are rejected before their methods can run.

Three original custom likelihood batches (logistic, shifted logistic and Laplace)
pass the same numerical and missing-kind checks. The Laplace panel also retains
all 65 original NaN-intermediate warnings. Native warning detection uses vector
expression semantics and a cached summary-operation graph; it preserves those
diagnostics without repeating the underlying statistical calculations. These
panels do not certify every custom expression, integration failure or advanced
optional-method contract.

A subtle original missing-row rule is retained: assigning a K-component prior
into M selected rows recycles its weights column by column. Normal conditional
means/SDs undergo the same permutation; uniform CDFs retain their interval order.
A non-coprime M/K panel checks repetition as well as permutation. Summary fallback
and raw CDF missingness are intentionally separate operations.

## Current local adapter checkpoint

The five Student-t moment failures are repaired. A standalone positive-Gamma
implementation uses Fullerton's public-domain SLATEC table; beta probability
evaluation retains representable complements near one and uses stable recurrences
and a compensated series for the affected shapes. Ordinary moment normalization
and log-tail normalization retain separate arithmetic schedules. All 1,476
bounded component rows (650 uniform, 196 truncated normal, 630 Student-t) now pass
the original bounds and unchanged oracles.

A new diagnostic panel exposed a limit of using R itself as a 16-ULP probability
oracle: some original tiny-tail beta probabilities differ from independent
90-digit truth by up to 274 ULP. The added leaf tests therefore check Gamma/beta
normalization against the original values, source CDF parity at `1e-13` relative,
and mathematical CDF accuracy at 16 ULP against an independent rational-beta and
integrated-binomial-series calculation. This changes no existing component or
public-method gate. The original diagnostic values remain in the repository.

The public ashr adapter passes the default, sign-s-value and positive-threshold
workflow on the existing 1,200-gene paired fixture. Maximum mean/SD differences
from original DESeq2 are below `5e-14`; the existing `1e-3` posterior and `0.02`
s-value bounds are unchanged. A rejected first adapter trial differed by
`0.02194` in mean because the general fitter checked all-coordinate KKT
convergence. Original mixSQP checks only positive coordinates. The compatibility
batch now explicitly selects that rule and the original default active-set
iteration limit; the general fitter retains its stronger convergence policy.
An independent 1,200-observation oracle and an analytical boundary example
protect both choices.

Threshold queries use the original result MLEs/SEs, a normal likelihood and
alpha zero, even when fitting inputs, alpha or likelihood are overridden. Ten
independent query panels cover that policy, recycling, missingness and
unnormalized priors. CDF-only requests avoid evaluating unused moment expressions.
The R importer reads already forced scalar captures through a read-only R 4.6
binding inspection entrypoint; active or unforced delayed bindings never run.
Older R retains the conservative literal/stored-value importer.

Actual freshly built bridge/R-free outputs match every bit for all 23 complete
location/custom batches, ten query batches and 630 Student-t component rows,
with R calculations blocked at the numerical boundary. The public workflow also
blocks R statistical calls, including original ashr fitting/CDF helpers, and
checks restricted logistic and canonical Student-t likelihoods. R object
presentation may still order labels and count storage sizes internally.

Still open: other built-in likelihood families, non-SQP optimizers, estimated
modes, remaining controls/output subsets, the full warning/integration-failure
matrix, standalone apeglm, and the combined authority/performance/memory/artifact
audit. Unsupported development options fail explicitly without an R fallback.
Full suites and release checks remain deferred until this work is combined.
Accepted checkpoint `153154a` and the published 0.4.0 release remain intact.

Reproduce the new focused checks after a local build/install:

```sh
cargo test -p rsdeseq2 --test bounded_posterior --test location_shrinkage \
  --test adaptive_result --test mixture_fit --test source_convergence
Rscript scripts/check_location_shrinkage_identity.R "$DRIVER" "$R_LIBRARY" "$FRESH_LOCATION_OUT"
Rscript scripts/check_adaptive_result_identity.R "$DRIVER" "$R_LIBRARY" "$FRESH_QUERY_OUT"
Rscript scripts/check_native_ashr_workflow.R "$R_LIBRARY" "$FRESH_WORKFLOW_OUT"
```

The `generate_adaptive_result_reference.R` generator uses only original ashr.
`generate_mixture_convergence_reference.R` reads immutable input bits and invokes
original ashr; its input may be the committed convergence fixture itself, with
a fresh output path. `generate_student_component_leaves.R` and
`generate_beta_mathematical_reference.py` reproduce the two independent leaf
panels. Candidate numerical outputs never enter those oracle calculations.
