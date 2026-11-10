# Numerical tests and parity checks

Local default ashr and negative-binomial apeglm migration checks are recorded
separately in [ashr development](shrinkage-component-development.md) and
[apeglm development](native-apeglm-development.md). They combine independent
original-package numerical references with direct native/R-free output-bit
comparisons and blocked R statistical calls. Advanced optional methods and the
final combined suite/performance/memory/artifact gates remain open.
The later [general-regression checkpoint](native-regression-development.md)
covers translated likelihoods and analytic `nbinomR`, with independent
high-precision density checks where original R loses tail precision.
The [beta-binomial checkpoint](native-beta-binomial-development.md) adds three
logit methods, saturated-row optimizer checks and explicit native warning
comparisons. Source warning counts and failure-reason wording can differ near
floating-point plateaus; these differences remain visible in the report.
The [random-start checkpoint](native-random-start-development.md) adds frozen
uniform/normal/state tables, 164 posterior profiles and public-session state
checks, including error handling and integer-NA state-word bits.

## Accepted benchmark runs

`benchmark_gtex_native_wald.py` defaults to one sequential pass. An accepted run
requires `--validated-results-root` and verifies that directory against the
reviewed `--validated-output-hashes` manifest (default: the v0.3.1 native-output
manifest). It also verifies the frozen 44 input/reference hashes. For example:

```sh
python3 scripts/benchmark_gtex_native_wald.py --repeats 1 \
  --validated-results-root results/benchmarks/release_0.3.1/gtex \
  --output-dir results/benchmarks/current-accepted
```

Use the documented pinned R environment. A finite but inaccurate candidate
cannot pass the native output hash check. Missing finite coverage, missingness
drift, altered binaries, changed baselines and mutated manifests fail the run.
`acceptance.json` is accepted only after successful completion. Companion stage
and diagnostic parity checks remain required when claiming preserved routes.

Explicit `--exploration` permits reporting without that acceptance baseline;
its completed report remains unaccepted for optimization/release evidence.
Historical timing files retain their original provenance and are not relabeled.

## Changes from the test-quality review

The numerical implementation and existing full reference files are unchanged.
The cleanup addresses failure detection and repeated setup:

- All 61 remaining DESeq2 reference tests require their committed fixtures.
  The shared reader rejects missing/empty files and duplicate columns. Shared
  numerical assertions reject unexpected infinities and NaNs.
- Two documentation-only checks were removed: a claim that native dispersion
  was future work, and a results-header-only check. Native dispersion and
  complete numerical result checks remain.
- Six R optimizer cases are checked in with a reproducible generator. They
  cover three real GTEx extremes in both compatibility and analytic-gradient
  modes, replacing a test that silently returned without a neighboring repo.
  Existing objective/parameter tolerances are retained.
- CLI failures check their intended error messages. Rlog option validation and
  all three pairs of incompatible beta-prior modes use shared table-driven
  fixtures. Weighted smoke tests now use nonuniform weights. Numerical routing
  tests compare size-factor, normalization-factor and observation-weight output
  with the library; weighted Wald and rlog must also differ from unweighted
  output. Independent DESeq2 oracles continue to check the library's arithmetic.
- Thread-count checks use verified one- and four-thread Rayon pools to exercise
  weighted and unweighted native Wald and LRT,
  including dispersions, coefficients, covariance, fitted means, optimizer
  diagnostics, convergence, Cook's output and complete result tables. Finite
  numerical fields must be identical; matching NaNs are checked explicitly.
- R ashr/apeglm assertions use per-gene absolute bounds and matching NA masks,
  with no infinities. Their method-specific tolerance constants are retained.
  The native bridge must be loaded, and a helper regression checks that a
  large gene cannot conceal an inaccurate tail gene.
- The current GTEx parity/timing/generation tools have failure-boundary tests:
  duplicate or missing genes/columns, malformed input, non-finite numbers,
  missingness drift, zero-error reduction, all-NA coverage, frozen hashes,
  binary/output mutation, memory units, stale counts and deterministic blocked
  and unblocked random splits. Prior frozen-gene scorer tests remain.

## CLI routing consolidation

The follow-up to `d2fdac0` replaces 44 CLI smoke/shape tests with 11 numerical
case groups: seven command-specific normalization tables and four Wald/LRT
contrast tables. These execute 68 explicit routing cases. Distinct invalid
input and ambiguous-alias tests remain. The CLI executable has 54 tests in total,
down from 87; the case tables retain the input routes and add numerical checks.

- Normalization tables compare every exported value, gene/sample label and
  missingness with direct library inputs. A control subset, frozen geometric
  means, supplied size factors, gene/sample normalization factors and weights
  must change numerical output relative to the default by more than 1e-6.
  Labeled input files deliberately have shuffled rows and columns. Both blind
  and design-aware transforms must match their respective library results and
  differ from one another.
- Contrast fixtures have eight samples and three coefficients with distinct
  effects. Requesting a nondefault coefficient therefore detects ignored
  selectors. Canonical names, cleaned reserved aliases, weighted numeric/list
  contrasts, factor levels and cleaned references check complete numerical
  result tables and relevant contrast metadata. Reserved-name fixtures repeat
  name-dependent routes; numeric-index coverage stays in the ordinary fixture.
- Shared assertions compare in-memory answers directly with exported cells,
  independently of the TSV writer. They reject unexpected non-finite cells and
  missingness differences. Existing independent R oracles validate the library's
  numerical implementation.
- The old thread-hint-only check is folded into the actual one/four-thread pool
  comparison, retaining size-factor estimation, base means and variance checks.
  The one-line optimizer-default flag assertion is folded into a behavioral
  test that verifies default fallback convergence and explicit disabling.

Focused files under `tests/cli/` separate numerical assertions, normalization
options and contrast options. Source plus helpers is smaller than the original
CLI suite. The replaced-test mapping, executed cases and validation are recorded
in `docs/data/cli_test_consolidation_2026-10-16_*`.

## R wrapper consolidation

The follow-up to `12b00fa` replaces 15 repetitive wrapper test blocks with five
behavior groups covering 18 explicit cases and 30 complete numerical result
tables. Including shared helpers, the R test suite shrinks from 4,678 to 4,136
lines. The result, Wald and LRT test files share small primitive fixtures;
object extraction and factor metadata now have focused test files.

- List fields, matrix mcols, S4 slots, model-matrix attributes and their aliases
  retain distinct cases. Unnamed count matrices still obtain gene labels from
  mcols; shuffled model-matrix rows still align to sample names. The nested S4
  assay/bracket-access case remains separate.
- Each table checks all six numerical fields at an absolute tolerance of 1e-12,
  gene/column labels, finite values and disabled independent-filter flags.
  Expected effects and standard errors are explicit; Wald tails and BH values
  come from base R. Coefficient alias fixtures request a nondefault coefficient
  whose answers differ from the default.
- Paired Wald/LRT cases share transformed/cleaned factor setup while retaining
  their original method-specific inputs, both reference alias spellings,
  coefficient names/weights and metadata. Direct constructors and object
  packaging remain distinct. All-zero contrasts check Wald's zero statistic
  and unit p-value, while LRT statistics/p-values remain supplied values.
  Selecting a coefficient by name on LRT mcols also preserves LRT inference.
- All 90 other test bodies and all 99 error assertions remain unchanged. Both
  files containing live DESeq2 comparisons are byte-identical, including the
  normalization, transformations, Cook's replacement and shrinkage parity
  checks. Existing numerical tolerances and frozen reference files are retained.

The source package passes `R CMD check` with 1,137 assertions and no failures,
warnings or skips. Runtime R/Rust sources and the native executable are
unchanged, so this consolidation does not change quantification precision,
speed or memory behavior. Test counts are not performance benchmarks.

The replaced-test mapping, explicit cases and validation are recorded in
[`data/r_wrapper_test_consolidation_2026-10-16_cases.json`](data/r_wrapper_test_consolidation_2026-10-16_cases.json)
and
[`data/r_wrapper_test_consolidation_2026-10-16_validation.json`](data/r_wrapper_test_consolidation_2026-10-16_validation.json).

## Remaining CLI options and R primitive routes

The follow-up to `c308904` removes eleven default-versus-`native=TRUE`
comparisons from normalization and Cook's tests. These invoked the same native
implementation twice. Explicit numerical expectations and the registered-DLL
bridge checks remain. The poscounts fixture now checks hand-computed values;
its answers were also checked against DESeq2 1.52.0 during this review.

Sixteen remaining CLI success/shape tests become six numerical groups with
31 explicit cases. Together with the previous routing tables, the CLI suite
now exercises 99 numerical cases in 44 test functions (previously 54).

- Three beta-prior groups cover explicit expanded, generated factor and
  generated additive designs. The oracle supplies explicit matrices, collapse
  groups, dispersions and base means directly to the library. Existing
  coefficient/contrast, normalization and weighting routes remain; additional
  cases select nondefault coefficients and use weighted positive/negative lists.
- Frozen rlog checks both blind modes, every transformed cell and shuffled
  intercept labels. Its inputs distinguish ordinary rlog, another prior variance
  and the other blind mode. Threshold tests distinguish both an ignored
  threshold and an ignored alternative; t-tail tests distinguish normal tails,
  scalar/residual degrees of freedom and shuffled nonuniform per-gene values.
- Wald/LRT filtering uses four replicates per group, an actual Cook's outlier
  and separated gene means. It checks complete results plus scalar, rejection
  and lowess exports. Forced independent filtering and a zero Cook's cutoff
  must change results. Active filtering checks alpha and theta numerically;
  disabled cases still check BH/missingness and metadata.
- Export comparisons retain 1e-12 absolute/relative tolerances, labels and
  both directions of missingness, rejecting unexpected nonfinite values.
  Counterfactual checks require differences above 1e-6 for the relevant options.
  All 38 other CLI test bodies, including validation and replacement sidecars,
  are byte-identical to the baseline.

Seven R packaging/list/formula-selector blocks become three groups checking
14 complete result tables. They retain supplied-SE packaging, chi-square LRT
tails, direct and packaged list routes, contrast-over-name precedence, original
formula coefficient names/weights, supplied LRT statistics and filtering
metadata. The R suite has 91 blocks (previously 95). All 99 error assertions and
both live DESeq2 parity files remain unchanged; the source package passes
`R CMD check` with 1,252 assertions and no failures, warnings or skips.

The full Rust suite passes 1,016 tests, with the existing full GTEx replay still
explicitly ignored; all 31 Python tests, formatting and Clippy pass. Runtime
R/Rust sources and the rebuilt release executable are unchanged. No new
quantification speed or memory benchmark is claimed for these test-only edits.
The mapping and validation are recorded in
[`data/test_consolidation_followup_2026-10-16_cases.json`](data/test_consolidation_followup_2026-10-16_cases.json)
and
[`data/test_consolidation_followup_2026-10-16_validation.json`](data/test_consolidation_followup_2026-10-16_validation.json).

## GTEx without a full fit on every run

The default fixed-dispersion replay uses **512 genes × all 61 samples**, selected
from 56,200 genes by reference/count tails and quantiles, missingness patterns
and deterministic background coverage. It checks all six result fields,
including adjusted p-values and both directions of missingness.

The panel uses full-cohort size factors and final dispersions. R recomputes BH
adjustment for the subset; taking original full-data adjusted p-values would be
incorrect. A separate default test still estimates normalization on all genes.
The panel is a fixed-dispersion regression, so the native dispersion, local-fit,
weighted, optimizer and full-contrast tests remain necessary.

The GTEx test executable was observed at **1.45 s**, including full normalization,
versus **11.46 s** in the preceding refactor validation (about **7.9× shorter**).
These are cached debug test-harness times, not a quantification speed benchmark.
The selected count table is 101,218 bytes versus 10,644,813 bytes for the full
matrix. Full data remain checked in for offline replay and contrast generation.

See the [fixture README](../crates/rsdeseq2/tests/data/gtex_v8_real_wald/README.md)
for selection details, provenance and regeneration commands.

## Commands and CI

```bash
cargo test --workspace
python3 -m unittest discover -s scripts/tests -p 'test_*.py'
cargo test -p rsdeseq2 --test gtex_real_wald_parity full_gtex_wald \
  -- --ignored --nocapture
```

Ordinary Rust CI runs the compact panel, full normalization, native/fixed R
oracles, optimizer replay and Python integrity checks. The manual Parity workflow
runs real frozen-oracle gates and the full 56,200-gene fit. The R workflow runs
live DESeq2 wrapper comparisons. The full replay is the single explicitly ignored
Rust case in routine runs; missing fixtures are never reported as passed cases.

Before accepting numerical or performance changes, also replay all eleven
biological/random full-data contrasts using the commands in the
[GTEx benchmark documentation](gtex-native-wald-benchmark.md). The compact panel
complements that release gate; it does not establish parity for every gene,
contrast or estimator by itself.

## Validation evidence

The CLI follow-up's validation is recorded in
[`data/cli_test_consolidation_2026-10-16_validation.json`](data/cli_test_consolidation_2026-10-16_validation.json).
It checks the complete Rust/Python suites, all eleven native GTEx contrasts,
exact release output hashes and executable identity. Runtime sources and the
existing frozen reference files are unchanged.

The initial cleanup's validation evidence is recorded in
[`data/test_quality_cleanup_2026-10-16_validation.json`](data/test_quality_cleanup_2026-10-16_validation.json).
It includes the executable hash, test results, full-reference gates and paired
speed/memory measurements. The executable is byte-identical to the baseline;
no numerical runtime code or R wrapper functions changed.

## Standalone batch checkpoint — 2026-10-30

The [batch checkpoint](rust-only-result-batches-2026-10-30.md) supersedes the
removed R-backed implementation's acceptance claims. Native R 4.6.2 validates
1,488 assertions with no failures, warnings or skips and source-package Status
OK. Rust passes 1,045 tests, with the existing opt-in full-GTEx case excluded
from the ordinary run; the separate eleven-contrast gate passes. Backend/Python
support suites pass one/45 tests. MSRV, root/backend Clippy/formatting, strict
documentation, crate verification and read-only vendor checks pass locally.
Remote platform jobs remain required for publication.

All 22 GTEx result/diagnostic exports remain byte-identical to `b0855d0`; 121
per-contrast and nine global precision gates pass. Direct standalone/wrapped
comparisons now cover 87 bit-equality cases, including mixed primitive batches
and 65-request public batches. R tests compare complete single/batch tables with
`identical()` and verify validation reuse, subset/edit handling, per-request
alpha/df/thresholds, Cook's masking, Wald/LRT, replacement and legacy states.
C ASan/UBSan and Rust ASan pass the malformed-call child and numerical controls;
instrumentation exclusions are recorded in the checkpoint.

BH/filtering, inference probabilities, contrast projection and contrast-zero
rules now use the standalone Rust engine. Remaining R calculations are listed in
the [migration inventory](rust-only-wrapper-migration.md); no reference or
precision bound has been relaxed to close a gate. Earlier checkpoint numbers and
measurements remain in their linked historical reports.

## Shared Cook/refit checkpoint — 2026-10-30

The [Cook/refit checkpoint](rust-only-cooks-policy-2026-10-30.md) passes 1,049
Rust tests, one backend test and 45 Python checks. The final R source package
has Status OK, 1,510 assertions, no failures/warnings/skips and 262 source files
matching the working tree. Direct standalone/wrapped comparisons cover 101
bit-equality cases. C ASan/UBSan and Rust ASan pass malformed-buffer and valid
shared/writable-storage controls. MSRV, formatting, Clippy, strict documentation,
vendor checks and crate verification pass. Remote publication jobs remain pending.

Independent row-policy checks cover the eligibility correction for additive or
mismatched designs, automatic stored formula context, exact thresholds, last
legacy ties, three-larger-count exceptions, missingness, replacement/refit behavior and
validation before mutation. Runtime guards reject R outer products/intermediate
BH in fitted presentation. All 22 GTEx exports remain byte-identical to
`a7b256f`, with 121+9 unchanged precision gates; those runs disable Cook's, so
the dedicated row-policy checks establish the corrected eligibility behavior.

One native paired timing pass shows +2.08% elapsed time and no optimization round.
All 18 paired R input/fit/output payloads are byte-identical. Targeted two-group
ten-result time improves 14.84% and peak RSS falls 27 MiB, while the full R
process sum changes +0.19%. Slower controls and remaining migration work are
retained in the checkpoint and completion audit.

A final installed-source audit also identifies remaining exact-tie and strict
formula-origin differences from real DESeq2. These checks preserve the legacy
wrapper tie rule for this migration; they do not establish full DESeq2 Cook's
parity for those edges. The checkpoint and completion audit keep the required
independent correction open.

The local [posterior analysis checkpoint](native-posterior-analysis-development.md)
adds forty original-apeglm profiles and sixteen fixed-mode grid queries. Native
and wrapped full analysis outputs agree exactly with R statistical calculations
blocked. Original mesh interval allowances and three deliberate source-bug
corrections are documented separately from unchanged MAP/SD/sign precision gates.
