# Numerical tests and parity checks

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
