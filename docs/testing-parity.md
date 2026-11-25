# Testing and numerical parity

Independent original R references check accuracy, and separate standalone and
wrapped runs check interface equality. The [0.6.0 development report](reports/hardening-progress.md)
records provenance for the local checks.

<a id="ci-status"></a>

## Validation status

The local 0.6.0 acceptance run checked source revision `a0808b2`.

| Check | Result |
| --- | --- |
| macOS arm64 release panel | 128 host R tests and 37 portable tests passed |
| Linux arm64 release panel | 145 host R tests and 37 portable tests passed |
| Golden outputs | 8,032 fields matched in each profile on both ARM targets |
| Full GTEx | Both profiles passed the 56,200-gene replay and all 130 eleven-contrast precision gates, with zero significance flips |
| R source package on macOS arm64 | R 4.5.3 and R 4.6.2 each passed 2,874 expectations with zero skips. Both passed native installation and SIGINT recovery checks |
| Rust source archives | Four archives reproduced byte for byte and passed 226 offline library tests |

The R 4.6.2 package check reported the Rust runtime's `_abort` import.
R 4.5.3 reported `Status: OK`. The [acceptance report](reports/hardening-progress.md#batched-acceptance-after-the-portable-arithmetic-changes)
records the checked artifacts and exact interface comparisons.

The [release procedure](release.md#ci-and-artifact-binding) defines the
hosted gates for publication. Native Intel macOS and Windows x64 runtime
checks are outside the local development acceptance above.

## Numerical coverage

- Independent normalization, dispersion, GLM, Wald/LRT, contrast, Cook's,
  filtering, transformation and shrinkage fixtures check both intermediates and results.
- Complete matrix/formula interface panels compare fits and later queries,
  including replacement rows, missing values, warnings, attributes and diagnostics.
- Formula tests cover term order, factors, offsets, expressions, polynomial
  groups, repeated variables and bounded combinatorial expansion.
- Shrinkage tests cover mixture fitting, posterior ownership, selected outputs,
  optimizer policy, likelihood expressions, far tails, random states and work limits.
- Invalid-input and CLI tests use numerical cases where dropped flags change
  results, plus separate process-failure, output-alias and completion-manifest checks.
- Frozen reference files are required. Missing, empty or malformed expectations
  fail; they never turn a numerical comparison into a silent skip.

Large Rust suites use one entry executable and focused child modules. R tests
share data/encoders while keeping independent original-package assertions.
A source-call audit supplements runtime checks with R statistical calls blocked;
it does not replace numerical comparisons.

## Golden output hashes

Refactoring that must not change results is checked with the
[golden harness](../crates/rsdeseq2/tools/validation/README.md#golden-output-harness):
one SHA-256 per output field over raw `f64` bits, for four inputs and the main
workflows at 1 and 4 threads. Rust CI runs both the host R and portable profiles on Linux x86_64.
Each release target runs the portable hash and bounded-posterior checks
alongside the host R numerical and CLI panels.

```sh
python3 scripts/references/fetch_pasilla_validation.py
cargo run --release --locked -p rsdeseq2 --example golden -- --check
```

## Compact and full GTEx checks

Ordinary tests fit a fixed-dispersion panel of **512 genes × all 61 samples**,
selected from count/reference tails, quantiles, missingness and deterministic
background coverage. It checks all six result fields and both directions of
missingness. Full-cohort size factors and final dispersions are supplied; the
subset's BH-adjusted p-values are recomputed for that subset. The compact
and full fixed-dispersion replays both enforce a log10 error of at most
`1e-5` for p-values and adjusted p-values.

The opt-in full replay covers all genes and independently checks their
normalization. Release CI runs it under the host R and portable profiles.
Restore the frozen full tables by SHA-256 from Git before running the full
replay or generating new contrasts. See the
[fixture description](../crates/rsdeseq2/tests/data/gtex_v8_real_wald/README.md).

```sh
cargo test --locked --workspace
python3 -m unittest discover -s scripts/tests -p 'test_*.py'
python3 scripts/references/stage_gtex_full_fixture.py --fetch
for profile in r-host portable; do
  RSDESEQ2_NUMERIC_PROFILE="$profile" cargo test --release --locked -p rsdeseq2 \
    --test integration gtex_real_wald_parity::full_gtex_ -- --ignored --nocapture
done
```

Significance agreement on the eleven contrasts (log10 errors of small p-values,
Spearman correlation of the statistic, top-k overlap and significance flips,
which fail the check) comes from
`python3 scripts/validation/gtex_significance.py --output-dir results/checks/significance`;
the compact panel test also fails on any flip at padj 0.05 or 0.1.

The eleven full biological/random contrasts provide a separate acceptance panel;
see [reproducibility](reproducibility.md#eleven-contrast-acceptance).

## R package and interface checks

Use the pinned environment and its restored dependencies. From the repository root:

```sh
mkdir -p results/checks/r
scripts/release/build_r_source.sh results/checks/r
(cd results/checks/r && \
  R CMD check --no-manual --no-build-vignettes rsdeseq2_0.7.0.tar.gz)
cargo build --release --locked -p rsdeseq2 --example evaluate_r_interface
Rscript --vanilla scripts/validation/check_rust_r_identity.R \
  target/release/examples/evaluate_r_interface results/checks/r/rsdeseq2.Rcheck \
  results/checks/primitive-identity
Rscript --vanilla scripts/validation/check_formula_fit_identity.R \
  target/release/examples/evaluate_r_interface results/checks/r/rsdeseq2.Rcheck \
  results/checks/formula-identity
RSDESEQ2_TEST_RSCRIPT="$(command -v Rscript)" \
RSDESEQ2_TEST_RLIB="$PWD/results/checks/r/rsdeseq2.Rcheck" \
  python3 -m unittest discover -s scripts/tests -p test_profile_r_workflows.py
```

Use fresh output directories. The [R workflow](../.github/workflows/r.yml) lists the complete
matrix/fresh-fit, expression, factor, frozen-local, mixture, posterior, warning
and live-reference identity panels; it is the executable combined gate.

The Windows x64 job installs the source archive checked by the Linux
reference lane with Cargo offline. It runs
`scripts/validation/check_r_native_install.R LIBRARY` to check registered
routine arities, normalization, recovery after a rejected count buffer and
agreement between one- and two-worker fits. It also checks estimated and
frozen rlog with per-cell factors and an all-zero row against fits on the
active rows.

## Platform arithmetic and native safety

Exact formula panels are compared under the macOS R profile and the portable
profile. Fresh SQP and uniform-posterior outputs remain separate from committed fixtures. Windows and
Intel macOS generate independent uniform CDFs at the existing bounds; other
posterior fields still use their original fixtures. See [platform arithmetic](platform_arithmetic.md).

```sh
python3 scripts/validation/check_posterior_interruption.py \
  --rscript "$(command -v Rscript)" --r-library results/checks/r/rsdeseq2.Rcheck
python3 scripts/validation/check_native_boundary_sanitizers.py \
  --rscript "$(command -v Rscript)" --r-library "$R_LIBS_USER" \
  --output-dir results/checks/native-boundary
```

Sanitizers instrument C with ASan/UBSan and Rust with ASan, including malformed
calls and replacement/formula-owner forced unwinds. R, dependencies and Rust's
standard library remain uninstrumented; leak detection is disabled. Compiler
instrumentation is verified from fingerprints, including Cargo's current layout.
Actual SIGINT, callback thread affinity, shared worker cancellation and reuse
after cancellation have separate tests. Grids poll cancellation every 256
likelihood calls; R polling is throttled to 50 ms and contained by `R_ToplevelExec`.

## Performance acceptance

Use [matched benchmarks](benchmarks.md) after numerical changes. Accepted runs
require frozen input/reference hashes and reviewed native result/diagnostic
exports before timing. Finite wrong values, changed gene coverage, missingness
drift or altered exports fail acceptance. Exploration mode is report-only.
Documentation changes do not require repeating quantification measurements.
