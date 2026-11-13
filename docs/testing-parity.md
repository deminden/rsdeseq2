# Testing and numerical parity

Version 0.5.0 is checked against independent original-R references and actual
standalone/wrapped execution. Accuracy and interface equality are separate gates:
matching the same Rust engine twice does not prove agreement with DESeq2.

## CI status

All nine executed jobs passed on source commit
[`2d86916`](https://github.com/deminden/rsdeseq2/commit/2d86916c089f55bcddb33aee1af9b079f2b9384e)
on 2026-11-12 UTC. Final release preparation updates documentation and package
discovery metadata; numerical code, tests and retained reference values are
unchanged. The tag-triggered release workflow validates the final source again.

| Workflow | Passing checks |
| --- | --- |
| [Rust](https://github.com/deminden/rsdeseq2/actions/runs/43331557396) | Vendor/version integrity, formatting, Clippy, full workspace and Python tooling tests |
| [R](https://github.com/deminden/rsdeseq2/actions/runs/43331557387) | Both pinned R 4.6.2 and current R: package checks and exact standalone/wrapped comparisons. Pinned R also checks fresh formula/SQP/uniform references, interruption/reuse, sanitizers and actual R payload gates |
| [Five platforms](https://github.com/deminden/rsdeseq2/actions/runs/43331557382) | Linux x86_64/aarch64, macOS x86_64/aarch64 and Windows x86_64 numerical/CLI panels in release profile |
| [Full GTEx](https://github.com/deminden/rsdeseq2/actions/runs/43331570676) | All 56,200 genes, 61 samples and six Wald result fields |

[Machine-readable evidence](data/ci_0.5.0.json) binds run/job IDs and downloaded
log hashes to that commit. The full-only Parity dispatch intentionally skips its
ordinary-oracle job because Rust CI already runs those tests. No checks were
weakened and no reference bounds were widened.

Both R package checks retain one known compiled-symbol warning for Rust runtime
termination entry points. Their tests pass; this is not a warning-free package
or CRAN-acceptance claim.

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

## Compact and full GTEx checks

Ordinary tests fit a fixed-dispersion panel of **512 genes × all 61 samples**,
selected from count/reference tails, quantiles, missingness and deterministic
background coverage. It checks all six result fields and both directions of
missingness. Full-cohort size factors and final dispersions are supplied; the
subset's BH-adjusted p-values are recomputed for that subset. A separate default
test estimates normalization over all genes.

The one opt-in full replay remains necessary for all-gene coverage. The compact
panel does not replace native dispersion, weighted, local-fit or optimizer tests.
Full counts are also retained for offline contrast generation. See the
[fixture description](../crates/rsdeseq2/tests/data/gtex_v8_real_wald/README.md).

```sh
cargo test --locked --workspace
python3 -m unittest discover -s scripts/tests -p 'test_*.py'
cargo test --locked -p rsdeseq2 --test gtex_real_wald_parity full_gtex_wald \
  -- --ignored --nocapture
```

The eleven full biological/random contrasts provide a separate acceptance panel;
see [reproducibility](reproducibility.md#eleven-contrast-acceptance).

## R package and interface checks

Use the pinned environment and its restored dependencies. From the repository root:

```sh
repo_root="$PWD"
mkdir -p results/checks/r
(cd results/checks/r && \
  R CMD build "$repo_root/r-pkg/rsdeseq2" && \
  R CMD check --no-manual --no-build-vignettes rsdeseq2_0.5.0.tar.gz)
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

Use fresh output directories. `.github/workflows/r.yml` lists the complete
matrix/fresh-fit, expression, factor, frozen-local, mixture, posterior, warning
and live-reference identity panels; it is the executable combined gate.

## Platform arithmetic and native safety

The pinned Linux lane independently regenerates seven formula panels under
R 4.6.2/OpenBLAS Haswell and compares every non-provenance field. Fresh SQP and
uniform-posterior outputs remain separate from committed fixtures. Windows and
Intel macOS generate independent uniform CDFs at the existing bounds; other
posterior fields still use their original fixtures. See [platform arithmetic](linux_numerical_parity.md).

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
