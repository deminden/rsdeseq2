# Release

No production release should claim DESeq2 parity until stage-by-stage reference
tests cover the implemented statistical pipeline.

## Version 0.3.0

This minor release combines the validated numerical correction in `9acb44b`
with the completed R wrapper integration. The Rust crate, R package, R backend,
vendored core, and lockfiles are aligned at 0.3.0.

### Automatic dispersion means

The CLI, top-level Rust workflows, and R backend now select DESeq2's automatic
dispersion mean fit: linear means for unweighted designs with as many distinct
sample rows as coefficients, and GLM means otherwise. Supplied observation
weights select GLM means, including uniform weights. The explicitly forced
GLM estimator retains its behavior.

The frozen comparison covers 11 GTEx v8 Kidney Cortex contrasts: two biological
comparisons, the original blocked contrast, five stratified random splits, and
three balanced random splits. Each contrast uses all 56,200 genes. All 121
per-contrast checks and all nine global maximum-error checks pass. The smallest
global improvement is 596.7-fold; the smallest targeted per-contrast
improvement is 664-fold. Gene coverage and finite/NA patterns are preserved,
and blocked controls and base means retain their previous error limits.

This evidence uses R 4.6.1 / DESeq2 1.52.0, parametric dispersion trends, no
beta prior, and disabled Cook's filtering/replacement and independent
filtering. The complete measurements, remaining error bounds,
input/reference hashes, and reproduction commands are in the
[GTEx dispersion correction report](gtex-dispersion-improvement.md).
The core numerical behavior matches `9acb44b`; subsequent precision
experiments were reverted. Earlier release benchmarks below remain historical.

### R wrapper integration

Fitted `DESeqDataSet` objects receive Rust fitted means and Cook's-distance
assays, fitted row metadata, and per-coefficient MLE and standard-error
columns. Cook's outlier replacement and refitting run for eligible replicate
samples, with a default minimum of seven replicates per design cell.
`minReplicatesForReplace` accepts an integer of at least three or `Inf` to
disable replacement. Original counts and Cook's distances are preserved;
replacement assays and flags are attached when rows are replaced.

Dataset result tables carry Bioconductor column descriptions and result
metadata, including independent-filtering and replacement/refit information.
Matrix and list fits retain their documented result interface. Seeded
integration tests cover fitted means, Cook's distances, replacement/refitting,
and result metadata. See [R wrapper status](r-wrapper.md) for the supported
scope and remaining limitations.

### Current parity, speed and memory evidence

The [11-contrast native Wald benchmark](gtex-native-wald-benchmark.md) contains
current parity results, matched Rust/DESeq2 process timings and peak resident
memory, and reproducible contrast-generation and benchmark commands. Detailed
evidence lives in the documentation; the main README links to the current results.

### Local release validation

Validated on 2026-10-14 with R 4.6.1 / DESeq2 1.52.0:

- Rust workspace: 1,064 tests passed; formatting and Clippy passed.
- R Rust bridge: formatting and Clippy passed.
- R source package: `R CMD check --no-manual --no-build-vignettes` reported
  `Status: OK`, with 860 passing test assertions and no skipped tests.
  `apeglm` 1.34.0 was installed in the pinned `deseq2-r461` R environment;
  DESeq2 remained at 1.52.0. The final check required all suggested packages.
- Python script tests: 18 passed.
- Full GTEx replay: 121/121 per-contrast and 9/9 global checks passed; all
  22 result/diagnostic files matched `9acb44b` exactly.
- Committed contrast generator: all 22 count/design files reproduced exactly.
- Matched native Wald timing: all 11 contrasts completed; time ratios were
  1.24–1.52× (DESeq2 time divided by Rust time), one process run per tool and
  contrast. Timed Rust results matched the validated outputs exactly.
- Matched native Wald peak RSS: all 11 contrasts completed; memory ratios were
  5.33–6.87× (DESeq2 peak RSS divided by Rust peak RSS), or
  81.2–85.4% lower peak RSS for Rust. All Rust memory-run results matched
  `9acb44b` exactly; frozen inputs and references were unchanged.
- Rust crate packaging and verification passed with the locked dependencies.

The current release comparison uses the complete 11-contrast GTEx suite.
Historical 100-row benchmark measurements below retain their original release
provenance; no new measurement of that older study is claimed here.

## Version 0.2.6

This release adds dataset-based R `vst()`,
`varianceStabilizingTransformation()`, and `rlog()` workflows backed by the
Rust transformations, plus `lfcShrink()` support for normal, `ashr`, and
`apeglm` methods on the documented inputs. Frozen trend reuse, weighted fast
VST subset behavior, and transformed-factor coefficient aliases are covered by
integration tests. The Rust dependency refresh includes `rcompat-locfit`
0.1.5, `rcompat-lbfgsb` 0.2.2, `clap` 4.6.7, `serde` 1.0.229,
`serde_json` 1.0.154, and `thiserror` 2.0.21. The Rust crate, R package,
vendored core, and their lockfiles are aligned at 0.2.6. A committed GTEx v8
fixture checks fixed-dispersion Wald parity across all 56,200 genes for 61
Kidney Cortex samples; maximum absolute differences from DESeq2 1.52.0 are
below `3.8e-7` for the tested result fields. The v0.2.5 benchmark
results below remain historical; no 0.2.6 release-to-release benchmark is
claimed here.

## Pre-Release Checklist

Keep the Rust crate, R package, and R-package Rust stub on the same release
version. Benchmark comparisons may retain an older version as an explicit
baseline.

Before preparing a release commit:

- Keep DESeq2 source-inspection copies under ignored `external/`.
- Keep Rust `target/`, generated archives/native objects, and generated
  parity/benchmark outputs out of git.
- Keep the README scoped to validated workflows and avoid claiming full DESeq2
  workflow parity.
- Run Rust formatting, linting, tests, package checks, and R wrapper checks.
- Run the versioned 100-row parity benchmark and record its measurements. The
  v0.2.5 median, mean, and maximum absolute errors are
  `6.063612945084174e-10`, `1.260818774570247e-4`, and
  `1.5259081158007781e-3`. The frozen v0.2.4 baseline measured
  `1.4637657972313423e-4`, `3.793917566690452e-4`, and
  `3.0938714191082184e-3`, respectively. The v0.2.5 median and mean are
  241401.589x and 3.00909032x lower; the maximum is 50.6796531% lower, with a
  v0.2.5-to-v0.2.4 ratio of `0.493203469`. In the fixed set, 89/100 rows
  improved and 78/100 improved by at least 10x. The run includes compensated
  accumulation of log counts when computing per-gene geometric means for
  ratio size factors.
- Review staged files for accidental large artifacts or generated reference
  outputs.

For 0.3.0, replay the frozen 11-contrast GTEx suite with the release binary:

```bash
cargo build --release -p rsdeseq2 --locked
python3 scripts/check_gtex_dispersion_improvement.py \
  --output-dir results/benchmarks/release_0.3.0/gtex
```

Use the existing reference files and their recorded hashes. The release
checker must pass all 121 per-contrast and nine global checks. Compare its
numeric results and fit diagnostics with the saved `9acb44b` outputs as well;
the version bump must preserve those values exactly. New reference
generation is a fixture maintenance step; release validation replays the
frozen evidence.

Before any release:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo package -p rsdeseq2 --locked
scripts/benchmark_rsdeseq2.sh --genes 1000 --samples 8 --repeats 1
Rscript scripts/generate_deseq2_references.R
python3 -m unittest discover -s scripts/tests -p 'test_*.py'
python3 scripts/score_frozen_worst_genes.py \
  --fixture docs/data/wald_frozen_worst100_r461.tsv \
  --diagnostics results/benchmarks/frozen_worst100_diagnostics.tsv \
  --report-only
```

The DESeq2 reference-generation script requires an R environment with
Bioconductor DESeq2 installed. Generated references should be reviewed before
they are committed. The frozen scorer additionally requires the saved
frozen-benchmark real-data diagnostics generated with a release binary and
`--diagnostics-limit 69045`; see
[reproducibility.md](reproducibility.md#versioned-high-error-benchmark).

## Crates.io Release

Crates.io publishing is handled by `.github/workflows/publish-crates.yml`.
Configure the repository secret `CARGO_REGISTRY_TOKEN` with a crates.io API
token that can publish `rsdeseq2`.

GitHub Release binaries are handled by
`.github/workflows/release-binaries.yml`. The workflow builds `rsdeseq2` for
Linux x86_64, Linux aarch64, macOS x86_64, macOS aarch64, and Windows x86_64,
packages each binary with the README and license, and uploads the archives to
the matching GitHub release.

To release version `X.Y.Z`:

1. Set `crates/rsdeseq2/Cargo.toml` to `version = "X.Y.Z"`.
2. Run the release checklist above locally.
3. Commit the release changes.
4. Push tag `vX.Y.Z`.

The GitHub Actions publish job runs formatting, clippy, the full workspace test
suite, and `cargo package -p rsdeseq2 --locked` before calling
`cargo publish`. It refuses to publish if the pushed tag version does not match
the crate version.
The binary release job performs the same tag/version check before building and
attaching release assets.
