# Release

No production release should claim DESeq2 parity until stage-by-stage reference
tests cover the implemented statistical pipeline.

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
below `2.4e-7` for the tested result fields. The v0.2.5 benchmark
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
