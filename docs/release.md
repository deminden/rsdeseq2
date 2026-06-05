# Release

No production release should claim DESeq2 parity until stage-by-stage reference
tests cover the implemented statistical pipeline.

## Pre-Release Checklist

Before preparing a release commit:

- Keep local DESeq2 inspection clones under ignored `external/`.
- Keep Rust `target/`, generated archives/native objects, and generated
  parity/benchmark outputs out of git.
- Keep the README scoped to validated workflows and avoid claiming full DESeq2
  workflow parity.
- Run Rust formatting, linting, tests, package checks, and R wrapper checks.
- Review staged files for accidental large artifacts or generated local
  references.

Before any release:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo package -p rsdeseq2 --locked
scripts/benchmark_rsdeseq2.sh --genes 1000 --samples 8 --repeats 1
Rscript scripts/generate_deseq2_references.R
```

The DESeq2 reference-generation script requires an R environment with
Bioconductor DESeq2 installed. Generated references should be reviewed before
they are committed.

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
