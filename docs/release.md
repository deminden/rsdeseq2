# Release

No production release should claim DESeq2 parity until stage-by-stage reference
tests cover the implemented statistical pipeline.

## First WIP Commit Checklist

The first public commit may be published as an experimental foundation if it is
clearly labeled work in progress and does not claim DESeq2 compatibility beyond
the tested primitive stages.

Before the first commit:

- Keep the local DESeq2 inspection clone under ignored `external/`.
- Keep Rust `target/`, generated archives/native objects,
  and generated parity/benchmark outputs out of git.
- Keep the README warning that the package is not production-ready.
- Run Rust formatting, linting, and tests.
- Review staged files for accidental large artifacts or generated references.

Before any release:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
scripts/benchmark_rsdeseq2.sh --genes 1000 --samples 8 --repeats 1
Rscript scripts/generate_deseq2_references.R
```

The DESeq2 reference-generation script requires an R environment with
Bioconductor DESeq2 installed. Generated references should be reviewed before
they are committed.
