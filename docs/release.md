# Release

No production release should claim DESeq2 parity until stage-by-stage reference
tests cover the implemented statistical pipeline.

## First WIP Commit Checklist

The first public commit may be published as an experimental foundation if it is
clearly labeled work in progress and does not claim DESeq2 compatibility beyond
the tested primitive stages.

Before the first commit:

- Keep the local DESeq2 inspection clone under ignored `external/`.
- Keep Rust `target/`, R check directories, package tarballs, native objects,
  and generated parity/benchmark outputs out of git.
- Keep the README warning that the package is not production-ready.
- Run Rust formatting, linting, and tests.
- Run the R source-tree tests and `R CMD check` for the scaffold package.
- Review staged files for accidental large artifacts or generated references.

Before any release:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
Rscript -e 'for (f in list.files("r-pkg/rsdeseq2/R", pattern="[.]R$", full.names=TRUE)) source(f); testthat::test_dir("r-pkg/rsdeseq2/tests/testthat", reporter="summary")'
R CMD build r-pkg/rsdeseq2
R CMD check --no-manual --no-build-vignettes rsdeseq2_*.tar.gz
Rscript scripts/generate_deseq2_references.R
```

The DESeq2 reference-generation script requires an R environment with
Bioconductor DESeq2 installed. Generated references should be reviewed before
they are committed.
