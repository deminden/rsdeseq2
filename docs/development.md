# Development

## Repository Structure

- `crates/rsdeseq2`: Rust numerical core and minimal CLI.
- `r-pkg/rsdeseq2`: R package scaffold.
- `scripts`: reference generation and future benchmark scripts.
- `docs`: algorithms, compatibility, reproducibility, and release notes.
- `results/parity`: generated DESeq2 reference outputs.
- `results/benchmarks`: benchmark outputs.

The structure follows the broad organization of `rsfgsea`: Rust-first core,
language wrappers, scripts, docs, CI, and validation outputs.

## Commands

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```

Run the R wrapper tests from the source tree:

```bash
Rscript -e 'for (f in list.files("r-pkg/rsdeseq2/R", pattern="[.]R$", full.names=TRUE)) source(f); testthat::test_dir("r-pkg/rsdeseq2/tests/testthat", reporter="summary")'
```

Run the R package check, including the registered C `.Call` bridge:

```bash
R CMD build r-pkg/rsdeseq2
R CMD check --no-manual --no-build-vignettes rsdeseq2_*.tar.gz
```

Generate R references:

```bash
Rscript scripts/generate_deseq2_references.R
cargo test -p rsdeseq2 --test dispersion_reference
cargo test -p rsdeseq2 --test results_reference
cargo test -p rsdeseq2 --test wald_reference
cargo test -p rsdeseq2 --test lrt_reference
```

The reference tests are skip-safe when
`crates/rsdeseq2/tests/data/deseq2_reference/` has not been generated. Once the
R script is run, they compare the implemented Rust stages against the generated
DESeq2 fixture files. Full DESeq2 result references are written for future
dispersion parity, while the current Wald/LRT golden checks use supplied
dispersions and `DESeq2:::fitNbinomGLMs` to match the Rust fixed-dispersion
scope. The script's default output is intended to be green; use
`Rscript scripts/generate_deseq2_references.R --include-known-gaps` only when
you intentionally want exploratory fixtures for currently divergent internal
weighted paths.

Run benchmarks:

```bash
cargo bench -p rsdeseq2
```

## Coding Conventions

- Keep the numerical core independent from R and Python bindings.
- Prefer explicit structs and enums over string options.
- Return `DeseqError::UnsupportedFeature` for unimplemented stages.
- Avoid panics in library code.
- Add hand-computable tests before reference tests.
- Keep compatibility behavior documented when inferred from DESeq2 docs or
  source.

## Adding a New Statistical Stage

Each stage should expose intermediate output in `DeseqFit` or a stage-specific
struct. When a field maps directly to DESeq2 row metadata, also update the
diagnostic alias view in `diagnostics.rs`. Add tests in three layers where
practical:

1. Hand-computable toy test.
2. DESeq2 golden-reference comparison.
3. Stress/property test for edge cases.
