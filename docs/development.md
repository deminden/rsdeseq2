# Development

## Repository Structure

- `crates/rsdeseq2`: Rust numerical core and minimal CLI.
- `r-pkg/rsdeseq2`: experimental R package scaffold and R CI surface.
- `scripts`: reference generation and future benchmark scripts.
- `docs`: algorithms, compatibility, reproducibility, and release notes.
- `results/parity`: generated DESeq2 reference outputs.
- `results/benchmarks`: benchmark outputs.

The current implementation is Rust-first: a Rust core, minimal Rust CLI,
scripts for parity fixtures, docs, CI, validation outputs, and an experimental
R package scaffold. R/DESeq2 is used only as an external reference generator
for tests.

## Commands

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
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
scope. The default generated set is intended to be green and includes
DESeq2-backed beta-prior variance, refit, and estimated-prior refit anchors,
weighted fixed-dispersion Wald/LRT, unweighted GLM-mu mean-trend MAP/Wald/LRT anchors
with and without Cox-Reid, weighted GLM-mu Cox-Reid gene-wise anchors, and the
current weighted GLM-mu mean and local-trend MAP/Wald/LRT anchors, including
result-row BH-adjusted p-value and compact result-table checks for those
matched Wald/LRT branches, plus the current unweighted GLM-mu local-trend
MAP/Wald/LRT result-table fixture and the unweighted/weighted GLM-mu Cox-Reid
local-trend MAP/Wald/LRT anchors.

Run benchmarks:

```bash
cargo bench -p rsdeseq2
scripts/benchmark_rsdeseq2.sh --genes 1000 --samples 8 --repeats 1
```

## Coding Conventions

- Keep statistical computation in Rust; see the implementation plan for future
  wrapper constraints.
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
