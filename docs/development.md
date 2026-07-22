# Development

## Repository Structure

- `crates/rsdeseq2`: Rust numerical core and minimal CLI.
- `r-pkg/rsdeseq2`: R access layer and package CI.
- `scripts`: reference-generation, parity-analysis, and benchmark scripts.
- `docs`: algorithms, compatibility, reproducibility, and release notes.
- `results/`: ignored generated parity, fixture, and benchmark outputs.

The repository is Rust-first: a Rust core, CLI, R access layer,
scripts for parity fixtures, docs, CI, and validation outputs. R/DESeq2 is used
as an external reference generator for tests.

The minimum supported toolchain is Rust 1.97.1, declared by `rust-version` in
`Cargo.toml`; the workspace uses the Rust 2024 edition. Development and CI use
that toolchain or newer. R package CI uses R 4.6.1, while numerical reference
fixtures retain their recorded R versions.
Public slice-like APIs prefer `RangeBounds` where range inputs are accepted, so
callers can use legacy range syntax and the `core::range` types.

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

The optimizer stress reference specifically requires R 4.6.1:

```bash
OPENBLAS_NUM_THREADS=1 \
  Rscript scripts/generate_lbfgsb_synthetic_stress_fixtures.R
```

The `rcompat-lbfgsb` 0.2.1 result is 512/512 exact endpoints, objectives, and
evaluation counts. The 0.1.6 baseline matched endpoint plus objective in
493/512 cases and the objective alone in 507/512 cases at practical tolerances,
0/512 exactly, and 311/512 for exact evaluation counts.

The reference tests skip when
`crates/rsdeseq2/tests/data/deseq2_reference/` is absent. When present, the
tests use these files to compare implemented Rust stages with generated DESeq2 outputs. Wald/LRT
golden checks use supplied dispersions and `DESeq2:::fitNbinomGLMs` to match the
fixed-dispersion scope. The generated reference set includes
DESeq2-backed beta-prior variance, refit, and estimated-prior refit references,
weighted fixed-dispersion Wald/LRT, unweighted GLM-mu mean-trend MAP/Wald/LRT references
with and without Cox-Reid, weighted GLM-mu Cox-Reid gene-wise references, and the
weighted GLM-mu mean and local-trend MAP/Wald/LRT references, including
result-row BH-adjusted p-value and compact result-table checks for those
matched Wald/LRT branches, plus the unweighted GLM-mu local-trend
MAP/Wald/LRT result-table fixture and the unweighted/weighted GLM-mu Cox-Reid
local-trend MAP/Wald/LRT references.

Run benchmarks:

```bash
cargo bench -p rsdeseq2
scripts/benchmark_rsdeseq2.sh --genes 1000 --samples 8 --repeats 1
```

After producing the top-1,000 benchmark diagnostics described in
[reproducibility.md](reproducibility.md#versioned-high-error-benchmark),
reproduce the reported measurements and run the scorer tests:

```bash
python3 scripts/score_frozen_worst_genes.py \
  --fixture docs/data/wald_frozen_worst100_r461.tsv \
  --diagnostics results/benchmarks/frozen_worst100_diagnostics.tsv \
  --report-only
python3 -m unittest discover -s scripts/tests -p 'test_*.py'
```

## Coding Conventions

- Keep statistical computation in Rust; see the implementation plan for
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
