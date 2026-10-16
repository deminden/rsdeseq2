# Large-file refactor

The 2026-10-16 refactor starts from `v0.3.1`
(`a683dda9dc26d1e085e87f4beed87bb74ae35059`). It reorganizes the three largest
integration-test suites, the R results wrapper, formula handling, IRLS, beta
fitting, and the R Rust bridge. Numerical operation bodies, public Rust paths,
R function definitions and C export signatures are preserved.

## Organization

| Area | Before | Entry file after | Topic files |
| --- | ---: | ---: | ---: |
| Native Wald integration tests | 6,131 lines | 274 lines | 12 |
| Design integration tests | 5,843 lines | 38 lines | 11 |
| Results integration tests | 4,259 lines | 71 lines | 11 |
| R results wrapper | 2,389 lines | 259 lines | 9 including the entry file |
| Formula handling | 2,308 lines | 139 lines | 10 including the entry file |
| IRLS | 2,374 lines | 410 lines | 9 plus the entry file |
| Beta fitting | 1,879 lines | 170 lines | 8 plus the entry file |
| R Rust bridge | 1,555 lines | 193 lines | 7 plus the entry file |

The [complete inventory](data/refactor_large_files_2026-10-16_inventory.tsv)
records entry and largest-topic file lengths, with topic-file counts excluding
the entry file. New test topic files are at most
907 lines; new production topic files are at most 497 lines, apart from the
516-line private IRLS regression-test module.

The existing three integration-test executables are retained. Their 216 cases
are grouped into child modules, with shared fixtures in the original entry
files. Test names gain a topic-module prefix; their original function names
and assertion bodies are unchanged.

Formula code follows the existing private implementation layout: entry points,
model frames, types, offsets, transform parsing, scaling, polynomial arguments
and arithmetic each have a focused file. IRLS and beta code similarly separate
fitting, optimizer fallback, arithmetic, linear algebra and validation. These
Rust parts use `include!` in their original module scope and item order, so
public paths and private helper access remain unchanged. Function bodies and
inline attributes are retained.

The R wrapper separates result dispatch, fits, object adapters, filtering,
contrast parsing, validation, statistics, all-zero handling and contrast-name
resolution. All files contain function definitions, with no added load-time
calls or dependence on file collation order. The Rust bridge keeps common
buffer/error handling in its entry file and groups existing C exports by
operation. The R package's vendored Rust core is regenerated from the main
crate with `scripts/sync_r_vendor.sh`.

## Validation

- Rust source tokens match the original IRLS, beta, formula and bridge code
  after expanding the new include files. Reconstructing the three test suites
  also produces identical tokens, including their assertions and fixtures.
- All 87 R results-wrapper definitions have identical parsed syntax, including
  formal arguments and bodies, with source-location attributes excluded.
- Rust workspace: **1,064 tests passed**, with no failures or ignored tests.
  The R Rust bridge's diagnostic-schema unit test also passed.
- R 4.6.1 package check: **Status: OK**, **860 assertions passed**, with no
  failures, warnings or skipped tests. The existing pinned local test library
  supplies `testthat`; DESeq2 remains at 1.52.0.
- Core and R bridge formatting and all-target Clippy passed.
- Rust crate packaging and verification compilation passed, including the
  split source and test files. All 126 core source files match the R vendor.
- The eleven full GTEx contrasts pass **121 per-contrast and nine global
  precision gates**. All **22 result and diagnostic exports** match 0.3.1
  byte for byte; frozen references and inputs are unchanged.

## Speed

One sequential pass across eleven contrasts, all 56,200 genes per contrast,
on Apple M5/macOS with Rust 1.99.0 and one fitting and BLAS thread measured **39.644 s** for
the saved local 0.3.1 release build and **36.687 s** for the refactored build.
The summed-time ratio is **1.0806×** (baseline / refactored), with no
aggregate slowdown. The median per-contrast ratio is **1.0079×**, a
roughly **0.8%** difference in the typical pair.

| Contrast | 0.3.1 baseline (s) | Refactored (s) | Baseline / refactored |
| --- | ---: | ---: | ---: |
| `condition_blocked` | 4.152 | 3.861 | 1.0755× |
| `sex_female_vs_male` | 2.711 | 2.661 | 1.0190× |
| `male_age4_vs_age3` | 2.272 | 2.255 | 1.0079× |
| `blocked_random_rep01` | 3.942 | 3.930 | 1.0031× |
| `blocked_random_rep02` | 3.989 | 3.989 | 1.0001× |
| `blocked_random_rep03` | 4.003 | 3.983 | 1.0049× |
| `blocked_random_rep04` | 4.300 | 4.355 | 0.9874× |
| `blocked_random_rep05` | 3.944 | 3.968 | 0.9941× |
| `unblocked_random_rep01` | 2.632 | 2.509 | 1.0491× |
| `unblocked_random_rep02` | 4.966 | 2.510 | 1.9788× |
| `unblocked_random_rep03` | 2.730 | 2.667 | 1.0238× |

An isolated slow baseline run for `unblocked_random_rep02` affects the summed
times. These single observations do not estimate repeated-run uncertainty;
the comparison supports preserved speed for this workflow. The refactor keeps
execution logic and inline attributes intact. Tool order rotates across
contrasts; all fitting processes run sequentially, with compilation and other
fits stopped. Timings include startup, input loading, fitting and result-table
output. Separate diagnostic validation is outside the timed processes.

## Recorded evidence

- [File inventory](data/refactor_large_files_2026-10-16_inventory.tsv)
- [Validation counts and binary hashes](data/refactor_large_files_2026-10-16_validation.json)
- [Per-contrast time and memory](data/refactor_large_files_2026-10-16_speed_memory.tsv)
- [All 33 process measurements](data/refactor_large_files_2026-10-16_runs.tsv)
- [Hardware and thread provenance](data/refactor_large_files_2026-10-16_provenance.json)
- [Frozen input/reference hashes](data/refactor_large_files_2026-10-16_input_sha256.tsv)
- [All 22 result/diagnostic hash comparisons](data/refactor_large_files_2026-10-16_gtex_sha256.tsv)

Both binaries embed the same Rust compiler commit
`b940084d7eb6a299eb4bfeb8e34901bc051e7ac4`.
The timing baseline is the saved local release-preparation binary with SHA-256
`f759f546185f4c80b7b49001983fe8cc53c41c0cb595f9fcac81a58bc2b06b1f`. The candidate binary SHA-256 is
`72a85c9c7ddbcd19b5c2afdc9dfc2ebd88a7ae588a1e5112d59cbce964ab5043`. The same committed benchmark script runs
both binaries and DESeq2; the table above compares the two native builds.

## Reproduce

Save a release binary built from `v0.3.1` before applying the refactor, then
build the candidate with the same toolchain and locked dependencies. Existing
contrast generation and frozen-reference setup are documented in the
[GTEx benchmark](gtex-native-wald-benchmark.md#reproduce).

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace
cargo fmt --manifest-path r-pkg/rsdeseq2/src/rust/Cargo.toml -- --check
cargo clippy --manifest-path r-pkg/rsdeseq2/src/rust/Cargo.toml \
  --all-targets --locked -- -D warnings
cargo build --release -p rsdeseq2 --locked
python3 scripts/check_gtex_dispersion_improvement.py \
  --workers 1 --output-dir results/refactor_large_files/gtex
python3 scripts/benchmark_gtex_native_wald.py \
  --binary target/release/rsdeseq2 \
  --baseline-binary /path/to/saved/rsdeseq2-0.3.1 \
  --baseline-label 0.3.1 \
  --validated-results-root results/refactor_large_files/gtex \
  --rscript /path/to/deseq2-r461/bin/Rscript \
  --repeats 1 --rust-threads 1 --blas-threads 1 \
  --output-dir results/refactor_large_files/speed_memory
```

Build the R source package and run `R CMD check --no-manual
--no-build-vignettes` in the pinned R environment with all suggested packages
available. Run builds, tests and numerical replay before timing; keep other
fitting and compilation stopped for the complete timing pass.
