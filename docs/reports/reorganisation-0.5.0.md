# 0.5.0 source layout and local preparation

This preparation reorganises the validated `28057cc` source and aligns own
package versions at 0.5.0. Publication and tags remain on hold. The checks below
were completed locally on 2026-11-12, using native macOS ARM, R 4.6.2 and the
pinned `deseq2-r462-native` comparison environment.

## Responsibilities

- `crates/rsdeseq2/` remains the canonical Rust core; `r-pkg/rsdeseq2/` remains
  the self-contained R access package and its generated vendor.
- `scripts/` groups references, validation, benchmarks, environment installation,
  release tooling and tests. [Tooling map](../../scripts/README.md).
- `docs/` contains current guides, `reports/` current empirical validation,
  `history/` superseded plans/checkpoints, `releases/` versioned notes, and
  `data/` immutable supporting evidence. [Documentation map](../README.md).
- `config/r-reference/` contains operational pinned R environment/source locks.
- Internal Rust validation drivers live under `tools/validation/` in the crate.
  Explicit Cargo example targets preserve tool names and executable paths while
  keeping them out of the installed CLI set.
- The root quickstart inputs and redundant CI smoke step are removed. Existing
  numerical CLI integration coverage remains; docs use caller-supplied tables.
- Old local R-check output, an old source tarball and the optimisation journal
  are preserved under ignored `results/local/`. Build caches, benchmark datasets
  and saved validation evidence are retained.

## Invariants

Production Rust source and independent frozen expected values retain their
original bytes. Calculations, tolerances and operation order are unchanged.
Metadata versions and tooling paths change deliberately. Historical JSON/TSV
provenance keeps its original paths, hashes and source identities.

`.DS_Store` is ignored and explicitly excluded from Cargo and R source packages,
including nested paths. The package inventory check includes existing local OS
metadata and verifies that it cannot enter distribution archives.

## Local acceptance

| Check | Result |
| --- | --- |
| Full Rust workspace | 1,158 tests passed; one existing opt-in test ignored |
| Full GTEx opt-in replay | Passed: 56,200 genes × 61 samples; 48,683 finite Wald rows |
| R source-package check | 2,870 assertions passed, no test warnings/skips/failures; one known package warning |
| Python tooling, including actual R payload comparisons | 49 tests passed; no skips |
| Standalone/R identity checks | All nine selected panels passed with exact floating-point comparisons |
| Formatting and all-target Clippy | Passed; Clippy treats warnings as errors |
| Minimum Rust 1.97.1 | Locked workspace/all-target compilation passed |
| R environment restoration check | All 82 pinned package versions verified |
| Vendor, notices and own package versions | Match; Rust/R/backend/vendor and lockfiles agree at 0.5.0 |
| Rust distribution | 0.5.0 crate built and its unpacked package compiled successfully |
| R distribution | 0.5.0 source archive built; package check completed |
| Package inventories | Neither archive contains `.DS_Store` or build caches |
| Paths and documentation | No stale moved paths or broken local Markdown targets; whitespace checks passed |

The standalone/R panels cover primitives, matrix fresh fits, formula fits,
formula expressions, native formulas/factors, frozen-local curves and normal and
uniform posterior summaries. They were invoked by absolute script path from
outside the checkout, exercising relocated helper and fixture discovery.

The full GTEx replay retained the existing bounds. Maximum absolute differences
against frozen original-DESeq2 values were `5.588e-9` for baseMean, `2.386e-7`
for log2 fold change, `4.140e-10` for SE, `3.734e-7` for the statistic,
`4.955e-8` for p-value and `8.221e-7` for adjusted p-value. Maximum relative
size-factor error was `5.551e-15`.

R's source-package check retains its established compiled-symbol warning for
Rust runtime `abort`/`exit` entry points. Its test suite passed; the package check
is not described as warning-free or as CRAN acceptance.

## Precision and speed scope

SHA-256 comparison with `28057ccbafae23b1455a7442d1eeff9200b958cd`
confirms byte-identical contents for all 286 canonical production Rust files,
39 R wrapper files and 46 native bridge source files. All 556 frozen reference,
measurement and operational-lock files retain their original hashes, including
the two locks moved into `config/r-reference/`. Reference README links are
updated separately. The Rust integration-test changes only redirect shared
validation-module paths; no numerical assertions or tolerances change.

This preparation changes layout, tooling imports/path discovery, package
metadata and documentation. It does not change statistical arithmetic or
numerical execution order. The existing artifact-bound speed and memory
comparisons remain in the [current benchmark](gtex-native-wald-benchmark.md) and
[correction report](review-fixes-2026-11-12.md). No new performance estimate is
claimed for 0.5.0, and historical measurements are not relabelled. The local
checks above do not substitute for a future Linux/multi-target CI run on the
new committed layout.

Ignored detailed outputs are retained under `results/reorganisation-0.5.0/`.
The checked local source archives are `target/package/rsdeseq2-0.5.0.crate` and
`results/reorganisation-0.5.0/package/rsdeseq2_0.5.0.tar.gz`. No tag, push,
release upload or registry publication was performed during this preparation.
