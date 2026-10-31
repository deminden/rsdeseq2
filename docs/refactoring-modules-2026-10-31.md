# Module refactoring — 2026-10-31

This change reorganizes the candidate after `9d4482d`. It does not complete the
[Rust-only wrapper migration](rust-only-wrapper-migration.md) or constitute a
new published release.

## Changes

- VST fast-subset preparation and local integration/FMM spline calculations
  have separate private modules, with their existing public exports retained.
- Lowess calculation is separated from independent-filter selection. Cook's
  replacement/refit planning and robust trimmed variance are private modules.
- The CLI's Wald and LRT branches are separate handlers under the existing
  staging/publish transaction.
- R's normal, ashr and apeglm shrinkage handlers and expanded-design preparation
  are separate source files; their function bodies are byte-identical to the
  previous combined file.
- The large IRLS and contrast integration-test files are organized by behavior.
  All 46 IRLS and 51 contrast test functions remain present.

The moved VST, lowess and Cook's kernel source matches the prior source apart
from module visibility and formatting; arithmetic expressions and their order are
unchanged. Public Rust paths and R exports are retained. The R package vendors
exactly the same reorganized core source.

## Numerical and package checks

The full Rust suite passes 1,072 tests, the backend passes three, and the
checked R source package passes 2,615 assertions with no failures, warnings or
skips. Root and backend all-target Clippy checks with warnings denied, Rust
1.97.1 checks, vendor integrity and the Python suite (46 tests, one platform
skip) pass. The two reorganized integration suites retain and pass all 97 test
functions.

All 121 per-contrast and nine global frozen GTEx precision gates pass; all 22
result/diagnostic exports are byte-identical to the preceding 0.4.0 candidate.
The independent Pasilla Wald/LRT comparison passes its bounds and all eight
result, diagnostic, coefficient and SE exports are exact. Direct standalone/R
checks pass 189 primitive, 188 expression, 54 native-formula, 51 factor-frame
and nine expanded-factor cases. Formula APIs with different public label order
are compared after their explicit validated column mapping.

The crate archive includes 499 verified source, test and example files; the R
source archive contains 333 verified source members, excluding its
build-normalized DESCRIPTION. The tested native CLI has SHA-256
`ce1194830a7eb4f1055beeda4053a12335b33f78de3a07c2dc0e0d3aca2a8e18`.
Earlier sanitizer checks apply to their recorded pre-refactor artifacts; no fresh
sanitizer run is claimed for this source-only rearrangement.

## Performance status

The paired native timing run was stopped after extreme external contention made
its result unsuitable for checking unchanged speed. Nine of eleven contrasts
completed, with 28 timed processes retained; the R workflow timing panel did
not start. The observed machine load reached 60 or more, with a separate
virtual machine and Rust compilers consuming CPU. In the fourth blocked split,
the previous binary took 21.10 s and the candidate 7.69 s; in the fifth, those
figures were 10.25 s and 26.05 s. These are not evidence of an isolated speed
change. No unchanged case was retimed and no optimization round is claimed.
**Speed and memory parity remain unverified by this measurement.**

Evidence is retained in the [validation record](data/refactor_modules_2026-10-31_validation.json),
[GTEx precision gates](data/refactor_modules_2026-10-31_gtex_precision.tsv),
[global gates](data/refactor_modules_2026-10-31_gtex_global_precision.tsv),
[22 export identities](data/refactor_modules_2026-10-31_native_output_identity.tsv),
[partial timings](data/refactor_modules_2026-10-31_partial_native_timings.tsv),
[background observation](data/refactor_modules_2026-10-31_background_observation.json)
and [interruption record](data/refactor_modules_2026-10-31_interruption.json).
The [previously prepared 0.4.0 bundle](release-0.4.0.md) predates this refactor.
