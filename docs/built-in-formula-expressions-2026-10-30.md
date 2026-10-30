# Built-in formula expressions and custom-function policy

The user elected to disallow custom formula functions for now, superseding the
earlier request to translate custom R function bodies. Formula calls now undergo
a non-executing check before any model-frame evaluation. Unknown functions,
computed calls, arbitrary namespaces, aliases and bindings that mask the
supported built-ins fail with an actionable error. Active bindings and delayed
definitions are rejected without forcing them. Explicit `base::` and
`stats::` qualification resolves the established function. This policy applies
to full/reduced fits, fitted-state validation, transformations, normal-shrinkage
preparation and direct model-matrix construction. A caller can supply a numeric
design matrix after preparing a custom transformation outside rsdeseq2.

## Standalone evaluator foundation

`design::expression::evaluate` accepts owned primitive bindings and an expression
tree. It has no callback, R pointer, R runtime or executable closure. Numerical
and categorical storage retains matrix shape, names, factor codes/levels,
polynomial coefficients and scaling metadata. The supported built-in subset
includes nested arithmetic, scalar math, coercion/missingness, conditional
selection, matrix construction, factors, scaling and univariate polynomials.
Custom calls are rejected even inside an unused conditional branch. Dimensions,
integer codes, metadata and recursion depth are validated before indexing,
including metadata constructed by a nested built-in. This is an unreleased
development interface; its checked subset is not a general R interpreter.

The evaluator reuses the standalone polynomial, scalar-mean and scaling kernels.
Scaling has one arithmetic implementation: the existing finite fitting API
rejects unusable final covariates, while expression evaluation permits missing
intermediate values that a subsequent built-in may replace. The previous 85
independent finite-transform references remain unchanged. Conditional selection
evaluates only needed branches and preserves R's selected output type and matrix
recycling order for the tested atomic cases.

Factor ordering needs collation metadata. ICU4X `icu_collator` and `icu_locale`
2.4.0 provide standalone compiled collation data; no R library is loaded.
These were the latest stable releases checked for this work (published
2026-10-22, declared MSRV 1.88). Both root and R-backend lockfiles carry the same
additions. The default is Unicode root collation; the API also accepts `C` byte
ordering and Unicode locale identifiers. This does not certify every host R
locale or ICU version. See the official
[collator](https://docs.rs/icu_collator/2.4.0/icu_collator/) and
[locale](https://docs.rs/icu_locale/2.4.0/icu_locale/) documentation.

## Independent checks and current limits

`scripts/generate_formula_expression_reference.R` pins R 4.6.2 / DESeq2 1.52.0,
does not load rsdeseq2 and refuses to overwrite a reference. Its 129 cases (114 successes and 15 rejections) include
negative controls, missingness, NA/NaN distinctions, integer overflow, nested
transforms, declared/Unicode factor levels, recurrence reuse, conditional branch
selection and metadata. Successful outputs are compared by every double bit,
storage type and retained attribute. The original 86 draft expectations were
preserved as cases were appended. Regeneration produced identical files.

This is an evaluator foundation, **not completion of C19**. Production R still
uses `stats::model.frame` and `stats::model.offset` behind the guarded adapter.
The new evaluator has not yet replaced that numerical path or the older native
formula parser. The R call guard is metadata validation, not numerical interface
identity evidence. The existing 189 interface cases cover earlier shared
operations, not this unconnected evaluator.

Before integration, finish the owned variable-shape output bridge with unwind
cleanup, common frame/offset assembly, legacy-parser adoption, multivariate and
non-finite raw polynomial semantics, warning propagation and broader argument/
class/locale handling. The fixture generator intentionally captures values while
muffling R warnings; matching these fixtures does not establish warning parity
or general R expression semantics. Frozen-local VST and optional posterior/
s-value migration also remain open in the
[Rust-only migration plan](rust-only-wrapper-migration.md).

To regenerate independently, supply a fresh output path:

```sh
Rscript scripts/generate_formula_expression_reference.R /fresh/path/expressions.json
cargo test --locked --test formula_expression --test numeric_transform
```

Use the pinned comparison environment (`deseq2-r462-native`) for generation.
The generator refuses an existing output path; its overwrite control was also
checked without modifying the original file. Existing release manifests and
numerical tolerances are unchanged.

## Packaged validation

The final R source package passes `R CMD check` with status OK: 2,014 assertions,
zero failures, warnings or skips. Root Rust passes 1,067 tests (one existing
opt-in replay remains ignored), the backend passes two, and all 46 Python
harness tests pass with the required Python/R runtimes. Root/backend Clippy deny
warnings, Rust 1.97.1 checks pass, and the vendor check is clean. No bridge ABI or
C ownership code changes in this checkpoint; the earlier C18 sanitizer evidence
is retained rather than claimed as a new sanitizer run.

All 189 existing standalone/wrapped identity cases pass against the final checked
package. The eleven full GTEx contrasts pass the 121 per-contrast and nine global
accuracy gates; all 22 result/diagnostic exports remain byte-identical to
`c181885`. This is evidence that existing quantification remains unchanged; the
unconnected expression evaluator's additional 129 cases are a separate check.

The checked archive contains 300 verified source members, and the frozen crate
contains 463 verified source/test/example members. Artifact hashes and the full
paired records are in [validation data](data/formula_expression_2026-10-30_validation.json).
The final R archive is `e7a0b9d52a588bd49e383ab8c34cdb8bac97d4a4338487c4a69e0be7f8b587d3`.
The final CLI is `b9ff7c0bdd6bb9c6ff7c753c1c703a123d8fc3899dc1d1fde71b96aa83a7c223`.

## Paired measurements

One final paired pass compares `c181885` with the corrected frozen build, using
one worker and one BLAS thread. Each study/workflow runs once per tool:

| Panel | Baseline total | Candidate total | Baseline maximum RSS | Candidate maximum RSS |
| --- | ---: | ---: | ---: | ---: |
| Native, eleven GTEx studies | 36.410119 s | 35.710168 s | 301.96875 MiB | 301.875 MiB |
| R, twenty workflows | 89.244561 s | 88.350224 s | 1,141.65625 MiB | 1,127.375 MiB |

Native elapsed time falls 1.92241% (throughput +1.96009%). This is below the 2%
criterion, so no optimization round is counted. Original DESeq2 totals
192.604199 seconds and peaks at 1,764.421875 MiB; its process-time ratio to the
candidate is 5.39354×. All 22 native exports and all twenty complete raw RDS
payloads are byte-identical, with no metadata-migration exception.

R aggregate elapsed time falls 1.00212%. Its maximum RSS falls 14.28125 MiB,
largely in the large multi-results case. Keep the individual costs: the large
no-refit case rises 945.78125 → 946.328125 MiB (+0.546875), the small multi-results
case rises 2.6875 MiB, and the largest rlog case takes an additional 0.05533 seconds.
Every case stays below the 2 GiB budget. A single pair records the observed
costs; it does not isolate the causal effect of each component or prove a
small speed/memory delta. The new evaluator itself is not in these production
R workflows yet. Full [R rows](data/formula_expression_2026-10-30_r_workflows.tsv)
and [native timings](data/formula_expression_2026-10-30_native_timings.tsv) retain
all cases rather than selected favorable ones. The main README is unchanged.

Earlier observations are retained separately: the first native pair totals
36.155928 → 35.922305 seconds (throughput +0.65036%), with maximum RSS
301.953125 → 302.03125 MiB. The first R guard pair totals
90.626819 → 88.198230 seconds / maximum RSS 1,153.0625 → 1,126.84375 MiB.
All twenty raw payloads are exact. Subsequent binding-lookup hardening received
its own R-only pair, also with twenty exact payloads. Their complete rows are in
[first native timings](data/formula_expression_2026-10-30_first_native_timings.tsv),
[first R guard](data/formula_expression_2026-10-30_first_r_guard_workflows.tsv) and
[hardened R guard](data/formula_expression_2026-10-30_hardened_r_guard_workflows.tsv).

The last standalone correction validates constructed output metadata before a
parent indexes it; two independent malformed-matrix rejection cases were added
without changing the previous 127 expectations. The rebuilt artifacts differ,
so the final corrected build receives its own single native/R pair. No unchanged
candidate or outlier is rerun to select a better result. The earlier observations
are superseded for final artifact acceptance, not discarded. None of these
checks establishes full Rust-only formula integration or a general R evaluator.
