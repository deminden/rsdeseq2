# Large-file refactoring follow-up

This refactor starts from `bfd0de8`, the standalone Rust filtering and
normalization checkpoint. It reorganizes the six requested large files while
preserving their numerical implementations, public APIs, fixtures and
assertions. Package versions remain 0.3.1.

## Organization

| Area | Before | Entry file after | Topic files | Largest topic |
| --- | ---: | ---: | ---: | ---: |
| Native LRT tests | 3,354 | 204 | 7 | 772 |
| Wald tests | 3,081 | 117 | 9 | 545 |
| R result tests | 1,789 | 240 | 5 plus entry | 521 |
| Builder workflows | 1,559 | 7 | 5 | 680 |
| Contrast handling | 1,388 | 18 | 5 | 358 |
| Dispersion trends | 1,660 | 32 | 7 | 568 |

The [inventory](../data/refactor_large_files_2026-10-29_inventory.tsv) records the
entry paths and file lengths.

Rust integration tests retain their executables and shared entry-file fixtures.
Child modules group fitting, contrast, formula, replacement and metadata cases.
Individual test function names are unchanged; their full names gain a module
prefix. R result tests are grouped by selectors, filtering, alignment, object
inputs and all-zero contrasts. Selecting `filter = "^results"` runs the whole R
group.

Builder entry points are private modules under
`core/implementation/builder/workflows/`, organized into default/named,
contrast, replacement, LRT and LRT-replacement workflows. Contrast handling has
separate request types, model-frame lookup, resolution, all-zero and alias
modules. Trend models and their fitting implementations stay together in
parametric, local and mean modules, with shared row selection, checked
arithmetic, the common fit type and numerical tests in separate modules.
Existing public paths re-export the same functions and types. Shared helpers
use visibility restricted to their parent implementation module. The R vendor
is regenerated from the Rust sources.

## Precision and coverage

The source comparison checks 181 production function/method bodies and their
signatures, plus 19 data definitions. Both sides receive the same Rust formatter
before comparison; documentation attributes, private helper visibility and
optional trailing parameter commas are excluded from signature/definition
comparison. Numerical operation bodies are identical. All 111 integration-test
and helper bodies, 15 trend unit-test bodies and 41 complete R test expressions
are unchanged, including fixtures and assertions.

The full Rust suite passes 1,030 tests, with the existing optional full-GTEx test
ignored in the normal suite. The R 4.6.2 source-package check passes all 1,406
assertions with no failures, warnings or skips. Clippy, minimum Rust 1.97.1,
strict Rust documentation, formatting, vendor synchronization and crate
packaging/verification pass. The final checked R archive matches all 226
packaged source files, allowing only R's generated DESCRIPTION fields and
whitespace normalization. Twenty standalone/wrapped primitive cases also
retain IEEE-754 bit equality; their R statistical-call and native-import guards
pass.

All **121 per-contrast and nine global precision gates** pass against the
frozen DESeq2 references with their existing bounds.
All 22 GTEx result/diagnostic exports and all 22 timed native results from the
11 candidate/baseline pairs match the checkpoint byte for byte. The frozen DESeq2 references
and original release hash manifest are preserved. Exact refactor comparison
uses the `bfd0de8` output manifest, which already includes the documented
canonical BH arithmetic correction. This comparison isolates file
reorganization from that earlier mathematical change.

## Native speed and memory

One sequential pass across all eleven 56,200-gene GTEx contrasts on Apple M5,
with one fitting/BLAS thread and native R 4.6.2, measured **36.756 s before**
and **36.706 s after**. The baseline/after ratio is **1.0014×**; total elapsed
time differs by **0.14%**. Median per-contrast ratio is 1.0069×. Maximum native
RSS was **302.34 MiB before** and **302.00 MiB after**.

The [complete comparison](../data/refactor_large_files_2026-10-29_speed_memory.tsv)
and [individual runs](../data/refactor_large_files_2026-10-29_timings.tsv) include
DESeq2 timing and memory. Tool order rotates across contrasts, and fits run
sequentially with compilation stopped. Timings include startup, input parsing,
fitting and result serialization. Diagnostic preflights are outside timing.
These single observations support preserved throughput and memory; they do not
establish repeated-run uncertainty or an optimization gain. The last contrast
has a slower candidate observation (3.110 versus 2.701 s); the aggregate remains
within 0.14%. No optimization round is counted.

The [provenance](../data/refactor_large_files_2026-10-29_provenance.json),
[acceptance record](../data/refactor_large_files_2026-10-29_acceptance.json) and
[checkpoint output hashes](../data/refactor_large_files_2026-10-29_output_sha256.tsv)
bind measurements to the two saved executables and unchanged references.

## Wrapped workflows

All **14 paired workloads** have exactly equal fitted state, result values and
metadata (`tolerance = 0`, including attributes). All fourteen serialized
payload files also have identical SHA-256 hashes. These cover Wald, LRT, small
and 56,200-gene repeated queries, normal/ashr/apeglm shrinkage, VST, actual
replacement, enabled replacement without refitting, and four rlog sizes up to
4,096 genes × 96 samples.

Their summed process times were **60.105 s before** and **60.349 s after**
(baseline/after **0.9960×**, a **0.41%** elapsed-time difference). Maximum RSS
was **1,320.67 MiB before** and **1,320.39 MiB after**, within the existing
2 GiB workflow budget. Individual single-pass fluctuations are retained in the
[workflow summary](../data/refactor_large_files_2026-10-29_workflows.tsv) and
[raw runs](../data/refactor_large_files_2026-10-29_workflow_runs.tsv); these timings
include R startup and package loading. The
[workflow provenance](../data/refactor_large_files_2026-10-29_workflow_provenance.json)
records the measured installed-package hashes. No fits or compilation overlap
the paired timing runs.

The [validation record](../data/refactor_large_files_2026-10-29_validation.json)
records test counts, source/body checks and artifact hashes. Detailed local
logs and checked packages are retained under
`results/codebase_improvements/refactor-*`.

## Scope

This checkpoint validates the requested refactor. The remaining
[Rust-only wrapper migration](rust-only-wrapper-migration.md), including the
remaining result and shrinkage calculations in R, is still in progress.
The [October 16 report](refactoring-large-files.md) retains evidence for the
earlier restructuring.
