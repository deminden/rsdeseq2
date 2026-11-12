# Native formula parsing uses the shared numerical engine

This is a historical checkpoint. Its factor/class limitations were subsequently
addressed in the [0.4.0 candidate](../releases/0.4.0.md); the statements below
describe the earlier parser and its measurements.

The native string-formula parser now translates its numeric transforms and
matrix columns into the same Rust operations used by the primitive interface
and R wrapper. The parser retains its existing grammar, aliases, coefficient
labels and column order. Factor-transform/class semantics remain a separate,
unfinished part of the migration.

## Numerical corrections and consolidation

Validated scalar, binary and unary operations, integer powers, elementary
functions, numeric/factor coercions, scaling and polynomials call
`design::expression`. Model offsets use `expression::frame::offsets` with the
same ordered addition and positive-zero seed. Compatibility coefficient labels
are metadata; they no longer select a different arithmetic implementation.

This corrects independently demonstrated differences in the previous parser:

- Integer powers used `powi`, which disagreed with R's power operation on some
  finite values. The shared evaluator retains the independently checked R
  operation order.
- `log(x, base=2)` and `log(x, base=10)` now use the same specialized operations
  as the shared expression path. `log1p(-0)` returns R's positive zero.
- Repeated identical offset terms are counted once, matching R's term handling.
  A lexical primitive key normalizes whitespace, quoted symbols and equivalent
  decimal constants while retaining argument/operator structure. Distinct
  variables with identical values remain distinct offsets.

The expanded native matrix builder translates its existing column layout into
an internal `ModelMatrixPlan` constructor. Factor indicators and interaction
products are evaluated by the same plan as the R-style term interface.
Compatibility labels and interaction-column order are retained explicitly.
Treatment contrasts use implicit indicators, avoiding a new quadratic matrix
allocation when the compatibility plan only needs full factor indicators.
Other contrast policies retain their existing numerical construction.

## Independent and direct evidence

The new `tests/data/native_formula_r462.json` contains 54 original-R cases:
49 finite matrices/offsets and five rejected non-finite fitting inputs. It covers
irregular finite powers, signed zero, scalar/binary arithmetic, log bases,
coercions, raw/orthogonal polynomials, scaling, factor/numeric/higher products,
offset identity and duplicate handling. Two treatment cases use 129 declared
levels, including unused levels and a middle reference. They preserve the
separate existing 95-level limit for polynomial contrasts.

`scripts/references/generate_native_formula_reference.R` requires R 4.6.2 / DESeq2 1.52.0,
uses original `model.frame`, `model.matrix` and `model.offset`, and never loads
rsdeseq2. It writes only to fresh paths. Two regenerations are byte-identical;
the original 52 draft cases were retained when the two treatment cases were
appended. Existing frozen references and release manifests are unchanged.

`scripts/validation/check_native_formula_identity.R` compares the actual native parser with
the checked R package's Rust-backed frame/matrix/offset route and with the
independent reference. Matrices are compared bit for bit after an explicit
mapping of R column names to the compatibility parser's column order. Native
public labels are checked separately against frozen expected labels. This does
not claim that the two distinct formula metadata APIs expose identical names.

A control probe builds the production code from `6456690`, adding only the
test-only JSON driver. It disagrees with original R in 13 of the 54 cases:
44 matrix entries and 31 offset entries. The current candidate agrees on every
case and fitting-input rejection. The control's source archive and probe source
hashes are retained with its per-case differences; its production sources were
not edited.

The full Rust suite passes 1,071 tests; the backend passes three. The checked R
source passes 2,445 assertions with zero failures/warnings/skips. The previous
189 primitive-interface cases and 188 expression cases remain identical,
including their tested metadata. All 121 per-contrast and nine global GTEx
accuracy gates pass, and all 22 native result/diagnostic exports are byte-identical
to `6456690`. Clippy with warnings denied, Rust 1.97.1 checks for both graphs,
vendor integrity and all 46 Python harness tests pass.

A fresh disposable sanitizer build passes four formula and four replacement
unwind/recovery controls. C uses address/undefined-behavior instrumentation;
Rust uses address instrumentation. The R interpreter, preinstalled dependencies
and Rust standard library are not instrumented, and leak detection is disabled.
The ordinary checked library contains no test hooks or sanitizer symbols.

## Artifacts and measurements

The checked R source archive is SHA-256
`f9ab8b91b8ba6e452b866cdf3cb620e87ce643dd0105af46e71b196b1bf0be93`,
with 318 verified source members. The crate has 476 verified source/test/example
members. The native CLI is
`75c5e6884624970995e52a4bdca2fd4462f474b80e89d9717d5fbaef5e02fdce`;
the standalone driver is
`fb7cb44ce0daa39ac2b64f9492fb8fa6ecab26907f2e862d9ed2f113712f3f58`.

The [validation record](../data/native_parser_2026-10-30_validation.json) binds the
artifacts, checked source, test results, exact outputs and measured package
identities. One final paired pass compares the candidate with `6456690` across
all eleven GTEx contrasts and twenty R workflows. Earlier drafts received
correctness/build checks only, with no timing pass. No unchanged case was
retimed or discarded. The main README and release manifests remain unchanged.

The measurements use Apple M5/macOS 27.2 and `deseq2-r462-native`, R 4.6.2 /
DESeq2 1.52.0. Rust and BLAS each use one thread. Native timing includes process
startup, input parsing, the full Wald fit and export. RSS covers the entire
process, including libraries, input and native allocations.

| Paired panel | Previous elapsed | Candidate elapsed | Change | Previous maximum RSS | Candidate maximum RSS |
| --- | ---: | ---: | ---: | ---: | ---: |
| Native, 11 contrasts | 36.269 s | 35.782 s | −1.34% | 302.141 MiB | 302.063 MiB |
| R, 20 workflows | 88.022 s | 87.686 s | −0.38% | 1,128.078 MiB | 1,127.594 MiB |

Native throughput improves 1.362%, below the 2% total-quantification threshold;
**no optimization round is counted**. Original DESeq2 takes 189.625 s in this
pass, a 5.299× process-time ratio, with maximum RSS 1,877.000 MiB. All twenty
complete R payloads are byte-identical without a metadata exception, and every
R workflow stays below 2 GiB.

Slower observations remain in the record. The native blocked fourth split takes
longer for the candidate. R's largest rlog adds 0.063 s and 3.953 MiB; small
replacement adds 0.037 s and small LRT adds 0.028 s. The native panel supplies
explicit design matrices, so these aggregate measurements do not establish the
specific speed effect of formula parsing. The independent matrix/offset panel
provides that path's correctness evidence.

Tracked evidence:

- [Previous-parser disagreements](../data/native_parser_2026-10-30_previous_parser_disagreements.json).
- [Native summary](../data/native_parser_2026-10-30_native_summary.tsv) and
  [individual timings/RSS](../data/native_parser_2026-10-30_native_timings.tsv).
- [R workflows](../data/native_parser_2026-10-30_r_workflows.tsv) and
  [stage timings](../data/native_parser_2026-10-30_r_stages.tsv).
- [Native export identity](../data/native_parser_2026-10-30_native_output_identity.tsv)
  and [complete R payload identity](../data/native_parser_2026-10-30_r_payload_identity.tsv).
- [Per-contrast precision](../data/native_parser_2026-10-30_gtex_precision.tsv),
  [global precision](../data/native_parser_2026-10-30_gtex_global_precision.tsv) and
  [R session](../data/native_parser_2026-10-30_R-session.txt).

## Reproduction

```sh
cargo build --release --locked --example evaluate_r_interface
Rscript scripts/references/generate_native_formula_reference.R /fresh/path/native_formula.json
Rscript scripts/validation/check_native_formula_identity.R STANDALONE_DRIVER R_LIBRARY FRESH_OUTPUT_DIR
cargo test --locked --test native_formula --test numeric_transform --test model_matrix
```

Use the pinned original-R environment for reference generation and the installed
candidate library for the direct comparison. The existing full GTEx and paired
R workflow harnesses provide the unchanged whole-workflow acceptance gates.
Install the previous revision and candidate into separate libraries for paired
replay. The validation record preserves their measured file hashes. After both
panels completed, the development library was updated from the checked archive;
installed-package frame/offset identity and custom-call rejection smoke checks
passed.

## Remaining scope

The compatibility parser still needs shared factor-transform/class handling.
In particular, native `ordered()`/`as.ordered()` currently use treatment coding,
and nested `relevel(ordered(...))` is incorrectly accepted. Default factor
handling of unused levels and broader argument/method/warning semantics require
their independent parity checks and explicit corrections. The present numeric
and matrix consolidation does not certify those paths or complete C19.

Frozen-local VST, optional ashr/apeglm posteriors/s-values, the initial redundant
adjustment, C21's authoritative result path and C23's final combined audit remain
open in the [migration plan](rust-only-wrapper-migration.md). Custom formula
functions remain rejected before execution.
