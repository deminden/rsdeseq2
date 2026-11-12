# Shared Rust formula frames and offsets

The wrapper now evaluates formula expressions, selects missing-value rows,
prepares categorical/response values and sums model offsets through the
standalone Rust engine. Production `.R` files no longer call
`stats::model.frame`, `stats::model.offset`, or construct character/logical
factors numerically in the model-matrix adapter. Custom formula functions remain
rejected before evaluation. No R numerical fallback or callback is introduced.

R retains symbolic `terms()` parsing, variable lookup, identifier alignment,
prediction-call construction from returned coefficients, class labels and
Bioconductor object construction. These steps translate data and metadata;
they do not recompute transformed values, scaling parameters or offsets.

## Shared operations

`design::expression::evaluate_many` initializes bindings/collation once.
`frame::evaluate` combines expression evaluation and row policy in one native
transaction. `frame::prepare_model_values` performs character/logical factor
coding and response conversion; `frame::offsets` uses the same ordered addition
and scalar-zero seed for standalone and wrapped calls. An absent offset returns
NULL by metadata inspection without a numerical round trip.

The independent extension adds integer/variadic sums, trimmed means,
character/factor/ordered comparisons, named concatenation, multivariate
polynomials and recurrence reuse. Raw polynomials permit non-finite intermediate
values, using the same arithmetic as the existing finite fitting API; that API
continues to reject invalid final covariates. Existing finite polynomial/scaling
references and fitting rank thresholds are unchanged. Integer carriers have one
missing sentinel and no signed zero, while ordinary doubles retain their bits.

Missingness is checked across every matrix column. Omission/exclusion also
reproduces R's full-row subsetting of attributed columns when no rows are
missing. This matters because subsetting drops names attached to dimension-label
vectors. Multivariate polynomial labels carry grid-row names; the primitive
interface now retains those separately. Named matrix elements follow row-major
storage internally and are reordered only when marshalled to/from R.

## Ownership and bounds

The versioned RFX2 transport contains only owned primitive data and metadata.
C completes input allocation, UTF-8 translation and ALTREP materialization
before entering Rust. Rust decodes owned values and returns an opaque owner of
encoded output bytes. It retains no R address and executes no R callback.

The R unwind continuation exists before Rust publishes an owner. After the Rust
frame returns, C allocates output objects under `R_UnwindProtect`; its cleanup
consumes and clears the owner on normal return and errors. No persistent R
external pointer or finalizer is needed. Native tests check publication only on
success, immutable byte queries, nonempty-slot rejection and idempotent drop.

Transport payloads and expression-created arrays are bounded. Transport permits
1,024 structural levels so list/argument wrappers do not reduce the evaluator's
128 semantic levels. The isolated R boundary tests accept 128 nested expression
levels and reject 129. These limits prevent a small symbolic request from
requesting an enormous sequence, replication, matrix or polynomial allocation;
they are separate from the representative-workflow memory budget.

## Evidence

The original 129 expression and 85 finite-transform references remain unchanged.
The new `formula_expression_extended_r462.json` contains 59 independent R 4.6.2 /
DESeq2 1.52.0 cases; regeneration produces identical files and existing paths
cannot be overwritten. Numeric calculations in the generator use original R,
not rsdeseq2. A shared test-only encoder records storage types, every double bit
and retained attributes, including names on dimension labels.

`scripts/validation/check_formula_expression_identity.R` compares the actual C/R bridge
with the standalone Rust conformance executable on all 188 cases, then checks
both against the independent references. Only JSON object field order is
canonicalized. No numerical tolerance or missingness exception is used.
Whole-frame tests additionally compare thirty live formulas under all four
standard missing-value actions, with `num.eq=FALSE` and `single.NA=FALSE`.
A runtime control blocks R numerical transforms/frame functions while the
wrapper constructs matching frames, matrices and offsets.

The checked candidate passes 1,069 root Rust tests, three backend tests and
2,445 R assertions, with zero R failures/warnings/skips. All 189 earlier direct
interface cases remain identical. The full eleven-contrast GTEx sweep passes
121 per-contrast and nine global accuracy gates; all 22 result/diagnostic exports
remain byte-identical to `19b5db4`. Clippy, Rust 1.97.1 checks for both graphs,
the vendor check and all 46 Python harness tests pass.

The disposable sanitizer build verifies C address/undefined-behavior and Rust
address instrumentation. Four forced formula allocation/decode errors release
their owners and recover exactly; the four existing replacement controls also
pass. Test hooks/counters are absent from the ordinary library. This does not
claim R/interpreter/dependency/standard-library instrumentation or complete heap
leak profiling. Leak detection remains disabled for the R process.

Earlier drafts are retained in `results/codebase_improvements`: one sanitizer
control initially compared formulas created in different environments, so its
positive recovery fixture was corrected to reuse one formula; source checks
also exposed named-label subsetting and integer/NULL corner cases. No timing
measurements were taken from those drafts, and no reference bounds were relaxed.

## Artifact and measurement record

The final source archive has 315 verified members and the crate has 470 verified
source/test/example members. The source SHA-256 is
`c38f279198235b89eba079feac680ac7d3e3df46382edc2a08fd4ed2d9e4a2a1`;
the CLI is `8e531702ee2ab73e0a1464233c7fd9da82eb66d974851f98439d352540ad4617`.
The [validation record](../data/formula_frames_2026-10-30_validation.json) binds
these artifacts to the checks, exact exports and measured package identities.
The main README and release manifests are unchanged.

One paired eleven-contrast native pass and one paired twenty-workflow R panel
compare the final candidate with `19b5db4`. There were no earlier timing drafts
for this checkpoint, and no unchanged case was rerun to replace a slower result.
The reference environment is `deseq2-r462-native`, R 4.6.2 / DESeq2 1.52.0 on
Apple M5/macOS 27.2; native Rust and BLAS each use one thread. Native timing
includes process startup, parsing, the full Wald fit and result export. RSS is
the maximum for the entire process, including libraries and input data.

| Paired panel | Previous total | Candidate total | Elapsed change | Previous maximum RSS | Candidate maximum RSS |
| --- | ---: | ---: | ---: | ---: | ---: |
| Native, 11 contrasts | 36.333 s | 36.270 s | −0.17% | 303.172 MiB | 302.359 MiB |
| R, 20 workflows | 90.708 s | 88.772 s | −2.13% | 1,128.484 MiB | 1,127.781 MiB |

Native throughput changes +0.173%, below the 2% total-quantification threshold;
**no optimization round is counted**. Original DESeq2 takes 190.333 s in the
same native panel, a 5.248× process-time ratio, with maximum RSS 1,764.453 MiB.
All native result/diagnostic hashes pass and all twenty complete R payloads are
byte-identical, without a metadata exception. Every R workflow remains below
the 2 GiB budget.

The aggregate timing differences are not evidence of a broad speed gain.
The native blocked fourth split has a slower baseline pair (4.555 → 3.964 s),
which dominates its small aggregate difference. The R large no-refit case is
8.567 → 6.507 s; the other nineteen cases together are 82.141 → 82.265 s
(+0.15%). Slower R controls include large multi-results (+0.171 s, +2.28%),
Normal shrinkage (+0.090 s, +4.11%) and small LRT (+0.079 s, +3.78%). Individual
memory costs also remain recorded: Normal shrinkage +4.797 MiB, large no-refit
+2.563 MiB and small Wald +2.047 MiB. These observations are retained for the
final combined audit; no case is discarded or retimed.

Tracked measurement evidence:

- [Native summary](../data/formula_frames_2026-10-30_native_summary.tsv) and
  [individual timings/RSS](../data/formula_frames_2026-10-30_native_timings.tsv).
- [R workflows](../data/formula_frames_2026-10-30_r_workflows.tsv) and
  [stage timings](../data/formula_frames_2026-10-30_r_stages.tsv).
- [R payload identities](../data/formula_frames_2026-10-30_r_payload_identity.tsv)
  and [native export identities](../data/formula_frames_2026-10-30_native_output_identity.tsv).
- [Per-contrast precision](../data/formula_frames_2026-10-30_gtex_precision.tsv),
  [global precision](../data/formula_frames_2026-10-30_gtex_global_precision.tsv)
  and [R session](../data/formula_frames_2026-10-30_R-session.txt).

Reproduction uses the existing full-contrast and workflow harnesses, comparing
the frozen `19b5db4` binaries/library with this checked source package. Reference
generation and direct identity checks are available without a benchmark run:

```sh
Rscript scripts/references/generate_extended_formula_reference.R /fresh/path/extended.json
Rscript scripts/validation/check_formula_expression_identity.R STANDALONE_DRIVER R_LIBRARY FRESH_OUTPUT_DIR
cargo test --locked --test formula_expression --test numeric_transform
```

The extended generator writes to a fresh path and refuses to replace an existing
reference. It requires the pinned original R/DESeq2 environment. The identity
script accepts the standalone `evaluate_r_interface` executable, an installed
candidate R library and a fresh output directory; it does not evaluate custom
function bodies. The new transport and frame tests run through the ordinary
Rust and R package suites, while
`scripts/validation/check_native_boundary_sanitizers.py` supplies the isolated instrumented
boundary check.

Replay should install the previous revision and candidate into separate R
libraries. The validation record preserves their measured file hashes. After
measurement, the development library was updated from the checked candidate
archive; an installed-package frame/offset identity and custom-call rejection
smoke check passed.

## Still open

This completes the production frame/offset numerical migration for the checked
primitive interface, not all of C19. The older native formula-string parser must
adopt the shared evaluator/term plan. General warning delivery and broader R
method/argument/class/collation semantics still need their parity audit; the
reference generator intentionally muffles R warnings to capture values. These
tests do not certify a general R interpreter or every ICU global configuration.

Frozen-local VST function/grid evaluation, optional ashr/apeglm posterior and
s-value calculations, initial redundant native adjustment, final shared-result
authority and the combined performance/artifact audit remain open in the
[migration plan](rust-only-wrapper-migration.md). No full Rust-only package,
release or completed-development-plan claim is made here.
