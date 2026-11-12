# Shared formula matrix assembly and expanded contrast weights

This checkpoint replaces production R `model.matrix` calls with a standalone
Rust matrix assembler. Expanded normal-shrinkage list contrasts also calculate
their mean-difference weights in Rust. The baseline is `6ace382`, whose native
and R numerical artifacts are unchanged from `af6e9b8`. Package versions remain
0.3.1. The wider [Rust-only migration](rust-only-wrapper-migration.md) remains
in progress: expression evaluation in `model.frame` still runs in R.

## Shared computation and interface behavior

`design::ModelMatrixPlan` consumes evaluated numeric/factor variables and formula
term codes. Numeric matrices are borrowed in either storage order. The plan
validates shapes, finite primitives, factor codes/levels, contrast matrices and
term codes. It constructs treatment, SAS, sum, Helmert, polynomial and explicit
matrix factor coding; expands interactions with the first variable's columns
varying fastest; and preserves multiplication order and signed zeros. There is
no rank test, rescaling or optimizer change in this operation.

The constructor retains declared unused levels, matrix-valued numeric terms,
backtick term labels, coefficient order and one-based term assignments. Without
an intercept, it promotes only the first participating factor/term cell to full
indicators. Empty formulas can yield a zero-column matrix before the fitting
layer's own validation. Optional response metadata omits a response-only term
while retaining numeric response interactions and the original term assignment
indices. Date storage is transported as numbers, like the previous R route.

The R adapter translates evaluated variables and term metadata, preserves sample
names and attaches contrast descriptions. Fitting, reduced designs, mutable-fit
validation, design-aware transforms and expanded shrinkage all call this shared
assembler. It also remains directly usable without R. Built-in contrast policies
above and explicit contrast matrices are supported; an arbitrary R contrast
function must be supplied as its explicit matrix rather than executed as an R
numerical callback. This compatibility constraint is distinct from the still
pending port of model-frame expression evaluation.

`contrasts::averaged_list_contrast` validates distinct, disjoint coefficient
indices and creates positive/negative group means. Either group can be empty,
but the union must be nonempty. The R expanded-shrinkage adapter aligns names
and calls that function; it no longer calculates `1 / length(group)` itself.
Ordinary result-list contrast weights retain their existing API semantics.

## Polynomial compatibility and independent references

Generic QR gave substantial high-level-count differences and was not accepted.
The final implementation uses relative-norm column deferral, rank-limited
Householder replay and a specified portable floating-point reduction schedule.
The [R QR source](https://svn.r-project.org/R/trunk/src/appl/dqrdc2.f) documents
the column-deferral policy. This basis-construction policy is separate from the
existing fitting rank tolerance; the latter is unchanged.

An initial sample including 64 and 95 levels matched after that implementation,
but the complete 2–95 sweep exposed four more cases: **53, 68, 85 and 91**.
Independent QR snapshots localized the differences to deferred-column ordering.
Dot products and norms already matched; the residual norm update needed the
specified fused operation. The correction matches **every output bit for all
94 supported level counts** on the R 4.6.2 / OpenBLAS 0.3.34 VORTEX reference
environment. Those four cases are now permanent regressions. No lookup table,
host-dependent frontend switch, R callback or BLAS import supplies the results.

The frozen panel contains **48 matrix cases** and **10 mean-contrast cases**.
It checks values, names and term assignments, covering additive/interaction/
nested/removed terms, numeric matrix terms, logical/character inputs, aliases,
offsets, explicit coding, unused levels, response terms and ordered factors
through 95 levels. Rust replays both input layouts; wrapper tests also retain
live R matrix-attribute comparisons. Malformed indices/codes/shapes and duplicate
or overlapping contrast selections reject before numerical use.

The independent generator pins R 4.6.2 and DESeq2 1.52.0, never loads rsdeseq2,
and refuses to overwrite references. Reproduce into fresh paths:

```sh
Rscript scripts/references/generate_model_matrix_reference.R /tmp/model-matrix-reference.json
Rscript scripts/references/generate_model_matrix_reference.R --all-levels /tmp/model-matrix-full.json
```

The second command emits the full sweep along with the common matrix cases.
The standalone conformance executable accepts each case directly as a JSON
input. For example, after building `evaluate_r_interface`:

```python
import json, pathlib, subprocess
reference = json.load(open('/tmp/model-matrix-full.json'))
case_file = pathlib.Path('/tmp/model-matrix-case.json')
for case in reference['cases']:
    case_file.write_text(json.dumps(case))
    result = json.loads(subprocess.check_output([
        'target/release/examples/evaluate_r_interface', str(case_file)]))
    assert result['names'] == case['expected_names']
    assert result['assignments'] == case['assignments']
    assert result['values'] == case['expected']
```

Original released fixtures, numerical tolerances and reference manifests are
unchanged. The high-count references certify the stated reference arithmetic
profile; additional platform CI has not run locally.

## Verification and source identity

The full Rust suite passes **1,063** tests, plus one backend test, followed by
the final test-only review's three targeted checks. The checked R
source package reports **Status OK**, **1,933 assertions**, and no failures,
warnings or skips. All **189** standalone/wrapped conformance cases pass. Runtime
guards reject R matrix/contrast calculations in the new cases alongside the
existing statistical-function guards; a positive control detects a prohibited
R model-matrix call. Formula-expression evaluation is outside this tested leaf.

Both Clippy checks, Rust 1.97.1 minimum-version checks, formatting, documentation,
vendor integrity and crate packaging/compilation verification pass. All 46
Python support tests are covered, including a separately enabled R-dependent
test. The checked R archive matches **284** current source members, allowing
only build-generated DESCRIPTION fields/DCF whitespace. All **448** packaged
crate source/test/example files match their working files.

Fresh C AddressSanitizer/UndefinedBehaviorSanitizer and Rust AddressSanitizer
checks pass for the synchronous bridge and its malformed-call/valid controls.
No pointer is retained across calls, and no R allocation/evaluation occurs while
Rust borrows input/output buffers. The R interpreter, installed R dependencies
and Rust standard library are not rebuilt with instrumentation; leak detection
is disabled. This records the bridge scope rather than whole-process safety.

All **121** GTEx per-contrast gates and **9** global gates pass, retaining minimum
target/global error reductions of 664×/597×. All **22** result/diagnostic exports
are byte-identical to the baseline. The final CLI hash equals the one used for
that replay, including after the last polynomial correction.

## Timing and memory

One paired native pass covers all eleven full GTEx contrasts, followed by one
paired pass over the same twenty R workflows as the preceding checkpoint.
No compilation, conformance or other fitting jobs overlap those measurements.
Complete-payload comparisons use the default zero-tolerance policy, with no
metadata migration exception. No unchanged timing candidate is rerun.

| Measurement | Baseline | Candidate | Before/candidate |
| --- | ---: | ---: | ---: |
| Native total process time, 11 contrasts | 35.56560 s | 36.01481 s | 0.98753× |
| Native maximum peak RSS | 302.01563 MiB | 302.04688 MiB | 0.99990× |
| R total process time, 20 workflows | 88.60918 s | 88.01108 s | 1.00680× |
| R maximum peak RSS | 1127.17188 MiB | 1126.64063 MiB | 1.00047× |

Native elapsed time increases **1.26%** and throughput falls **1.25%** in this
single paired observation. **No optimization round is counted.** The
`blocked_random_rep03` case changes 3.90371 → 4.54164 s and dominates the total
increase; it remains included. This observation is not discarded or rerun.
Original DESeq2 takes 188.94596 s, a **5.246×** DESeq2/Rust total-time ratio,
with 1764.56250 MiB maximum peak RSS. No native quantification speed improvement
is attributed to the new matrix assembler.

R total time falls **0.67%**, with all **twenty complete serialized payloads
byte-identical**, including fitted state and metadata. The full-size multi-result
control changes 7.92376 → 7.42265 s and contributes substantially to that total;
these one-pass timings do not establish a causal frontend speed gain.

| R workflow | Process s before → candidate | Peak RSS MiB before → candidate | Target stage s before → candidate |
| --- | ---: | ---: | ---: |
| Normal shrinkage, 1,200 × 16 | 2.21103 → 2.20805 | 861.79688 → 861.51563 | Shrinkage 0.116 → 0.115 |
| Ten separate results, 56,200 × 16 | 7.92376 → 7.42265 | 1127.17188 → 1126.64063 | Ten results 3.108 → 3.054 |
| Batched results, 56,200 × 16 | 6.15892 → 6.12337 | 1125.03125 → 1124.45313 | Ten results 1.785 → 1.786 |
| Two-group results, 56,200 × 16 | 4.77706 → 4.78460 | 1005.56250 → 1005.06250 | Ten results 1.394 → 1.376 |
| SF offsets, 56,200 × 16 | 4.79461 → 4.80456 | 899.04688 → 898.64063 | Fit 2.201 → 2.215 |
| Matrix offsets, 56,200 × 16 | 4.90582 → 4.85373 | 938.59375 → 937.96875 | Fit 2.220 → 2.234 |
| rlog, 2,048 × 64 | 5.23383 → 5.24563 | 747.10938 → 744.76563 | rlog 3.113 → 3.103 |
| rlog, 4,096 × 96 | 21.36570 → 21.38292 | 742.84375 → 739.87500 | rlog 18.984 → 18.983 |

Native maximum RSS rises **0.03125 MiB**. Maximum R RSS falls **0.53125 MiB**;
all twenty workflows remain below the registered 2 GiB budget. The preceding
checkpoint's observed large-batch and rlog memory increases are not considered
resolved by this small change against its accepted baseline. The final combined
performance/memory audit remains open, including this native timing cost.

The [native process comparisons](../data/model_matrix_2026-10-30_native_summary.tsv),
[all native timings](../data/model_matrix_2026-10-30_native_timings.tsv),
[R process metrics](../data/model_matrix_2026-10-30_r_workflows.tsv),
[R stages](../data/model_matrix_2026-10-30_r_stages.tsv),
[complete payload identities](../data/model_matrix_2026-10-30_r_payload_identity.tsv),
[full ordered-factor sweep](../data/model_matrix_2026-10-30_polynomial_full_sweep.tsv)
and [validation record](../data/model_matrix_2026-10-30_validation.json) preserve the
observations. Immutable runtime artifacts use `model-matrix-*` under
`results/codebase_improvements`.

After timing, the new generator's default/fullsweep modes were independently
reproduced: default reference bytes match and all 130 full-sweep/common matrices
match names, assignments and value bits. Test-only review removes a dead
approximate-comparison branch and adds invalid-response checks. The final crate
is repackaged/verified with those test changes; CLI, conformance executable and
checked R source remain byte-identical to their measured artifacts.

The main README and its original-DESeq2 comparisons remain unchanged. Earlier
checkpoint observations, including the preceding preparation checkpoint's
individual memory increases, remain part of the final combined audit.

## Remaining work

R still evaluates model-frame expressions, including numerical formula transforms.
The older standalone formula parser also needs to share the corresponding
canonical evaluation/assembly path. Frozen-local VST grid generation/function
evaluation and optional ashr/apeglm/posterior/s-value calculations still contain
R numerical work. Initial native row adjustments, C18's enabled-but-no-refit
replacement allocation and final combined speed/memory/source reconciliation
remain open. This checkpoint completes neither C19 nor the whole development
plan.
