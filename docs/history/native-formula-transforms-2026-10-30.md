# Standalone formula-transform numerical kernels

This checkpoint corrects native formula transformations and prepares their
numerical operations for the [Rust-only wrapper migration](rust-only-wrapper-migration.md).
Its baseline is `deda996`, whose runtime is unchanged from `8fec875`. Package
versions remain 0.3.1. **R model-frame expression evaluation is still pending**:
the new standalone APIs do not by themselves replace that evaluator.

## Numerical changes

`design::numeric_transform` exposes finite model-covariate operations independently
of the formula parser:

- `polynomial_basis` constructs raw powers or centered orthogonal polynomials,
  returning values and reusable `alpha`/`norm2` recurrence metadata. Supplied
  coefficients evaluate an existing basis on new samples without fitting a QR.
- `scale_covariates` borrows row-major or column-major matrices, accepts automatic,
  disabled or explicit column controls, and returns applied centers and divisors.
- `integer_covariate` truncates toward zero and rejects values that R would turn
  into missing integers, including its reserved minimum integer.

The native formula parser uses these operations for its existing supported
grammar. Its existing 1–16 polynomial-degree syntax limit remains; the standalone
polynomial API has no such fixed limit. Matrix-valued scaling and reused
polynomial coefficients are available through the standalone API, with parser
integration still required. These are finite fitting-covariate operations:
non-finite outputs reject rather than presenting a general R vector API with
warnings/missing values.

The old `as.integer` numeric path was an identity operation. It now truncates
fractions and checks range. For example, adding `dose + dose +
as.integer(exposure)` for dose `[1,2,4,8]` and exposure `[.5,1,1.5,2]` now gives
`[2,5,9,18]`, agreeing with original R, instead of `[2.5,5,9.5,18]`.

Orthogonal polynomials previously used uncentered powers and one-pass
Gram–Schmidt. They now center before relative-norm column-deferral QR. The
Householder construction and explicit reduction schedule are shared with
polynomial factor contrasts from the [matrix checkpoint](rust-only-model-matrix-2026-10-30.md).
The scalar polynomial mean includes a correction pass and overflow fallback;
the scaling column mean uses its separate single-reduction operation. Neither
is selected by frontend, and the fitting rank tolerance is unchanged.

Raw powers retain R's signed-zero and small integer-power rounding behavior.
Scaling accepts finite negative supplied divisors, uses `max(1,n-1)` for its
automatic divisor, and avoids computing a mean when centering is disabled.
Thus identity scaling of large finite values no longer fails because an unused
sum overflows. The initial explicit zero for the column-mean reduction also
preserves negative-zero subtraction behavior.

## Independent references and existing coverage

The new frozen panel contains **85 original-R cases**, comprising **69 finite
successes and 16 invalid fitting-output decisions**. It covers irregular and
repeated inputs, covariates shifted by 1e12, cancellation, degrees through 16,
long vectors, numerical rank, fractional raw powers, recurrence reuse, scaling
controls and metadata, signed zeros and integer boundaries. Every successful
value and numerical metadata field is compared as exact IEEE-754 bits. Scaling
replays both storage layouts; supported cases also exercise the production
native formula parser.

The independent generator pins R 4.6.2 and DESeq2 1.52.0, never loads rsdeseq2,
and refuses to overwrite references. Run it into a fresh path:

```sh
Rscript scripts/references/generate_numeric_transform_reference.R /tmp/numeric-transform-reference.json
```

Fresh generation is byte-identical to the frozen panel, and the overwrite guard
is verified. The stated bit agreement uses the recorded R 4.6.2 / OpenBLAS
0.3.34 VORTEX reference environment; no additional platform CI ran locally.

The new references are additive. No original frozen oracle, numerical tolerance
or release manifest changes. An old hand-written offset assertion encoded the
integer-coercion bug and is corrected; an obsolete rejection of negative scaling
is removed. Independent references cover both corrected behaviors.

The full Rust suite passes **1,065 tests**, plus one backend test. One existing
opt-in test remains ignored. The checked R source package reports **Status OK**,
**1,933 assertions**, and no failures, warnings or skips. All **189** existing
standalone/wrapped identity cases pass. The full previous matrix panel of
**130 cases**, including **all 94 ordered-level counts**, retains exact values,
names and term assignments after sharing the QR implementation.

All **121 per-contrast and nine global GTEx precision gates** pass, retaining
the 664×/597× minimum target/global error reductions. All **22** result and
diagnostic exports are byte-identical to the preceding checkpoint. Root/backend
Clippy and Rust 1.97.1 checks, strict documentation, vendor consistency and crate
packaging/compilation verification pass. The checked R archive matches **287**
source members, allowing only generated DESCRIPTION fields/DCF whitespace;
all **453** crate source/test/example members match their current files.

C ASan/UBSan and Rust ASan checks also pass for the existing synchronous native
bridge. The R interpreter, installed dependencies and Rust standard library are
not instrumented; leak detection is disabled. This checkpoint adds no FFI entry
point or R numerical callback.

## Single-pass speed and memory comparison

One paired pass covers all eleven full 56,200-gene GTEx contrasts, including
the biological and random-split designs. Native tools and BLAS each use one
thread. Diagnostic exports are checked separately outside the measured fit/
result-export scope; every timed result matches its reviewed hash. Builds,
conformance checks and unrelated fits do not overlap measurement.

| Measurement | Before | Updated | Before / updated |
| --- | ---: | ---: | ---: |
| Native total elapsed, 11 contrasts | 35.565 s | 35.662 s | 0.9973× |
| Native maximum peak RSS | 302.094 MiB | 299.906 MiB | 1.0073× |
| R workflow total elapsed, 20 cases | 89.130 s | 88.196 s | 1.0106× |
| R workflow maximum peak RSS | 1126.625 MiB | 1127.469 MiB | 0.9993× |

Native elapsed time increases **0.27%**, with maximum peak RSS down **2.19 MiB**.
Native throughput decreases 0.27%; **no 2% optimization round is counted**.
Original DESeq2 takes 190.109 s across the same contrasts, giving a **5.331×**
DESeq2/updated-Rust process-time ratio. These are single-pass observations, not
statistical-equivalence or causal-attribution claims. No case is discarded or
repeated to improve the measurement.

All **twenty complete R workflow payloads are byte-identical**, with default
strict comparison and no origin-metadata exception. Every workflow remains
below 2 GiB. R elapsed falls **1.05%**, largely from the biggest rlog case
(22.336 to 21.268 s); its RSS rises 2.45 MiB. The large two-group results case
instead increases 0.173 s and 8.27 MiB. Maximum panel RSS rises 0.84 MiB.
Per-case costs remain recorded, alongside the earlier preparation checkpoint's
still-open memory increases; this comparison does not resolve that audit.
Panel maximum RSS values can belong to different cases.

See the [native timings](../data/numeric_transform_2026-10-30_native_timings.tsv),
[native comparisons](../data/numeric_transform_2026-10-30_native_summary.tsv),
[R workflows](../data/numeric_transform_2026-10-30_r_workflows.tsv),
[R stages](../data/numeric_transform_2026-10-30_r_stages.tsv),
[complete-payload identities](../data/numeric_transform_2026-10-30_r_payload_identity.tsv),
[GTEx gates](../data/numeric_transform_2026-10-30_gtex_precision.tsv),
[matrix sweep](../data/numeric_transform_2026-10-30_matrix_sweep.tsv) and
[validation record](../data/numeric_transform_2026-10-30_validation.json).

After measurement, a one-line API documentation correction makes the scaling
divisor's square-root placement explicit. No executable statement changes.
The final source archive is rebuilt/checked, and final native runtime artifacts
are verified byte-identical to those measured. The measured R package has the
same numerical source, differing only in that vendored documentation comment.

## Remaining scope

This is a native numerical foundation, not completion of formula-interface
identity. The R wrapper still calls `model.frame`, which evaluates expressions
using R. The complete expression evaluator, integration of native term ordering
with `ModelMatrixPlan`, frozen-local VST and optional posterior/s-value routes
remain required. C18's enabled-but-no-refit allocation, initial redundant native
adjustment and the final combined performance/memory audit also remain open.

No library is extracted or published in this checkpoint. The numerical module
boundary permits later extraction after both interfaces use the same evaluator
and its contract is stable.
