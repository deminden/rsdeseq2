# Sparse ashr SQP compatibility correction

The [combined-suite audit](native-plan-audit-development.md) found a real
default-ashr discrepancy in the existing reordered-gene subset test. It is
corrected in the shared Rust engine, without an R numerical fallback or a
precision-bound change.

## Cause and implementation

Original mixsqp 0.3-54 performs floored EM assignment updates, then solves
direction-based active-set quadratic subproblems. Working-set transitions
consume the iteration budget; the default budget is `min(20, components + 1)`.
Its outer convergence check examines positive coordinates only. This can accept
a boundary solution despite an improving inactive component. A mathematically
better objective from a different iteration path therefore does not guarantee
original-output compatibility.

The earlier native implementation solved for target coordinates and checked
inactive multipliers after a completed move. It also lacked the source's
condition-based approximate linear solve. A successful Cholesky factorization
was treated as sufficient even when the reduced system was nearly singular.

`MixtureSqpOptions::ashr_compatibility` now selects the source EM floors,
gradient/Hessian operation order and active-set transition policy. Complete
location-scale shrinkage and the shared ashr control parser enable it by
default. Explicit active-set budgets remain supported. The general convex
fitter keeps its existing target-coordinate and full-KKT policy; both interfaces
use the same selected core options.

The direction solver calculates a reciprocal 1-norm condition estimate from its
Cholesky factor. Below machine epsilon, it uses a Rust SVD with the original
dimension-scaled singular-value cutoff. No R, Armadillo, BLAS or LAPACK runtime
is called for this solve. An independent nearly rank-one system checks the
distinction between positive definiteness and a trustworthy full-rank direction.
Ridge adjustment has a finite failure path rather than a non-progressing loop.
The MIT mixsqp notice ships with the core, R package and binary notices.

## Numerical evidence

The existing public 18-gene subset test retains its `1e-5` bound and now passes.
For identical primitive estimates/SEs, posterior mean and SD differences from
original ashr fall from approximately **0.0620 / 0.0380 to zero** on this case.
The public fitted workflow, whose unshrunk estimates differ from original by
about `3.9e-15`, has maximum mean/SD differences of approximately
**8.3e-16 / 3.1e-16**.

The additional original-only generator records nine complete and nine primitive
profiles: default/no-EM, explicit 14/20 and automatic active-set budgets, a short
SQP run, reversed observation order and full-data controls. All preserve the
existing SQP `2e-6` numerical bounds and original convergence decisions.
The complete original warning list was not asserted at this checkpoint. The
later [warning/platform follow-up](native-sqp-warning-platform-gates.md) detects
and corrects a missing mixsqp nonconvergence warning, and adds that assertion.
Maximum complete-profile mean/SD/s-value errors are about
`7.2e-15 / 3.3e-15 / 2.1e-15`. Six 256-gene public control profiles and three
1,200-gene default/sign/threshold workflows also pass their existing original
comparisons with R calculations blocked; the latter mean/SD errors remain below
`2.3e-14 / 4.4e-14`.
The generator loads original DESeq2/ashr/mixsqp only. It recreates the seeded
40-gene input independently of candidate package helpers. Previously committed
reference files are unchanged.

All 29 prior replay sets are retained. Nineteen remain completely exact,
including 31 matrix-fit requests/113 later queries, 33 formula requests/96 later
queries, general and beta-binomial regression, random starts, posterior analysis,
EM/cycling/interior-point methods and scalar likelihood kernels. The other ten
contain intentional SQP numerical/diagnostic changes; all replay warning lists
remain identical. The old immutable executable and complete before/after
payloads remain available. These results do not claim every previous SQP bit is
unchanged. All 29 sets preserve field sets, shapes, acceptance, finite/missing
patterns and categorical diagnostics, in addition to warning lists.

## Reproduction and scope

```sh
Rscript --vanilla scripts/generate_sparse_sqp_reference.R "$FRESH_REFERENCE_JSON"
cargo test --locked -p rsdeseq2 --test sqp_controls
cargo build --release --locked -p rsdeseq2 --example evaluate_r_interface
Rscript --vanilla scripts/check_location_shrinkage_identity.R \
  target/release/examples/evaluate_r_interface "$R_LIBRARY" "$FRESH_LOCATION_DIR" \
  crates/rsdeseq2/tests/data/sparse_sqp_r462.json
Rscript --vanilla scripts/check_native_ashr_optimizers.R \
  target/release/examples/evaluate_r_interface "$R_LIBRARY" "$FRESH_DIAGNOSTICS_DIR" \
  - crates/rsdeseq2/tests/data/sparse_sqp_r462.json
```

Both R CI lanes require the new complete-output and optimizer-diagnostic bit
checks with R calculations blocked. The dash skips cached public workflows;
the primitive comparison still runs. Independent original accuracy and actual
standalone/wrapped identity establish different requirements.

The full locked Rust workspace passes **1,154 tests**, with one existing opt-in
test ignored. The rebuilt R source package has **Status: OK**, with **2,867
assertions**, no failures, warnings or skips. Its 430 files match the working
tree, allowing only standard DESCRIPTION formatting/build metadata. Actual
interface checks and source identities are recorded in the
[local report](data/native-sparse-sqp-local.json).
The Student-t canonical syntax fingerprint is also constructed without implying
an undeclared numerical-package import; it remains syntax inspection only.

This closes the sparse SQP finding, not the final combined C19/C21/C23 audit.
Current total speed/RSS, frozen GTEx/held-out replay, required platform gates,
boundary sanitizers and final release-bundle reconciliation remain required.
No dependency version, L-BFGS-B/locfit implementation, optimization round or
publication changes are part of this correction.
