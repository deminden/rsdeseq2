# Current-source test and contract audit

This local audit starts from `6e85f68`. It adds the
[supported-contract inventory](../native-supported-contract.md), corrects a stale
test expectation, and identifies an unresolved original-ashr comparison failure.
It does not close C19/C21/C23 or certify a release.

The later [sparse SQP correction](native-sparse-sqp-development.md) resolves the
numerical finding below in shared Rust and passes the existing test unchanged.
This report retains the failing checkpoint and its original evidence.

## Formula normalization assertion

The complete-formula correction in `a83221d` intentionally consumes sample
scales into the effective normalization matrix and canonicalizes unused sample
factor metadata to ones. Its independent original/wrapped checks already cover
that behavior. The older `row_preparation.rs` test still expected the tiny
supplied sample factors in the formula fit, even though its complete result-row
and normalization-matrix comparisons passed.

The test now requires ones for the complete formula workflow and retains an
additional assertion that an explicit matrix fit preserves supplied sample
factor metadata. Exact fitted rows and effective matrix comparisons remain.
No production numerical source, existing oracle or precision allowance changes.

## Unresolved sparse SQP comparison

The existing `test-fitted-view.R` compares normal, ashr and apeglm shrinkage on
18 reordered genes from a 40-gene fit against original DESeq2. Normal and apeglm
pass; ashr posterior means and SDs fail their existing `1e-5` test. Maximum
absolute differences are approximately **0.0620** and **0.0380** respectively.
This is a release blocker, not a new permitted tolerance.

The same differences occur when original `ashr::ash` and the native engine
receive exactly the same primitive estimates and SEs. The underlying original
and native unshrunk estimates already agree to about `3.9e-15`, and SEs to
`1.0e-15`. Fitted-object subset mapping therefore does not explain this failure.

Both engines create the same 13-scale normal-mixture grid. Original ashr selects
one component; the native engine retains approximately 0.7205563 and 0.2794437
on two adjacent components. On the same row-normalized likelihood matrix their
unregularized log likelihoods are approximately -2.549833 and -2.540524.
At the original boundary solution an inactive component has a negative
objective gradient of about -0.003775. The original positive-coordinate-only
convergence check does not inspect that improving inactive direction.

Independent original-package probes distinguish the paths:

- Original default `numiter.em=20` and automatic active-set budget 14 return the
  single-component boundary solution in two SQP iterations.
- Original `numiter.em=0` or explicit `maxiter.activeset=20` retains the two
  components and agrees with the native solution to the displayed precision.
- Native explicit budget 14 still returns the two-component solution, so merely
  substituting the original default budget is insufficient.
- Original `tol.svd=0` leaves the default boundary result unchanged. The original
  default trace also falls back to the full matrix in this case.

Inspection of pinned `mixsqp` 0.3-54 sources identifies relevant policy differences:
its EM assignment probabilities include `1e-15` floors, and its active-set
routine solves for a direction using the current reduced gradient. It tests
bound multipliers when the direction is near zero, and counts those transitions
against the iteration budget. The existing native quadratic routine solves for
target coordinates and checks inactive gradients after a completed move.
These observations identify the next implementation work; they do not yet prove
which arithmetic/transition difference causes the final boundary choice.

The native solution has a better objective here, but that is not evidence of
passing the established original-output compatibility contract. Preserve the
existing failing test and isolate the source iteration/regularization path
before changing solver behavior. Do not change L-BFGS-B, loosen numerical
bounds, add an R fallback, or force a component based on this fixture.

## Checks and evidence limits

- Full locked Rust workspace: **1,153 passed**, one existing opt-in test ignored.
- Python 3.14.7: **48 cases**, one skip, no failures. The first launch used the
  unsupported macOS Python 3.9; that interpreter failure is retained separately.
- The production-R lexical guard passes all **39 files**, with eight negative
  controls. Source formatting, whitespace and read-only vendor checks pass.
- All **631 runtime source paths** match the accepted immutable build. The
  current changes affect tests and documentation only.
- A correctly loaded installed-package R run of the affected files passes
  built-in/count likelihood imports, posterior primitives and native boundary
  checks; the fitted-view test retains the two ashr numerical failures above.
  The initial ad hoc broad R launcher did not expose private helpers or the
  package to child processes. Its harness errors are not production failures,
  and that run is not reported as a successful full R suite.

The [local report](../data/native-plan-audit-local.json) records source identities,
commands, logs and diagnostic artifacts. No new total-time/RSS measurement,
Linux oracle, package/platform/sanitizer gate, release bundle, optimization
round or publication is claimed. Resolve the numerical finding before the
remaining final combined checks.
