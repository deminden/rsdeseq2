# Standalone mixture-weight fitting

After the [posterior foundation](posterior-primitives-2026-11-10.md), the shared
Rust engine can fit nonnegative mixture proportions from a matrix of component
log likelihoods, observation weights and Dirichlet parameters at least one.
The primitive R transport calls this same implementation. Public `lfcShrink`
still uses its existing optional packages: this checkpoint supplies a fitting
kernel, not the complete ashr/apeglm method and option adapters.

## Method and numerical contract

The implementation uses EM warmup followed by nonnegative sequential quadratic
programming (SQP). Each quadratic subproblem uses an active set and a Cholesky
solve of the free coordinates; a small diagonal correction handles numerical
rank deficiency. A feasible backtracking step updates the weights. This follows
the mathematical approach described by
[Kim et al.](https://arxiv.org/abs/1806.01412), implemented directly in Rust with
the existing nalgebra dependency. It uses no R callback or R math runtime and
does not change the GLM optimizer or the L-BFGS-B dependency.

Log-likelihood rows are normalized by their maximum to retain small likelihoods.
The Dirichlet penalty supplies diagonal pseudo-observations. Zero-weight data
rows do not contribute, including rows whose component log likelihoods are all
negative infinity. Invalid shapes, weights, priors and controls return errors.
Nonrepresentable objectives, gradients, Hessians or solutions fail explicitly.
Both borrowed matrix layouts use the same logical arithmetic order.

The standalone defaults fit an unregularized likelihood (`likelihood_offset=0`).
The compatibility panel explicitly supplies `1e-6`, the pinned ashr mixSQP
adapter's offset inside row-mixture logarithms. The offset also applies to
Dirichlet pseudo-observations. Warmup, iteration limits, dual tolerance, active
iterations and the zero threshold are explicit options. Supplied initial weights
are validated and normalized; the future ashr adapter must also preserve its
own initialization preprocessing. A zero-likelihood initial population fails.

The reported dual residual is one-sided dual feasibility before final weight
normalization. The optimizer's weight sum is returned separately. This preserves
the distinction between the stopping criterion and an absolute stationarity
norm, especially with a nonzero likelihood offset. An iteration limit or failed
line search remains visible as `converged=false` with an explanatory message.
No tolerance or oracle was widened to turn a trial into a passing comparison.

## Finding from the rejected trial

An initial L-BFGS-B trial optimized the unregularized convex objective. It
achieved almost identical likelihoods, but one shifted-prior posterior differed
by `5.354e-6`, beyond the trial's predeclared `2e-6` bound. Adding the original
likelihood offset and refining that optimizer reduced the difference to
`2.165e-6`, still outside the same bound. Neither trial was accepted.

The offset and the original EM/SQP stopping procedure both matter for matching
the established output. The native SQP implementation passes the unchanged
references and bounds, with maximum observed mean/SD differences of
`1.110e-14`/`9.479e-14`. The unsuccessful trials and their logs remain in the
ignored working evidence and private optimization journal. No production fit
or public shrinkage route was switched to either trial.

## Independent accuracy and direct identity

The [reference generator](../../scripts/references/generate_mixture_fit_reference.R) uses only
original R 4.6.2 and ashr 2.2-63 with mixSQP's established default controls.
Its [frozen panel](../../crates/rsdeseq2/tests/data/mixture_fit_r462.json) contains
eight prior-fitting profiles and **1,920 observations**: null-dominated,
mixed signals, nonuniform/zero observation weights, a null-biased Dirichlet
penalty, shifted components, duplicate and near-collinear components, and
duplicate components with positive fitted mass. Local regeneration is
byte-identical.

The [Rust tests](../../crates/rsdeseq2/tests/mixture_fit.rs) compare weighted
likelihoods and posterior means/SDs against those independent fits. They retain
the original trial bounds of `1e-5` on weighted log likelihood and `2e-6` on
posterior moments. Analytic categorical/Dirichlet and prior-only modes check
known solutions; other cases verify malformed inputs, zero-weight exclusions,
visible iteration limits and exact row/column-major equivalence. Component
weights can be non-identifiable, so independent tests compare likelihoods and
posteriors rather than requiring a particular split between duplicate columns.

The [interface check](../../scripts/validation/check_mixture_fit_identity.R) requires every
fitted weight, likelihood and diagnostic bit to match the R-free executable and
actual native bridge. R `sqrt`, `exp`, `log`, `mean`, `sum`, `optim` and `dnorm`
are blocked during bridge calls, with a positive guard control.

An initial identity attempt exposed 1,901 input-bit mismatches from reparsing
17-digit fixture strings through two decimal readers. The corrected check
exports the actual R input bits, including controls, to the standalone process.
It preserves the independent decimal fixture and accuracy tests. This establishes
identical engine input rather than mistaking decimal-reader differences for a
different statistical kernel. All eight complete fits now match every output
bit. This transport correction does not change a production numerical path.

The [validation record](../data/mixture_fit_validation_2026-11-10.json) binds the
source tree, generator/fixture and checked R source/library and standalone driver.
Local checks pass: 1,083 Rust tests, three backend tests, 2,659 R assertions with
`R CMD check` status OK, 47 Python tests with one platform skip, both Clippy and
MSRV 1.97.1 graphs, and vendor integrity. All 328 archived R/C/Rust files match
the checked source. Earlier posterior and frozen-VST identity panels still pass.
Both R CI lanes include the new complete-fit identity check.

## Remaining adoption work

The full optional-method contract in the
[posterior inventory](posterior-primitives-2026-11-10.md#remaining-public-method-contract)
remains required. In particular, this kernel has not implemented all upstream
optimizer/control choices, low-rank SVD approximation, automatic prior grids,
mode/range handling, alternate component/likelihood families, pruning, missing
rows or warning semantics. It does not implement apeglm adaptation, MAP fitting
or posterior integration. Production adoption must retain the existing supported
methods/options and independent DESeq2 bounds alongside direct interface identity.

No public fitting/shrinkage path changes here, so unchanged GTEx/R workflows are
not retimed. There is no new speed or memory claim and no counted optimization
round. Final combined performance, memory, result authority and release artifacts
remain open after the remaining production ports.

## Reproduction

```sh
Rscript scripts/references/generate_mixture_fit_reference.R /tmp/mixture-fit.json
cmp /tmp/mixture-fit.json crates/rsdeseq2/tests/data/mixture_fit_r462.json
cargo test --locked -p rsdeseq2 --test mixture_fit
cargo build --release --locked -p rsdeseq2 --example evaluate_r_interface
Rscript scripts/validation/check_mixture_fit_identity.R \
  target/release/examples/evaluate_r_interface R_CHECK_LIBRARY /tmp/mixture-fit-identity
```

Use the pinned native R environment for regeneration and a checked candidate
library for the identity comparison. Keep the output directory fresh.
