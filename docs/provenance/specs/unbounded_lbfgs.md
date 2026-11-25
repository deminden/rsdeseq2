# Specification: unbounded L-BFGS with a strong-Wolfe line search

Status: written from allowed inputs only: Nocedal, J. and Wright, S. J.
(2006), Numerical Optimization, 2nd ed., Springer (chapter 3, line search
methods, Algorithms 3.5 and 3.6; chapter 7, L-BFGS, Algorithms 7.4 and
7.5), Liu, D. C. and Nocedal, J. (1989), "On the limited memory BFGS method
for large scale optimization", Math. Programming 45:503-528, the
documentation of R's `apeglm` and `RcppNumerical` packages for behaviour,
and black-box runs of R. No LBFGSpp, RcppNumerical, apeglm or R source was
consulted.

## Purpose

The heavy-tailed shrinkage engine (`crates/rcompat-shrinkage/src/shrinkage/heavy_tail/`)
finds each gene's posterior mode with an unbounded quasi-Newton method, as
apeglm does for its negative-binomial prior fits. Two callers:

- `fit.rs`: `unbounded_lbfgs::minimize(f, start)`;
- `beta_binomial.rs`: `unbounded_lbfgs::minimize_with_tolerances(f, start, 1e-12, 1e-10)`.

```rust
pub(super) fn minimize(f: impl FnMut(&[f64]) -> (f64, Vec<f64>), start: Vec<f64>) -> OptimResult;
pub(super) fn minimize_with_tolerances(
    f: impl FnMut(&[f64]) -> (f64, Vec<f64>),
    start: Vec<f64>,
    function_tolerance: f64,
    gradient_tolerance: f64,
) -> OptimResult;
```

`f` returns the objective and its gradient. `OptimResult` and `OptimCounts`
come from `rcompat_lbfgsb`; fill `par`, `value`, `counts` (function and
gradient evaluations) and `convergence` the way the callers read them
(read `fit.rs`, `beta_binomial.rs` and `analysis.rs` for which fields are
used and how `convergence` is interpreted; 0 means converged).

## Behaviour

`minimize` uses memory m = 6 correction pairs, at most 300 iterations,
gradient tolerance 1e-8 and function tolerance 1e-8 (documented settings of
the apeglm negative-binomial prefit).

Iteration (Nocedal and Wright, Algorithm 7.5 with the two-loop recursion of
Algorithm 7.4):

1. Evaluate f and g at the start. The first direction is -g; the first trial
   step is 1/||g|| (so the first move has unit length); later trial steps are 1.
2. Line search along the direction d satisfying the strong Wolfe conditions
   with c1 = 1e-4 and c2 = 0.9 (Algorithms 3.5 and 3.6: bracketing phase,
   then a zoom phase that picks the trial step inside the bracket by
   interpolation, with safeguards that keep it away from the bracket ends).
   Stop the search with an error status after a bounded number of trials
   (20 is a reasonable choice) or if the step becomes non-finite.
3. Store s = x_new - x and y = g_new - g when y's = s'y > 0 (skip the
   update otherwise), keeping the newest m pairs; scale the initial inverse
   Hessian by gamma = s'y / y'y of the newest pair.
4. Convergence (checked after each accepted step):
   - gradient: ||g|| <= gradient_tolerance * max(1, ||x||);
   - function: relative decrease over one past iteration,
     |f_prev - f| / max(|f|, 1) <= function_tolerance, unless black-box
     runs show another form;
   - iteration limit: report non-convergence.
5. Non-finite objective or gradient values at a trial point are treated as a
   failed trial (shrink the step).

Details that the books leave open (the zoom interpolation formula, the
exact convergence expressions, what is returned after a failed line
search) are to be chosen so that the existing apeglm-style reference tests
pass; they may be pinned down with black-box runs of `apeglm` in R
(`/Users/dendem/miniforge3/envs/deseq2-r462-native/bin/Rscript`) on inputs
you choose, recording any probe script under `scripts/provenance/`.

## Acceptance

- Unit tests: Rosenbrock (from (-1.2, 1)), a strictly convex quadratic in
  10 dimensions, and a logistic-regression likelihood converge to the
  known optimum at the stated tolerances; non-finite objectives are
  rejected cleanly.
- Every existing workspace test passes unchanged, in particular the
  heavy-tailed and beta-binomial posterior references
  (`crates/rsdeseq2/tests/data/*posterior*_r462.json`,
  `regression_*`, `random_start_posterior_r462.json`) under the default and
  `portable` profiles.
- Golden manifests may change; report which fields change.
