# Specification: ashr-compatible policy of the mixture SQP solver

Status: written from allowed inputs only: Kim, Y., Carbonetto, P.,
Stephens, M. and Anitescu, M. (2020), "A fast algorithm for maximum
likelihood estimation of mixture proportions using sequential quadratic
programming", Journal of Computational and Graphical Statistics
29(2):261-273 (with its supplement); Nocedal and Wright (2006), Numerical
Optimization, chapter 16 (active-set methods for convex QP); the R help
pages of `mixsqp::mixsqp` and `ashr::ash` (controls and their defaults);
and black-box runs of R. No mixsqp, ashr, Armadillo or LAPACK source was
consulted.

## Purpose

`crates/rcompat-shrinkage/src/shrinkage/mixture_fit/sqp.rs` fits mixture
proportions by SQP. When `MixtureSqpOptions::ashr_compatibility` is true
(the setting the ashr-style interface uses), it calls three items of the
submodule `sqp/compatibility.rs`, which is to be rewritten:

```rust
impl Problem {
    pub(super) fn ashr_hessian_into(&self, weights: &[f64], hessian: &mut DenseMatrix);
    pub(super) fn ashr_expectation_step(&self, weights: &[f64]) -> Result<Vec<f64>, Error>;
}
pub(super) fn quadratic(
    h: &DenseMatrix,
    c: &[f64],
    start: &[f64],
    max_iterations: usize,
    options: MixtureSqpOptions,
) -> /* see the call site in sqp.rs */;
```

Read `sqp.rs`, `mixture_fit.rs` and `mixture_fit/controls.rs` (current code,
allowed) for the `Problem` fields (likelihood matrix, row weights, prior
counts, total mass), the objective and gradient conventions, and what the
caller does with the result of `quadratic`.

## Behaviour

- Objective (paper, section 2): minimise
  f(x) = -sum_i w_i log(sum_k L_ik x_k + eps) over the simplex, with eps the
  `eps` control of `mixsqp` (help page: a small non-negative number added
  inside the logarithms). The SQP works on the relaxed problem of section 3
  (x >= 0 with the sum constraint handled as in the paper).
- `ashr_hessian_into`: the Hessian of that objective at `weights`,
  sum_i w_i L_i L_i' / (L_i x + eps)^2, written into `hessian`.
- `ashr_expectation_step`: one EM update of the proportions (paper,
  section 2, the EM iteration used before SQP; `numiter.em` in the help
  page), x_k <- sum_i w_i L_ik x_k / (L_i x + eps) / sum_i w_i, combined with
  the prior counts the way `sqp.rs`'s own `expectation_step` does, and with
  proportions below a floor set to that floor. The floor value and exactly
  where it applies are to be fixed by black-box runs of `mixsqp` (for
  example by comparing proportions after `numiter.em` EM steps with
  `maxiter.sqp = 0`); record the probe under `scripts/provenance/`.
- `quadratic`: the active-set method for the convex QP subproblem
  min_p 1/2 p'Hp + c'p, p >= 0 (paper, Algorithm 2; Nocedal and Wright,
  Algorithm 16.3), starting from `start`, with at most `max_iterations`
  iterations, `options.active_tolerance` for deciding which coordinates are
  zero, `options.search_direction_threshold` for when a search direction is
  considered zero (`zero.threshold.searchdir` in the help page), and
  `options.ridge_multiplier` for regularising a singular reduced Hessian.
  Solving the reduced system must detect numerically singular directions
  (for example with a reciprocal condition estimate from Higham's 1-norm
  estimator, Higham 1988, ACM TOMS 14:381-396) and drop them, rather than
  return non-finite values.

The help page of `mixsqp` lists the controls and defaults
(`convtol.sqp`, `zero.threshold.solution`, `zero.threshold.searchdir`,
`suffdecr.linesearch`, `stepsizereduce`, `minstepsize`,
`identity.contrib.increase`, `eps`, `maxiter.sqp`, `maxiter.activeset`,
`numiter.em`); `controls.rs` already maps them.

Details that the paper and help pages leave open are to be chosen so that
the existing references pass, and may be fixed by black-box runs of
`mixsqp`/`ashr` in R (`/Users/dendem/miniforge3/envs/deseq2-r462-native/bin/Rscript`,
`verbose = TRUE` prints the per-iteration objective, step size and active
set size).

## Acceptance

- Every existing workspace test passes unchanged at its current tolerance,
  in particular `sqp_controls_r462.json`, `sparse_sqp_r462.json`,
  `mixture_fit_r462.json`, `mixture_convergence_r462.json`,
  `ashr_*_r462.json`, `adaptive_result_r462.json` and the R-wrapper
  location-shrinkage tests, under the default and `portable` profiles.
- Unit tests for the three items on small hand-checked problems.
- Golden manifests may change; report which fields change.
