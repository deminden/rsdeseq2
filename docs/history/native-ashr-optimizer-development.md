# Local native ashr optimizer development

This unreleased checkpoint adds native `mixEM`, `w_mixEM`, `mixVBEM` and
`cxxMixSquarem` fitting to the shared [location engine](shrinkage-component-development.md).
The accepted default SQP path stays separate. R imports controls and presents
results; none of these methods executes an R or C++ numerical routine.

## Algorithms and controls

`MixtureFitOptions::optimizer` selects SQP or one of four EM policies. The
independent squared-extrapolation implementation follows the published
[Varadhan/Roland algorithm](https://doi.org/10.1111/j.1467-9469.2007.00585.x)
and [Du/Varadhan description](https://arxiv.org/abs/1810.11163). No GPL ashr or
SQUAREM source implementation is copied, translated or linked into Rust.
Responsibilities, weighted sufficient statistics, Dirichlet updates, variational
expected log weights, ELBO/KL evaluation and safeguards run in the core. No new
runtime dependency is added.

Generated non-SQP grids use the original concentrated initialization before
asymmetric component expansion and range constraints. Supplied grids keep their
stored starting weights; zero starts receive the original positive floor.
Variational fitting starts with unit concentrations. `prior="unit"` and positive
numeric Dirichlet parameters below one are accepted only for variational fitting.
Convex SQP/EM validation remains strict. A positive prior mass can fit an entirely
zero observation-weight batch. All-zero likelihood columns are removed before
fitting and their output weights restored as zeros.

`SquaremOptions` exposes tolerance, evaluation budget, the three step formulas,
step bounds/growth and the objective-increase safeguard. The default tolerance is
`min(1e-7,0.1/n)` after observation exclusion. Zero budgets return the starting
weights and a nonconvergence diagnostic without a fixed-point call. The K=1 R
policy honors the step-formula choice; the compiled-style policy uses formula 3
regardless of that argument, as the original does. Both original K=1 methods
ignore `square`; the native interface retains that control behavior. Trace
controls are validated but native progress text is not emitted. The later [advanced optimizer checkpoint](native-ashr-advanced-optimizer-development.md)
adds higher-order cycling and native interior-point fitting.

The compiled-style fixed-point map retains extrapolated signed weights, unlike
R EM's positive-clipped map. Rejecting or clipping these candidates changes its
safeguards and stopping path. Its budget rule also reserves one probe. The
variational map continues digamma through noninteger negative extrapolations
before the stabilizing update restores positive concentrations; rejecting them
changes sparse-prior trajectories. Positive digamma uses the existing native
Bernoulli recurrence factored out of the beta-binomial score. Accepted beta-
binomial/random-start outputs remain exact after that reuse.

## Explicit corrections to original compiled diagnostics

Original ashr 2.2-63 has three independently demonstrated compiled-method issues:

* If the second fixed-point probe is stationary, it reports convergence but
  returns the old iterate. A four-observation categorical problem returns
  `(0.5,0.5)` instead of its exact Dirichlet-penalized solution `(4/7,3/7)`.
  Rust commits the stationary second update. The fixture retains the original
  output and records the independently computed correction separately.
* Its reported objective can refer to another iterate. Rust's compiled-style
  diagnostic evaluates the returned state afresh. Original `B`, a fresh original-
  state evaluation and the corrected-state reference remain separate fixture
  fields. No source value or numerical tolerance is changed to close this gate.
* Its empty diagnostic-vector read emits `subscript out of bounds (index 0 >=
  vector size 0)`. Safe native code does not reproduce that access or warning.
  Other fitting warnings and convergence decisions remain covered.

For the source's second-probe cases, the generator applies two original EM
updates and tests the original stopping criterion to record corrected weights.
An independent categorical solution protects against duplicating the same bug.
The first-probe policy still returns the prior iterate within its requested
update tolerance. These explicit corrections are separate from the unchanged
posterior precision gates and the final complete warning/authority audit.

## Evidence and reproduction

Original-only generation uses R 4.6.2, ashr 2.2-63 and SQUAREM 2026.1 in
`deseq2-r462-native`. Ninety-six complete 48-observation profiles cover all four
methods, generated/supplied and zero-start priors, asymmetric families, alpha,
weighted/prior-only fitting, unit/fractional priors, step formulas/bounds,
nonmonotone/unchecked objectives, early convergence and small budgets. Ninety-
seven captured primitive optimizer problems include the categorical correction.
Inputs are captured before the original optimizer; candidate values never enter
the generator. Original weights, summaries, convergence, iteration counts and
objective values remain frozen alongside separate correction references.

All weights, concentrations, posterior summaries, sign probabilities and log
likelihoods retain the existing `2e-6` batch bound. Iteration counts and
convergence decisions match. Independent categorical/Dirichlet solutions and
invalid control/likelihood/optimizer cases are also checked. Fifteen focused Rust
integration tests, core/backend Clippy and vendor integrity pass.

Actual R/native and R-free processes agree on every output bit for all 96 full
batches and all 97 optimizer/diagnostic batches. Candidate requests block R
probability, Gamma, optimizer and original ashr/SQUAREM routines. Four public
256-gene workflows compare with original ashr: maximum posterior mean/SD/s-value
errors are respectively `3.70e-7`, `1.20e-7` and `1.84e-8`, below the unchanged
public `1e-3`/`0.02` limits. Exact accepted-source replays cover default location,
adaptive queries, NB/general/beta-binomial/random-start and posterior-analysis
payload sets against the immutable `031e75d` driver.

[Machine-readable evidence](../data/native-ashr-optimizer-local.json) binds the actual
source and local artifacts. This checkpoint makes no throughput/RSS claim and
counts no optimization round. Full suites, combined precision/speed/memory,
sanitizers, platform checks and the source bundle remain final combined gates.
No GitHub publication is authorized.

```sh
Rscript scripts/references/generate_ashr_optimizer_reference.R "$FRESH_REFERENCE_JSON"
cargo test -p rsdeseq2 --test ashr_optimizers --test location_shrinkage --test mixture_fit --test source_convergence --test mixture_grid --test regression_posterior
cargo build --release --example evaluate_r_interface
Rscript scripts/validation/check_location_shrinkage_identity.R "$DRIVER" "$R_LIBRARY" "$FRESH_BATCH_OUT" "$FRESH_REFERENCE_JSON"
Rscript scripts/validation/check_native_ashr_optimizers.R "$DRIVER" "$R_LIBRARY" "$FRESH_OPTIMIZER_OUT" "$PAIRED_FITS_RDS"
```

The optional fourth optimizer-checker argument reuses an accepted count-fit
snapshot for the public comparisons. Fresh output paths preserve previous
artifacts. The later [advanced checkpoint](native-ashr-advanced-optimizer-development.md)
adds interior-point/higher-order fitting. Still open: remaining SQP controls, estimated modes, other likelihood families, output subsets and the
complete warning/error/authority contract. C19/C21/C23 are not closed here.
