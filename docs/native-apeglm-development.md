# Local native apeglm development

This unreleased checkpoint moves the default `lfcShrink(type="apeglm")`
negative-binomial Cauchy/Laplace workflow into the standalone Rust engine.
The R wrapper imports arguments, aligns rows and constructs `DESeqResults`.
Likelihoods, empirical prior adaptation, optimization, Hessian/SD calculation,
normal probabilities, scale conversion and s-values run in Rust. No R or Rcpp
statistical routine is called. The default method no longer requires the apeglm
package at runtime; original apeglm remains the independent comparison engine.

The full [Rust-only migration](rust-only-wrapper-migration.md) remains open.
Current native public methods are `nbinomCR` and MAP-only `nbinomC`.
General/restricted custom likelihoods, `nbinomR`, random-start methods,
beta-binomial methods, grid intervals, parameter uncertainty, contrasts and
alternate optimizer policies still require implementation. Unfinished methods
fail explicitly, with no R calculation fallback. This checkpoint is not a new
release, and GitHub publication remains unauthorized.

## Engine and numerical contracts

`fit_negative_binomial_posterior` accepts borrowed counts, design, dispersions,
optional weights and offsets. Offsets may be log values, positive normalization
factors or sample size factors. Size-factor logs are calculated once per batch;
the wrapper does not allocate a repeated offset matrix or a default all-one
weight matrix. Only one row's likelihood/curvature work is retained at a time.
Process-level memory and end-to-end performance are reserved for the final
combined audit; this allocation change makes no measured speed/RSS claim.

Cauchy priors apply to shrunken coefficients, with independent wide Gaussian
priors on nuisance coefficients. MLE adaptation uses a bracketed empirical
variance root and caps the resulting scale at one. The native bracketed root
retains the reference absolute tolerance, including lower-bound and missing-MLE
behavior. Log2 MLEs, thresholds and result summaries are converted in the core.
All-zero negative-binomial rows retain missing modes, SDs and sign rates.

The initial mathematically well-converged L-BFGS-B trial failed the existing
public MAP bound on one of 1,200 genes: maximum log2 error `0.010455605` exceeded
`0.01`. Its gradient was closer to the optimum, but original apeglm intentionally
uses a different iteration and stopping policy. The accepted default prefit uses
a pure Rust adaptation of the MIT-licensed LBFGSpp unbounded algorithm and its
Nocedal-Wright strong-Wolfe search, with the reference normalization, additive
objective constant and stopping settings. Existing `rcompat-lbfgsb` supplies the
fallback refinement. No GPL apeglm implementation is copied into the crate.
The [third-party notices](../THIRD_PARTY_NOTICES.md) retain LBFGSpp attribution
and its MIT license.

Laplace curvature uses central differences of the analytic gradient, retaining
the reference step and symmetric Hessian. Native convergence counts are explicit
diagnostics; they are not manufactured to match another optimizer's counts.

## Focused evidence

The original reference uses `deseq2-r462-native`: R 4.6.2, DESeq2 1.52.0 and
apeglm 1.34.0. Thirteen independent posterior profiles cover weights, offsets,
two/three coefficients, several shrunken coefficients, narrow/wide priors,
empirical adaptation, signed thresholds, MAP-only output and zero-count rows.
The thirteenth profile is a frozen sensitive row that detects the rejected
optimizer's premature compatibility claim. Six prior-only profiles cover the lower boundary, interior
root, near-boundary root, scale cap, missing MLEs and the upper endpoint.

The existing public posterior bound remains `0.01` log2 units. The core table
uses the same bound converted to natural log units; sign/threshold s-value gates
remain `0.02`/`0.05`. Prior variances/scales use separate tight reference checks.
Original reference values and bounds are not changed to accept a candidate.

The larger 1,200-gene default/sign/threshold workflow compares native public
results directly with original DESeq2. A positive-control firewall blocks R
statistical routines, including apeglm fitting/prior/s-value routines, during
candidate calls. Direct R/native versus a separate R-free process compares
posterior-column, adapted-prior, nuisance-index and convergence bits for the
thirteen profiles plus a normalization-factor conversion profile. A separate full-dataset
comparison also covers the public size-factor and log2-input conversion path.
R dataframe/S4 ordering and other object presentation remain permitted.

The [machine-readable local report](data/native-apeglm-local.json) binds the exact source, executables and
measured errors. Full suites, sanitizer/platform checks, combined speed/memory
acceptance, final documentation reconciliation and source-bundle rebuilding
remain deferred until the remaining methods are implemented.

## Reproduction

```sh
Rscript scripts/generate_cauchy_posterior_reference.R "$FRESH_REFERENCE_JSON"
cargo test -p rsdeseq2 --test heavy_tail
cargo build --release --example evaluate_r_interface
Rscript scripts/check_native_apeglm_workflow.R "$DRIVER" "$R_LIBRARY" "$FRESH_OUTPUT"
```

The reference generator loads only original apeglm. It includes frozen input
bits for the sensitive optimizer row; candidate numerical outputs do not enter
reference calculations. An optional fourth checker argument reuses a saved
paired count-fit snapshot, avoiding repeated DESeq fitting while validating a
changed shrinkage implementation. Fresh output paths protect earlier evidence.
