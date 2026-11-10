# Local beta-binomial shrinkage development

This unreleased checkpoint adds standalone beta-binomial logit shrinkage through
`betabinR`, `betabinCR` and MAP-only `betabinC`. It extends the
[general-regression engine](native-regression-development.md) and uses its shared
prior adaptation, validation, posterior summaries, scale conversion and s-values.
R imports inputs and presents results; likelihoods, gradients, optimization and
uncertainty calculations run in the same Rust engine as the standalone interface.

## Model and fitting policies

Select `RegressionLikelihood::BetaBinomial` and its matching `RegressionMethod`
in `fit_regression_posterior`. The parameter matrix has one row per response
and contains concentration followed by one trial-count column per sample.
Concentrations must be positive and finite; successes must lie between zero
and finite trial counts. Empty models and malformed matrices fail before fitting.
The R interface accepts the same matrix through `param` with the corresponding
`apeMethod` in `lfcShrink(type="apeglm")`.

The model uses a logit link, fits zero-success and zero-trial rows, and permits
intercept shrinkage. Supplied offsets are ignored after shape validation, as in
original apeglm. Gaussian nuisance and Cauchy shrinkage priors retain their
existing scale/adaptation semantics. `betabinR` starts at zero and uses the
analytic score with BFGS or bounded L-BFGS-B. `betabinCR` uses a fast native
prefit and computes uncapped analytic-score curvature; failed prefits or invalid
variance estimates trigger analytic refitting. MAP-only `betabinC` omits SD and
sign/threshold summaries.

The fast prefit retains the original probability caps, normalization anchor,
additive constant and distinct stopping settings. The existing MIT-licensed
LBFGSpp adaptation now accepts those settings without changing the default
negative-binomial settings. No GPL apeglm or R optimizer implementation is copied.
Gamma and digamma increments follow the recurrence and asymptotic equations in
[DLMF 5.5](https://dlmf.nist.gov/5.5) and [DLMF 5.11](https://dlmf.nist.gov/5.11),
with stable differences to avoid subtracting large, nearly equal values.

Saturated profiles exposed three missing common BFGS policies: periodic metric
resets, a fresh-metric convergence probe retaining the previous score, and an
iteration budget counted by accepted gradients including the initial gradient.
Independent original-R evaluation traces and frozen posterior profiles identify
these policies. Simply changing the relative-stop condition failed. The accepted
correction retains all previous general/negative-binomial references and bounds.

## Focused evidence

Original references use R 4.6.2, DESeq2 1.52.0 and apeglm 1.34.0 in
`deseq2-r462-native`. Thirty-nine profiles cover 156 posterior rows across all
three beta-binomial methods, including saturated responses, weights, small/large
concentrations, ignored missing offsets, zero trials, three coefficients, adapted
and intercept priors, all-Gaussian/empty nuisance sets, bounds and zero weights.
Nine independent likelihood/analytic-score profiles use a tight
`1e-9 + 1e-12 * abs(reference)` gate. Invalid concentration, success/trial, matrix
and empty-model inputs are rejected explicitly.

All posterior profiles pass the existing mode/SD gate of `0.01` log2 units, sign
quantity gate of `0.02` and threshold gate of `0.05`; original references and
bounds remain unchanged. Actual R/native and separate R-free processes agree
bit for bit on selected posterior columns, s-values, adapted priors, nuisance
indices, convergence status and warning strings for all 39 batches. R matrix,
gamma, density, optimizer and posterior calculations are blocked during candidate
calls, with a positive control proving the block is active.

Public checks on 256 genes reuse the accepted count fit and compare the three
new methods directly with original apeglm on identical posterior inputs. Maximum
log2 mode differences are about `4.74e-7` for `betabinR`, `6.63e-7` for `betabinCR`
and `1.28e-6` for `betabinC`. Existing 1,200-gene general, `nbinomR` and custom
likelihood checks also pass. The corrected BFGS policy lowers their maximum
mode differences to `2.16e-6`, `3.82e-7` and `2.16e-6`, respectively. The default
negative-binomial thirteen-profile, normalization-factor and public workflow
checks still pass. Exact errors and source/artifact identities are in the
[local report](data/native-beta-binomial-local.json).

Native optimizer warnings are transported as data and emitted by the R wrapper.
All 39 frozen profiles match the original ordered warning-event categories;
37 match the complete wording. Two profiles classify a floating-point plateau
as insufficient decrease rather than insufficient precision. In the larger
public fast-method comparison, Rust emits 21 line-search warnings while original
R emits 19; numerical differences remain below `1.3e-6` log2 units. These are
actual native solver outcomes, not warnings synthesized to reproduce source
counts. Native/wrapped warning strings agree exactly. Full warning-policy
reconciliation remains visible in the final authority audit.

## Reproduction and remaining scope

```sh
Rscript scripts/generate_beta_binomial_posterior_reference.R "$FRESH_REFERENCE"
cargo test -p rsdeseq2 --test regression_posterior --test heavy_tail
cargo test -p rsdeseq2 --lib analytic_likelihood_and_score_match_independent_original_apeglm
cargo build --release --example evaluate_r_interface
Rscript scripts/check_native_regression_workflow.R "$DRIVER" "$R_LIBRARY" "$FRESH_OUTPUT" "$PAIRED_FITS_RDS" "$FRESH_REFERENCE"
```

The generator loads only original apeglm. Fresh output paths protect earlier
evidence. The shared checker reuses saved original posterior comparisons and
accepted count fits when supplied, avoiding repeated unchanged fitting.

The later [random-start checkpoint](native-random-start-development.md) adds
`betabinC*` and `nbinomC*` with portable native RNG state. Grid intervals,
parameter uncertainty, contrasts and other
optimizer policies remain unfinished, alongside advanced ashr methods/options,
the final authority audit and combined precision/speed/memory/sanitizer/platform/
source-bundle gates. No unchanged quantification timing or RSS is remeasured for
this checkpoint; no speed gain or optimization round is claimed. All work stays
local and this is not a release-completion claim.
