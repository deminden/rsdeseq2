# Standalone posterior foundation

The Rust-only shrinkage migration now has fixed normal-mixture posterior and
ordered s-value primitives in `rsdeseq2::shrinkage`. The primitive R transport
calls these same kernels. This is a numerical foundation: public `lfcShrink`
still uses the established optional R packages for ashr/apeglm. Prior fitting,
other mixture/likelihood families, public options and production adoption remain
required before those methods satisfy the Rust-only contract.

## Implemented scope

`NormalMixture` validates finite component means, nonnegative SDs and nonnegative
weights, with at least one positive weight. Weights need not sum to one. It
conditions on a finite estimate and positive finite SE with representable
variance. Component probabilities normalize in log space; zero-weight components
retain their positions. The result exposes the mean, SD, component probabilities,
right-continuous CDF and directly evaluated upper tail. Point masses contribute
to the CDF at their location and to the upper tail strictly above the threshold.
Upper tails do not subtract a near-one CDF. Nonrepresentable likelihoods/moments
return errors. Missing-observation presentation belongs to the future full
method, rather than being silently invented by this finite-observation primitive.

`svalues` averages sorted local error rates and restores their original order.
Ties preserve input order, including signed zero. Missing entries remain missing
and do not enter the averages; finite rates outside `[0,1]` are errors. This
matches apeglm's ordering and the wrapper's thresholded ashr ordering. It does
not substitute for all of ashr's missing-value behavior or BH adjustment.

The implementation is original Rust using conjugate normal-mixture arithmetic
and the existing standalone normal-tail kernel. It does not load R, use R's math
runtime or call an R function. The existing public shrinkage routes and defaults
are unchanged, so no unchanged workflow is retimed and no new speed/memory or
optimization-round claim is made.

## Independent references and identity

The [generator](../../scripts/references/generate_normal_mixture_posterior_reference.R) uses
only original R 4.6.2, ashr 2.2-63 and apeglm 1.34.0. Local regeneration is
byte-identical. Eight priors cover a single normal, zero/shifted point masses,
centered/shifted mixtures, zero weights, a `1e-180` weight and widely differing
component scales. Each prior has thirteen estimates and five SEs: **520
observations**, with lower and directly evaluated upper tails at ten thresholds.
Eight s-value panels cover empty/single, sorted/repeated, missing/all-missing,
signed-zero and fractional rates. See the
[frozen fixture](../../crates/rsdeseq2/tests/data/normal_mixture_posterior_r462.json)
and [Rust tests](../../crates/rsdeseq2/tests/posterior_primitives.rs).

Moment checks use `2e-12 + 2e-12 * abs(reference)`; component probabilities and
both tails use `1e-300 + 2e-11 * abs(reference)`, retaining meaningful checks on
small probabilities. S-values use an absolute `2e-15` bound. These are new
primitive checks; existing public-method parity bounds and frozen references
are unchanged. Point-mass boundaries and malformed inputs have explicit tests.

The [direct interface check](../../scripts/validation/check_posterior_identity.R) requires
**every output bit** to match between the R-free executable and actual native
R transport, including all component weights and both tails. During bridge
calls it blocks R `sqrt`, `exp`, `log`, `mean`, `cumsum`, `sort`, `rank`, `pnorm`
and `dnorm`; a positive control verifies the guard. Both R CI lanes run this
check. It is a primitive identity check, not evidence that public `lfcShrink`
has already been migrated.

The [validation record](../data/posterior_primitives_validation_2026-11-10.json)
binds source/fixture/generator hashes and the checked artifacts. Local checks
pass: 1,080 Rust tests, three backend tests, 2,650 R assertions with source check
status OK, 47 Python tests with one platform skip, Clippy and MSRV 1.97.1 in
both graphs, and matching vendored source. All 326 archived R/C/Rust files match
the source used for checking. No unsafe/native buffer interface changed.

## Remaining public-method contract

The production wrappers forward more than default normal-mixture shrinkage.
The installed pinned packages and actual wrappers were inspected before this
foundation was written. Their public documentation is available from
[ashr](https://stephens999.r-universe.dev/ashr/doc/manual.html) and
[apeglm](https://bioconductor.org/packages/release/bioc/manuals/apeglm/man/apeglm.pdf).
The following is the migration inventory, not a new supported-options promise:

| Area | Existing behavior that the full port must address |
| --- | --- |
| ashr mixture/likelihood | Normal, uniform, half-uniform, one-sided uniform and half-normal families; supported normal/t or other likelihood combinations must retain the original accepted/rejected decisions. Normal mixtures with t likelihood already reject upstream. |
| ashr prior fitting | Automatic/supplied grids, point mass and prior choices, null weights, fixed/supplied mixtures, mode/range, alpha scaling, pruning threshold, observation weights, optimization method/control and convergence diagnostics. |
| ashr forwarded options | `...` can override estimates/SEs, mixture/likelihood and fitting/output arguments. The existing threshold calculation uses the original result estimates/SEs even if forwarded fitting arguments replace them; resolve this explicitly before changing behavior. |
| ashr missingness | Full result summaries use prior moments for excluded zero/infinite-SE or missing observations. Raw `cdf_post` handles missing and infinite SE differently. The installed package's q-value/s-value helper can recycle values and warn on nonfinite estimates; a replacement must audit these decisions instead of silently dropping them. |
| apeglm fit and adaptation | Coefficient selection, supplied counts/design/dispersion/offset/weights, adaptive or supplied prior controls, no-shrink coefficients, natural-log/log2 conversion, MAP search, posterior SD/integration, sign and threshold probabilities. |
| apeglm method/options | `general`, negative-binomial and forwarded beta-binomial method choices; fast MAP-only and random-start variants; optimizer/bounds, interval/integration controls, multiplier, link/sign policy and parameter uncertainty. Preserve actual upstream method restrictions. |
| callable arguments | Current forwarded likelihood objects/functions can contain executable R code. The final engine needs a standalone representation and validation policy; putting an R callback in the bridge does not satisfy the requirement. |
| output contract | Row/subset alignment, missingness, s-values/thresholds, independent-filter metadata, coefficient descriptions, convergence messages and prior information must accompany numerical parity. |

Next, implement prior fitting and the remaining likelihood/component kernels
behind independently frozen tests, then route the supported method matrix
through the shared engine. Keep independent DESeq2 precision checks separate
from Rust/wrapped bit identity. Only after production adoption should the
combined speed, memory, result-authority and final artifact gates close.

## Reproduction

```sh
Rscript scripts/references/generate_normal_mixture_posterior_reference.R /tmp/posterior.json
cmp /tmp/posterior.json crates/rsdeseq2/tests/data/normal_mixture_posterior_r462.json
cargo test --locked -p rsdeseq2 --test posterior_primitives
cargo build --release --locked -p rsdeseq2 --example evaluate_r_interface
Rscript scripts/validation/check_posterior_identity.R \
  target/release/examples/evaluate_r_interface R_CHECK_LIBRARY /tmp/posterior-identity
```

Use the pinned native R environment for regeneration and a checked candidate
library for the interface comparison. Keep the output directory fresh.
