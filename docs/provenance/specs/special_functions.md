# Specification: special functions and the gamma distribution

Status: written from allowed inputs only: published papers and the NIST
Digital Library of Mathematical Functions (DLMF), R help pages for
behaviour, and black-box runs of R 4.6.2. No statrs, Boost, R or Cephes
source was consulted.

## Purpose

`crates/rcompat-numeric/src/numeric/special.rs` provides the special
functions the pipeline uses: log-gamma in the negative-binomial likelihood,
digamma in the dispersion score, regularized incomplete gamma for
chi-square tails, beta functions in the shrinkage likelihoods, the normal
density and quantile, and a gamma distribution with cdf, survival function
and quantile. The earlier file was a port of a third-party crate; this one
must be written independently.

## Public API (unchanged signatures)

```rust
pub fn ln_gamma(x: f64) -> f64;
pub fn gamma(x: f64) -> f64;
pub fn gamma_ur(a: f64, x: f64) -> f64;
pub fn gamma_lr(a: f64, x: f64) -> f64;
pub fn digamma(x: f64) -> f64;
pub fn ln_beta(a: f64, b: f64) -> f64;
pub fn beta(a: f64, b: f64) -> f64;
pub fn standard_normal_pdf(x: f64) -> f64;
pub fn standard_normal_inverse_cdf(p: f64) -> f64;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GammaDistribution { shape: f64, rate: f64 }
impl GammaDistribution {
    pub fn new(shape: f64, rate: f64) -> Option<Self>;
    pub fn chi_squared(freedom: f64) -> Option<Self>;
    pub fn freedom(&self) -> f64;
    pub fn cdf(&self, x: f64) -> f64;
    pub fn sf(&self, x: f64) -> f64;
    pub fn inverse_cdf(&self, p: f64) -> f64;
}
```

## Behaviour

Domain and special values follow R's `lgamma`, `gamma`, `digamma`,
`pgamma`, `beta`, `lbeta`, `dnorm`, `qnorm` and `qgamma` (their help pages,
and R runs for edge values), with these contract points that callers rely on:

- `ln_gamma(x)` is log |Gamma(x)|; +Inf at 0 and negative integers.
- `gamma(x)`: NaN at 0 and negative integers; overflow gives +Inf.
- `gamma_lr(a, x)` = P(a, x), `gamma_ur(a, x)` = Q(a, x) = 1 - P(a, x),
  each computed directly (not as 1 minus the other) so that small tails keep
  relative accuracy down to about 1e-300. NaN arguments give NaN.
  They panic (`assert!`) unless a > 0 and finite and x >= 0 and finite;
  P(a, 0) = 0, Q(a, 0) = 1.
- `digamma(x)`: psi(x); NaN at 0 and negative integers.
- `ln_beta(a, b)`, `beta(a, b)`: panic unless a > 0 and b > 0. `ln_beta`
  must stay accurate when one argument is very large (avoid subtracting
  large log-gammas; DLMF 5.11 asymptotics or the ratio form).
- `standard_normal_pdf(x)` = exp(-x^2/2) / sqrt(2 pi).
- `standard_normal_inverse_cdf(p)`: panics outside [0, 1]; -Inf at 0,
  +Inf at 1. Use Wichura, M. J. (1988), "Algorithm AS 241: The percentage
  points of the normal distribution", Applied Statistics 37(3):477-484,
  coefficients as printed in the paper.
- `GammaDistribution::new(shape, rate)`: None unless shape and rate are
  finite and > 0. `chi_squared(k)` = new(k/2, 1/2); `freedom()` = 2*shape.
  `cdf(x)` = P(shape, rate*x) and `sf(x)` = Q(shape, rate*x) for x > 0;
  cdf = 0 and sf = 1 for x <= 0; cdf = 1 and sf = 0 for x = +Inf.
  `inverse_cdf(p)` panics outside [0, 1]; returns 0 at p = 0 and +Inf at
  p = 1; otherwise a safeguarded Newton or Halley iteration on the cdf from
  the Wilson-Hilferty starting value (Wilson and Hilferty, 1931), bracketed
  so it always converges, to full binary64 accuracy.

All elementary operations go through `crate::numeric::math::ProfileMath`
(`ln_p`, `exp_p`, `powf_p`, ...), so the `portable` profile gives identical
bits on every host and the R profiles use the host C library, as the rest
of the crate does. Do not call `f64::ln`/`exp` directly.

## Suggested methods (all published)

- log-gamma and gamma: a Lanczos approximation (Lanczos, 1964, SIAM J.
  Numer. Anal. B 1:86-96) with coefficients you derive yourself in high
  precision, or Stirling's series (DLMF 5.11.1) with upward recurrence
  (DLMF 5.5.1) and a separate treatment near 1 and 2 (for example the
  Taylor series of log-gamma about 1 and 2, DLMF 5.7.3) so that the result
  keeps relative accuracy near its zeros; reflection (DLMF 5.5.3) for
  x < 0.5.
- digamma: recurrence (DLMF 5.5.2) up to x >= 6 or more, then the
  asymptotic expansion (DLMF 5.11.2); reflection (DLMF 5.5.4) for x < 0.
  The crate already has `probability::digamma_positive` written from DLMF;
  reuse it where it applies.
- incomplete gamma: the power series (DLMF 8.7.1) for x < a + 1 and the
  continued fraction (DLMF 8.9.2) otherwise, with the prefactor
  x^a e^-x / Gamma(a+1) evaluated in logs, and Temme's uniform expansion
  (DLMF 8.12, Temme 1979/1996) or an equivalent for large a near x = a.
- ln_beta: log-gammas for moderate arguments; for a large argument use
  DLMF 5.11.12-style asymptotics for the ratio Gamma(b)/Gamma(a+b).

Any other published method is acceptable if it is cited in a comment.

## Acceptance

1. A reference generator `scripts/references/generate_special_functions_reference.R`
   (base R only, refuses to overwrite) writes
   `crates/rcompat-numeric/tests/data/special_functions_r462.tsv` with R
   values on grids that cover small, moderate and large arguments, both
   tails, negative non-integers, and the edge values above. Register it in
   `scripts/references/MANIFEST.tsv`.
2. Unit tests compare against it. Targets: relative error at most 4e-15
   for `ln_gamma` away from its zeros (absolute 1e-15 within 0.1 of 1 and
   2), `gamma`, `digamma` away from its zero near 1.4616, `beta`, `ln_beta`
   and the normal functions; at most 1e-13 relative for P and Q in the
   body and 1e-12 relative in tails down to 1e-300; 1e-12 relative for
   `inverse_cdf`. Report the largest observed error per function.
3. Every existing workspace test passes unchanged at its current
   tolerance, under the default (`r-macos` on the Mac) and `portable`
   profiles.
4. Golden manifests are expected to change (last digits). Do not rewrite
   them; report which fields change.
