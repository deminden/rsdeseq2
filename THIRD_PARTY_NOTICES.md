# Third-party notices

## Published coefficients

The `rcompat-numeric` crate incorporates no third-party source code. Two
files use coefficients printed in published papers:

- `src/numeric/special.rs`: the normal quantile uses the coefficients of
  Wichura, M. J. (1988), "Algorithm AS 241: The percentage points of the
  normal distribution", Applied Statistics 37(3):477-484.
- `src/probability/normal_log.rs`: the normal tail uses the rational
  approximations of Cody, W. J. (1969), "Rational Chebyshev approximations
  for the error function", Mathematics of Computation 23(107):631-637.

## Public-domain material

The `rcompat-shrinkage` crate incorporates no third-party source code under
a licence. Two files use public-domain numerical tables from netlib:

- `src/shrinkage/gamma_normalization.rs`: the Chebyshev series of
  W. Fullerton's SLATEC routine DGAMMA (<https://netlib.org/slatec/fnlib/dgamma.f>).
  The SLATEC library is in the public domain (SLATEC guide, section 4,
  <https://netlib.org/slatec/guide>).
- `src/shrinkage/log_f_quadrature.rs`: the 15-point Gauss-Kronrod nodes and
  weights of QUADPACK's DQK15I (<https://netlib.org/quadpack/dqk15i.f>),
  which is in the public domain.
