# Third-party notices

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
