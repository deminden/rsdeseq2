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
