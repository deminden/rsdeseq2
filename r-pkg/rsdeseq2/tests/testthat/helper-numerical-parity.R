# Method-specific compatibility bounds for ashr/apeglm are deliberately retained.
# Check every gene's absolute error: vector-wide scaling can hide a tail error.
.rsdeseq2_expect_shrinkage_parity <- function(actual, expected, tolerance) {
  same_length <- length(actual) == length(expected)
  same_na <- identical(unname(is.na(actual)), unname(is.na(expected)))
  no_infinity <- !any(is.infinite(actual)) && !any(is.infinite(expected))
  finite <- if (same_length) is.finite(actual) & is.finite(expected) else logical()
  max_error <- if (any(finite)) max(abs(actual[finite] - expected[finite])) else Inf
  testthat::expect_true(
    same_length && same_na && no_infinity && max_error <= tolerance,
    info = sprintf(
      "same length=%s; same NA mask=%s; no infinity=%s; max abs error=%.8g; bound=%.8g",
      same_length, same_na, no_infinity, max_error, tolerance
    )
  )
}
