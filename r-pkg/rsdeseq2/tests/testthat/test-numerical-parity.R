testthat::test_that("posterior parity checks cannot hide a tail error behind a large gene", {
  .rsdeseq2_expect_shrinkage_parity(c(1000000, 0), c(1000000, 0), 0.01)
  testthat::expect_failure(
    .rsdeseq2_expect_shrinkage_parity(c(1000000, 0.02), c(1000000, 0), 0.01)
  )
  testthat::expect_failure(
    .rsdeseq2_expect_shrinkage_parity(c(1, NA_real_), c(1, 0), 0.01)
  )
  testthat::expect_failure(
    .rsdeseq2_expect_shrinkage_parity(c(1, Inf), c(1, NA_real_), 0.01)
  )
})
