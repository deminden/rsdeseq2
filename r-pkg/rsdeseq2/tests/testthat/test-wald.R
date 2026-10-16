testthat::test_that("R wrapper Wald test validates primitive inputs", {
  testthat::expect_error(nbinomWaldTestRust())
  testthat::expect_error(nbinomWaldTestRust(baseMean = 1, beta = matrix(1, nrow = 2)))
})
