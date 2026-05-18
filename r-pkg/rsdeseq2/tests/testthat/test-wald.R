testthat::test_that("R wrapper Wald test is explicit future work", {
  testthat::expect_error(nbinomWaldTestRust())
})
