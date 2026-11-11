testthat::test_that('canonical count likelihoods import stored data and preserve numerical kernels', {
  testthat::skip_if_not_installed('ashr')
  y <- c(0,2,7)
  models <- list(ashr::lik_pois(y,c(.5,1,2)), ashr::lik_pois(y,c(.5,1,2),'log'),
    ashr::lik_binom(y,9), ashr::lik_binom(y,9,'logit'))
  for (lik in models) {
    payload <- .rsdeseq2_native_ash_likelihood(lik)
    testthat::expect_identical(payload$type,lik$name)
    x <- c(-1,0,1); a <- c(-2,-1,0); b <- c(-1,0,1)
    expected <- list(log_cdf=lik$lcdfFUN(x),log_density=lik$lpdfFUN(x),
      first=as.numeric(lik$etruncFUN(matrix(a,ncol=1),matrix(b,ncol=1))),
      second=as.numeric(lik$e2truncFUN(matrix(a,ncol=1),matrix(b,ncol=1))))
    actual <- .rsdeseq2_formula_request(list(operation='posterior_builtin_likelihood',likelihood=payload,
      arguments=x,lower=a,upper=b))
    for (field in names(expected)) testthat::expect_equal(actual[[field]],expected[[field]],tolerance=2e-8)
  }
})

testthat::test_that('count likelihood imports never run modified closures or deferred captures', {
  testthat::skip_if_not_installed('ashr')
  called <- FALSE
  lik <- ashr::lik_pois(c(0,2),link='log')
  lik$lpdfFUN <- function(x) { called <<- TRUE; stop('modified closure ran') }
  testthat::expect_error(.rsdeseq2_native_ash_likelihood(lik),'unsupported|one expression')
  testthat::expect_false(called)
  lik <- ashr::lik_pois(c(0,2),link='log'); env <- environment(lik$lpdfFUN)
  rm('y',envir=env); makeActiveBinding('y',function() { called <<- TRUE; c(0,2) },env)
  testthat::expect_error(.rsdeseq2_native_ash_likelihood(lik),'active')
  testthat::expect_false(called)
  lik <- ashr::lik_pois(c(0,2),link='log'); env <- environment(lik$lpdfFUN)
  delayedAssign('y1',{ called <<- TRUE; c(1e-5,2.00001) },assign.env=env)
  testthat::expect_error(.rsdeseq2_native_ash_likelihood(lik),'delayed|stored')
  testthat::expect_false(called)
  lik <- ashr::lik_binom(c(0,2),3,'logit'); environment(lik$lpdfFUN)$y1 <- c(1,2)
  payload <- .rsdeseq2_native_ash_likelihood(lik)
  testthat::expect_error(.rsdeseq2_formula_request(list(operation='posterior_builtin_likelihood',
    likelihood=payload,arguments=c(0,0),lower=c(-1,-1),upper=c(1,1))),'modified adjusted count')
  lik <- ashr::lik_pois(c(0,2)); lik$data$y <- c(1,2)
  testthat::expect_error(.rsdeseq2_native_ash_likelihood(lik),'metadata differs')
})
