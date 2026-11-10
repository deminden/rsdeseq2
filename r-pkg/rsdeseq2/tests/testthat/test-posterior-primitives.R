testthat::test_that('posterior bridge preserves point masses and validates complete requests', {
  request <- list(operation='posterior_normal_mixture',weights=c(.2,.8),means=c(0,0),
    sds=c(0,1),estimates=c(-2,0,2),standard_errors=c(.5,.5,.5),thresholds=c(-Inf,0,Inf))
  answer <- .rsdeseq2_formula_request(request)
  # Independent conjugate single-observation calculation, including the zero atom.
  null <- .2 * stats::dnorm(request$estimates,sd=.5)
  alt <- .8 * stats::dnorm(request$estimates,sd=sqrt(1.25))
  w <- alt/(null+alt)
  mu <- request$estimates/1.25
  variance <- .25/1.25
  testthat::expect_equal(answer$posterior_mean,w*mu,tolerance=2e-12)
  testthat::expect_equal(answer$posterior_sd,sqrt(w*(variance+mu^2)-(w*mu)^2),tolerance=2e-12)
  testthat::expect_equal(answer$cdf[c(2,5,8)],(1-w)+w*stats::pnorm(0,mu,sqrt(variance)),tolerance=2e-12)
  testthat::expect_equal(answer$upper_tail[c(2,5,8)],w*stats::pnorm(0,mu,sqrt(variance),lower.tail=FALSE),tolerance=2e-12)
  testthat::expect_equal(.rsdeseq2_formula_request(list(operation='posterior_svalues',rates=c(.2,NA,.1,.2))),c(.15,NA,.1,1/6),tolerance=2e-15)
  bad <- list(list(weights=c(0,0)),list(means=0),list(sds=c(0,-1)),
    list(standard_errors=.5),list(estimates=c(-2,NA,2)),list(standard_errors=c(.5,0,.5)),
    list(thresholds=NA_real_),list(weights='invalid'))
  for (change in bad) testthat::expect_error(.rsdeseq2_formula_request(utils::modifyList(request,change)))
  for (rate in list(-.1,1.1,Inf))
    testthat::expect_error(.rsdeseq2_formula_request(list(operation='posterior_svalues',rates=rate)))
})
