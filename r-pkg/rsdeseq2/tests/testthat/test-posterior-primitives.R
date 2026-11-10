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

testthat::test_that('native mixture fitting retains weighted prior modes and diagnostics', {
  request<-list(operation='posterior_fit_mixture',observations=2L,components=2L,
    log_likelihoods=c(0,-Inf,-Inf,0),observation_weights=c(3,1),dirichlet=c(2,2),initial_weights=NULL,
    options=list(max_iterations=1000L,gradient_tolerance=1e-8,em_iterations=20L,active_iterations=100L,zero_threshold=1e-8,likelihood_offset=0))
  fit<-.rsdeseq2_formula_request(request)
  testthat::expect_true(fit$converged)
  testthat::expect_equal(fit$weights,c(2/3,1/3),tolerance=1e-12)
  testthat::expect_equal(fit$penalized_log_likelihood,4*log(2/3)+2*log(1/3),tolerance=1e-12)
  testthat::expect_lte(fit$dual_residual,1e-8)
  for(change in list(list(observation_weights=1),list(dirichlet=c(0,1)),list(components=3L),
      list(options=list(likelihood_offset=-1)),list(options=list(active_iterations=0L))))
    testthat::expect_error(.rsdeseq2_formula_request(utils::modifyList(request,change,keep.null=TRUE)))
})
