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

testthat::test_that('selected adaptive outputs preserve conjugate values and omitted fields', {
  x <- c(-1,0,2); se <- rep(.5,3)
  prior <- structure(data.frame(pi=1,mean=0,sd=1),class='normalmix')
  run <- function(output,threshold=0) .rsdeseq2_native_ash(
    list(betahat=x,sebetahat=se,method='shrink',mixcompdist='normal',g=prior,
      fixg=TRUE,outputlevel=output),x,se,threshold,TRUE)
  mean <- .8*x; sd <- rep(sqrt(.2),3)
  for(level in c(1,4)) {
    answer <- run(level)
    testthat::expect_equal(answer$result$PosteriorMean,mean,tolerance=2e-12)
    testthat::expect_equal(answer$result$PosteriorSD,sd,tolerance=2e-12)
    testthat::expect_null(answer$result$svalue)
    testthat::expect_equal(answer$fitted_g$pi,1)
  }
  answer <- run(c('PosteriorMean','svalue'))
  testthat::expect_equal(answer$result$PosteriorMean,mean,tolerance=2e-12)
  testthat::expect_null(answer$result$PosteriorSD)
  testthat::expect_null(answer$fitted_g)
  rates <- stats::pnorm(-abs(mean)/sd)
  expected <- numeric(3);expected[order(rates)]<-cumsum(sort(rates))/seq_along(rates)
  testthat::expect_equal(answer$result$svalue,expected,tolerance=2e-12)
  answer <- run(5)
  testthat::expect_null(answer$result$PosteriorMean)
  testthat::expect_null(answer$result$PosteriorSD)
  testthat::expect_null(answer$fitted_g)
  rates <- ifelse(mean>0,stats::pnorm(.25,mean,sd),stats::pnorm(-.25,mean,sd,lower.tail=FALSE))
  expected[order(rates)]<-cumsum(sort(rates))/seq_along(rates)
  testthat::expect_equal(run(4,.25)$result$svalue,expected,tolerance=2e-12)
  testthat::expect_error(run(5,.25),'PosteriorMean and fitted_g')
})

testthat::test_that('estimated prior modes recover an analytical normal profile without R fitting', {
  x<-c(-2,1,3);se<-c(1,2,.5);truth<-10.25/5.25
  arguments<-list(betahat=x,sebetahat=se,method='shrink',mixcompdist='normal',mixsd=0,mode='estimate')
  answer<-.rsdeseq2_native_ash(arguments,x,se,0,FALSE)
  testthat::expect_equal(answer$fitted_g$mean,truth,tolerance=1e-8)
  testthat::expect_equal(answer$result$PosteriorMean,rep(truth,3),tolerance=1e-8)
  testthat::expect_equal(answer$result$PosteriorSD,rep(0,3),tolerance=1e-12)
  arguments$mode<-rep(1,3)
  testthat::expect_error(.rsdeseq2_native_ash(arguments,x,se,0,FALSE),'mode requires')
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

testthat::test_that('log normal tails retain underflowed probabilities through the native bridge', {
  x<-c(-Inf,-40,-8,0,8,40,Inf,NaN)
  actual<-.rsdeseq2_formula_request(list(operation='normal_log_probabilities',arguments=x))
  expected<-list(log_cdf=stats::pnorm(x,log.p=TRUE),log_upper=stats::pnorm(x,lower.tail=FALSE,log.p=TRUE))
  for(name in names(expected)){
    a<-actual[[name]];b<-expected[[name]];finite<-is.finite(b)
    testthat::expect_identical(is.na(a),is.na(b))
    testthat::expect_equal(a[!finite],b[!finite])
    testthat::expect_true(all(abs(a[finite]-b[finite])<=16*.Machine$double.eps*pmax(abs(b[finite]),1e-300)))
  }
  testthat::expect_true(is.finite(actual$log_cdf[[2L]]))
  testthat::expect_error(.rsdeseq2_formula_request(list(operation='normal_log_probabilities',arguments='invalid')))
})

testthat::test_that('uniform posterior transport preserves numerical boundaries and missing rows', {
  request<-list(operation='posterior_uniform_mixture',weights=1,lower=-10,upper=10,
    estimates=c(0,1),standard_errors=c(.5,.7),thresholds=c(-Inf,0,Inf))
  result<-.rsdeseq2_formula_request(request)
  testthat::expect_equal(result$posterior_mean,c(0,1),tolerance=2e-12)
  testthat::expect_equal(result$posterior_sd,c(.5,.7),tolerance=2e-12)
  testthat::expect_equal(result$cdf[c(2,5)],stats::pnorm(0,c(0,1),c(.5,.7)),tolerance=2e-12)
  testthat::expect_equal(result$log_likelihoods,rep(-log(20),2),tolerance=2e-12)
  missing<-utils::modifyList(request,list(lower=-1,upper=1,estimates=40,standard_errors=.01))
  missing<-.rsdeseq2_formula_request(missing)
  testthat::expect_true(all(is.nan(missing$posterior_mean)))
  testthat::expect_true(all(is.nan(missing$posterior_sd)))
  testthat::expect_true(all(is.nan(missing$cdf)))
  for(change in list(list(lower=11),list(weights=0),list(standard_errors=c(.5,-.7)),list(estimates=0)))
    testthat::expect_error(.rsdeseq2_formula_request(utils::modifyList(request,change)))
})
