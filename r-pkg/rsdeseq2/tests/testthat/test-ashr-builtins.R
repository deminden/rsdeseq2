testthat::test_that('canonical built-in likelihoods import only stored primitive parameters', {
  testthat::skip_if_not_installed('ashr')
  matrix<-cbind(c(.5,.7,1),c(1,2,3))
  input<-ashr::lik_normalmix(c(.7,.3),matrix)
  imported<-.rsdeseq2_native_ash_likelihood(input)
  testthat::expect_identical(imported,list(type='normalmix',weights=c(.7,.3),scales=matrix))
  answer<-.rsdeseq2_formula_request(list(operation='posterior_builtin_likelihood',likelihood=imported,
    arguments=rep(0,3),lower=rep(-1,3),upper=rep(1,3)))
  expected<-log(.7*stats::dnorm(0,sd=matrix[,1])+.3*stats::dnorm(0,sd=matrix[,2]))
  testthat::expect_equal(answer$log_density,expected,tolerance=2e-12)
  testthat::expect_equal(answer$first,rep(0,3),tolerance=2e-12)
  imported<-.rsdeseq2_native_ash_likelihood(ashr::lik_logF(2,2))
  answer<-.rsdeseq2_formula_request(list(operation='posterior_builtin_likelihood',likelihood=imported,
    arguments=0,lower=-2,upper=2))
  testthat::expect_equal(answer$log_cdf,-log(2),tolerance=2e-12)
  testthat::expect_equal(answer$log_density,-log(4),tolerance=2e-12)
  testthat::expect_equal(answer$first,0,tolerance=2e-12)
  testthat::expect_true(is.na(answer$second))
})

testthat::test_that('modified functions and active or delayed built-in captures never execute', {
  testthat::skip_if_not_installed('ashr')
  called<-FALSE
  input<-ashr::lik_normalmix(c(.7,.3),c(.5,2))
  input$lpdfFUN<-function(x){called<<-TRUE;stop('executed modified likelihood')}
  testthat::expect_error(.rsdeseq2_native_ash_likelihood(input),'unsupported|one expression')
  testthat::expect_false(called)
  input<-ashr::lik_normalmix(c(.7,.3),c(.5,2));env<-environment(input$lpdfFUN)
  rm('pilik',envir=env);makeActiveBinding('pilik',function(){called<<-TRUE;c(.7,.3)},env)
  testthat::expect_error(.rsdeseq2_native_ash_likelihood(input),'active')
  testthat::expect_false(called)
  input<-ashr::lik_normalmix(c(.7,.3),c(.5,2));env<-environment(input$lpdfFUN)
  delayedAssign('pilik',{called<<-TRUE;c(.7,.3)},assign.env=env)
  testthat::expect_error(.rsdeseq2_native_ash_likelihood(input),'delayed|stored')
  testthat::expect_false(called)
})

testthat::test_that('equal-length matrix captures preserve their column-wise global mixture', {
  testthat::skip_if_not_installed('ashr')
  weights<-matrix(c(.9,.1,.2,.8),nrow=2,byrow=TRUE)
  scales<-matrix(c(.5,1,2,3),nrow=2,byrow=TRUE)
  payload<-.rsdeseq2_native_ash_likelihood(ashr::lik_normalmix(weights,scales))
  testthat::expect_identical(payload$weights,weights)
  testthat::expect_identical(payload$scales,scales)
  answer<-.rsdeseq2_formula_request(list(operation='posterior_builtin_likelihood',likelihood=payload,
    arguments=rep(0,3),lower=rep(-1,3),upper=rep(1,3)))
  expected<-sum(as.numeric(weights)*stats::dnorm(0,sd=as.numeric(scales)))
  testthat::expect_equal(exp(answer$log_density),rep(expected,3),tolerance=2e-12)
  testthat::expect_equal(answer$log_cdf,rep(0,3),tolerance=2e-12)
})
