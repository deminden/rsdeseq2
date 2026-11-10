testthat::test_that('random-state import preserves stored integer words without calculations', {
  env <- new.env(parent=baseenv())
  inspect <- rsdeseq2:::.rsdeseq2_posterior_rng_state
  testthat::expect_null(inspect(env))
  state <- c(10403L,1L,NA_integer_,rep(1L,623L))
  env$.Random.seed <- state
  testthat::expect_identical(inspect(env),state)
  testthat::expect_identical(env$.Random.seed,state)
  for (bad in list(as.numeric(state),matrix(state,ncol=1L),structure(state,class='custom'))) {
    env$.Random.seed <- bad
    testthat::expect_error(inspect(env),'stored primitive integer vector')
  }
})

testthat::test_that('random-state import rejects executable and unwritable bindings', {
  env <- new.env(parent=baseenv())
  inspect <- rsdeseq2:::.rsdeseq2_posterior_rng_state
  touched <- FALSE
  makeActiveBinding('.Random.seed',function() {touched <<- TRUE; stop('must not execute')},env)
  testthat::expect_error(inspect(env),'active')
  testthat::expect_false(touched)
  rm('.Random.seed',envir=env)
  delayedAssign('.Random.seed',{touched <<- TRUE; rep(1L,626L)},assign.env=env)
  testthat::expect_error(inspect(env),'delayed')
  testthat::expect_false(touched)
  rm('.Random.seed',envir=env)
  env$.Random.seed <- c(10403L,624L,rep(1L,624L))
  lockBinding('.Random.seed',env)
  testthat::expect_error(inspect(env),'must be writable')
})
