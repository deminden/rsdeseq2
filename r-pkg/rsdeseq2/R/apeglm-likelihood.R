# Import likelihood syntax and stored scalar captures; never invoke a closure.
.rsdeseq2_apeglm_likelihood <- function(fun) {
  if (typeof(fun) != "closure") stop("general apeglm requires a likelihood function", call. = FALSE)
  standard <- quote({
    xbeta <- x %*% beta + offset
    mean.hat <- exp(xbeta)
    dnbinom(y, mu = mean.hat, size = 1/param, log = TRUE)
  })
  env <- environment(fun)
  if (identical(names(formals(fun)), c("y", "x", "beta", "param", "offset")) &&
      identical(body(fun), standard) && isNamespace(env) &&
      identical(unname(getNamespaceName(env)), "apeglm")) {
    return(list(type = "negative_binomial"))
  }
  list(type = "restricted", model = .rsdeseq2_restricted_likelihood_function(
    fun, 5L, "general apeglm likelihood", regression = TRUE))
}
