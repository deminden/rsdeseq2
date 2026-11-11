# Object/syntax translation only. Fitting, summaries and threshold arithmetic
# are performed by the standalone Rust engine through one native request.
.rsdeseq2_native_ash_prior <- function(g) {
  if (isS4(g) || typeof(g) != "list") stop("native ashr prior requires a mixture data object", call. = FALSE)
  classes <- class(g)
  kind <- classes[classes %in% c("normalmix", "unimix", "tnormalmix")]
  if (length(kind) != 1L) stop("unsupported native ashr mixture class", call. = FALSE)
  g <- unclass(g)
  field <- function(name) {
    value <- g[[name]]
    if (is.object(value) || !typeof(value) %in% c("integer", "double") || !is.null(dim(value))) {
      stop(sprintf("mixture field '%s' must be a primitive numeric vector", name), call. = FALSE)
    }
    value
  }
  list(component_class = kind, weights = field("pi"),
    means = if (kind == "unimix") NULL else field("mean"),
    sds = if (kind == "unimix") NULL else field("sd"),
    lower = if (kind == "normalmix") NULL else field("a"),
    upper = if (kind == "normalmix") NULL else field("b"))
}
.rsdeseq2_native_ash_likelihood <- function(lik) {
  if (isS4(lik) || typeof(lik) != "list") stop("likelihood requires a primitive list", call. = FALSE)
  lik <- unclass(lik)
  builtin <- .rsdeseq2_native_ash_builtin(lik)
  if (!is.null(builtin)) return(builtin)
  unwrap <- function(f, arguments = "x") {
    if (typeof(f) != "closure" || !identical(names(formals(f)), arguments)) return(NULL)
    body <- body(f)
    if (is.call(body) && identical(body[[1L]], as.name("{")) && length(body) == 2L) body <- body[[2L]]
    body
  }
  cdf <- unwrap(lik[["lcdfFUN"]])
  pdf <- unwrap(lik[["lpdfFUN"]])
  first <- unwrap(lik[["etruncFUN"]], c("a", "b"))
  second <- unwrap(lik[["e2truncFUN"]], c("a", "b"))
  standard_scope <- function() {
    functions <- lik[c("lcdfFUN", "lpdfFUN", "etruncFUN", "e2truncFUN")]
    if (!all(vapply(functions, function(f) typeof(f) == "closure", logical(1)))) return(FALSE)
    env <- environment(functions[[1L]])
    if (!all(vapply(functions, function(f) identical(environment(f), env), logical(1)))) return(FALSE)
    parent <- parent.env(env)
    isNamespace(parent) && identical(unname(getNamespaceName(parent)), "ashr") &&
      !any(vapply(c("my_etruncnorm", "my_e2truncnorm"), exists, logical(1), envir = env, inherits = FALSE))
  }
  if (identical(lik[["name"]], "normal") &&
      identical(cdf, quote(stats::pnorm(x, log = TRUE))) &&
      identical(pdf, quote(stats::dnorm(x, log = TRUE))) &&
      identical(first, quote(my_etruncnorm(a, b))) &&
      identical(second, quote(my_e2truncnorm(a, b))) && standard_scope()) return(list(type = "normal"))
  if (identical(lik[["name"]], "t") &&
      identical(cdf, quote(stats::pt(x, df = df, log = TRUE))) &&
      identical(pdf, quote(stats::dt(x, df = df, log = TRUE))) &&
      identical(first, quote(etrunct::e_trunct(a, b, df = df, r = 1))) &&
      identical(second, quote(etrunct::e_trunct(a, b, df = df, r = 2))) && standard_scope()) {
    df <- .rsdeseq2_frozen_binding("df", environment(lik[["lcdfFUN"]]), context = "Student-t likelihood")
    return(list(type = "t", df = df))
  }
  list(type = "restricted", model = .rsdeseq2_restricted_likelihood(lik))
}
.rsdeseq2_native_ash <- function(arguments, original_estimates, original_se, lfcThreshold, svalue) {
  if (!is.null(arguments$g)) arguments$g <- .rsdeseq2_native_ash_prior(arguments$g)
  if (!is.null(arguments$lik)) arguments$lik <- .rsdeseq2_native_ash_likelihood(arguments$lik)
  native <- .rsdeseq2_formula_request(list(operation = "ashr_shrinkage", arguments = arguments,
    query_estimates = original_estimates, query_standard_errors = original_se,
    lfc_threshold = lfcThreshold, svalue = svalue))
  for (message in native$warnings) warning(message, call. = FALSE)
  prior <- native$fitted_prior
  fitted <- if (is.null(prior)) NULL else switch(prior$component_class,
    normalmix = data.frame(pi = prior$weights, mean = prior$means, sd = prior$sds),
    unimix = data.frame(pi = prior$weights, a = prior$lower, b = prior$upper),
    tnormalmix = data.frame(pi = prior$weights, mean = prior$means, sd = prior$sds, a = prior$lower, b = prior$upper))
  if (!is.null(prior)) class(fitted) <- prior$component_class
  list(result = list(PosteriorMean = native$posterior_mean, PosteriorSD = native$posterior_sd,
      svalue = native$svalue), fitted_g = fitted)
}
