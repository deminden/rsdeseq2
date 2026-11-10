# Inspect language and stored primitives only. A dispersion closure, active
# binding or delayed expression is never executed to discover its values.
.rsdeseq2_frozen_binding <- function(name, env) {
  while (!identical(env, emptyenv())) {
    if (exists(name, envir = env, inherits = FALSE)) {
      if (bindingIsActive(name, env)) stop("active frozen dispersion bindings are unsupported", call. = FALSE)
      value <- do.call(base::substitute, list(as.name(name), env))
      if (is.call(value) || is.symbol(value)) {
        stop("frozen dispersion constants must be stored primitive values in a local environment; delayed expressions are unsupported", call. = FALSE)
      }
      return(value)
    }
    env <- parent.env(env)
  }
  stop(sprintf("unknown frozen dispersion binding '%s'", name), call. = FALSE)
}

.rsdeseq2_saved_local_calls <- function(env) {
  trusted <- list(baseenv(), asNamespace("base"), asNamespace("stats"),
    asNamespace("DESeq2"), parent.env(asNamespace("DESeq2")))
  for (name in c("exp", "predict", "data.frame", "log")) {
    scope <- env
    while (!identical(scope, emptyenv()) && !exists(name, envir=scope, inherits=FALSE)) scope <- parent.env(scope)
    if (!any(vapply(trusted, identical, logical(1), y=scope))) {
      stop("masked saved local prediction functions are unsupported", call.=FALSE)
    }
    canonical <- if (name == "predict") stats::predict else get(name, baseenv(), inherits=FALSE)
    if (!identical(get(name, envir=scope, inherits=FALSE), canonical)) {
      stop("modified saved local prediction functions are unsupported", call.=FALSE)
    }
  }
}

.rsdeseq2_frozen_local_model <- function(fun) {
  if (!is.function(fun) || is.primitive(fun) || length(formals(fun)) != 1L ||
      identical(names(formals(fun)), "...")) {
    stop("frozen local dispersion functions require one mean argument", call. = FALSE)
  }
  expression <- body(fun)
  env <- environment(fun)
  argument <- names(formals(fun))[[1L]]
  standard <- quote(exp(predict(fit, data.frame(logMeans = log(means)))))
  if (identical(argument, "means") && identical(expression, standard)) {
    .rsdeseq2_saved_local_calls(env)
    fit <- .rsdeseq2_frozen_binding("fit", env)
    if (!is.list(fit) || !inherits(fit, "locfit") || isS4(fit)) {
      stop("saved local dispersion fit must be a locfit object", call. = FALSE)
    }
    fit <- unclass(fit)
    transform <- fit[["trans"]]
    if (!is.function(transform) || is.primitive(transform) || length(formals(transform)) != 1L ||
        identical(names(formals(transform)), "...") ||
        !identical(body(transform), as.name(names(formals(transform))[[1L]]))) {
      stop("saved local prediction requires an identity transform", call.=FALSE)
    }
    evaluation <- fit[["eva"]]
    if (!is.list(evaluation) || isS4(evaluation)) stop("invalid saved local evaluation data", call. = FALSE)
    evaluation <- unclass(evaluation)
    return(list(type = "saved", mi = fit[["mi"]], knots = evaluation[["xev"]],
      coefficients = evaluation[["coef"]], polynomial = evaluation[["pc"]],
      scale = evaluation[["scale"]], style = fit[["sty"]], derivative = fit[["deriv"]]))
  }
  # A single-expression block/return is syntax, not an R evaluation step.
  if (is.call(expression) && identical(expression[[1L]], as.name("{"))) {
    if (length(expression) != 2L) stop("frozen dispersion bodies must contain one expression", call. = FALSE)
    expression <- expression[[2L]]
  }
  if (is.call(expression) && identical(expression[[1L]], as.name("return"))) {
    if (length(expression) != 2L) stop("frozen dispersion return requires one expression", call. = FALSE)
    expression <- expression[[2L]]
  }
  formula <- structure(call("~", expression), class = "formula", .Environment = env)
  .rsdeseq2_validate_formula_calls(formula)
  names <- setdiff(all.vars(expression), argument)
  constants <- stats::setNames(lapply(names, function(name) {
    value <- .rsdeseq2_frozen_binding(name, env)
    if (!(is.numeric(value) || is.logical(value)) || is.object(value)) {
      stop("frozen dispersion constants must be primitive numeric or logical values", call. = FALSE)
    }
    value
  }), names)
  list(type = "expression", argument = argument,
    expression = .rsdeseq2_formula_ast(expression), constants = constants)
}

.rsdeseq2_frozen_local_curve <- function(model, counts, sizeFactors, normalizationFactors) {
  grid <- .Call("rsdeseq2_frozen_local_grid", counts, sizeFactors, normalizationFactors, PACKAGE = "rsdeseq2")
  .rsdeseq2_formula_request(list(operation = "frozen_local_grid", model = model, means = grid))
}
