# Translate syntax and captured scalar primitives. Never call the R functions.
.rsdeseq2_restricted_likelihood_function <- function(fun, arity, label) {
  if (!is.function(fun) || is.primitive(fun) || length(formals(fun)) != arity ||
      "..." %in% names(formals(fun))) {
    stop(sprintf("%s requires %s scalar arguments", label, arity), call. = FALSE)
  }
  expression <- body(fun)
  if (is.call(expression) && identical(expression[[1L]], as.name("{"))) {
    if (length(expression) != 2L) stop(sprintf("%s must contain one expression", label), call. = FALSE)
    expression <- expression[[2L]]
  }
  if (is.call(expression) && identical(expression[[1L]], as.name("return"))) {
    if (length(expression) != 2L) stop(sprintf("%s requires one return expression", label), call. = FALSE)
    expression <- expression[[2L]]
  }
  allowed <- c("(", "I", "+", "-", "*", "/", "^", "%%", "%/%", "!", "&", "|",
    "<", "<=", ">", ">=", "==", "!=", "log", "log2", "log10", "log1p", "exp",
    "expm1", "sqrt", "abs", "sign", "floor", "ceiling", "trunc", "ifelse", "pmin", "pmax")
  inspect <- function(x, depth = 0L) {
    if (depth > 128L) stop(sprintf("%s nesting exceeds 128", label), call. = FALSE)
    if (!is.call(x)) return(invisible(NULL))
    head <- x[[1L]]
    name <- if (is.symbol(head)) as.character(head) else if (
      is.call(head) && length(head) == 3L && identical(head[[1L]], as.name("::")) &&
      identical(head[[2L]], as.name("base")) && is.symbol(head[[3L]])) {
      as.character(head[[3L]])
    } else "computed call"
    if (!name %in% allowed) stop(sprintf("unsupported %s call '%s'", label, name), call. = FALSE)
    for (i in seq_len(length(x) - 1L)) inspect(x[[i + 1L]], depth + 1L)
    invisible(NULL)
  }
  inspect(expression)
  env <- environment(fun)
  formula <- structure(call("~", expression), class = "formula", .Environment = env)
  .rsdeseq2_validate_formula_calls(formula)
  arguments <- names(formals(fun))
  captured <- setdiff(all.vars(expression), arguments)
  constants <- stats::setNames(lapply(captured, function(name) {
    value <- .rsdeseq2_frozen_binding(name, env, context = label)
    if (is.object(value) || !typeof(value) %in% c("integer", "double", "logical") ||
        length(value) != 1L || !is.null(dim(value))) {
      stop(sprintf("%s constants must be stored numeric or logical scalars", label), call. = FALSE)
    }
    value
  }), captured)
  list(arguments = arguments, expression = .rsdeseq2_formula_ast(expression), constants = constants)
}

.rsdeseq2_restricted_likelihood <- function(likelihood) {
  if (isS4(likelihood) || typeof(likelihood) != "list") stop("custom likelihood must be a primitive list", call. = FALSE)
  likelihood <- unclass(likelihood)
  name <- likelihood[["name"]]
  if (is.object(name) || typeof(name) != "character" || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("custom likelihood requires one nonempty name", call. = FALSE)
  }
  optional <- function(field) {
    value <- likelihood[[field]]
    if (is.null(value)) return(NULL)
    .rsdeseq2_restricted_likelihood_function(value, 2L, paste("custom likelihood", field))
  }
  list(name = name,
    log_cdf = .rsdeseq2_restricted_likelihood_function(likelihood[["lcdfFUN"]], 1L, "custom likelihood log-CDF"),
    log_density = .rsdeseq2_restricted_likelihood_function(likelihood[["lpdfFUN"]], 1L, "custom likelihood log-density"),
    first = optional("etruncFUN"), second = optional("e2truncFUN"))
}
