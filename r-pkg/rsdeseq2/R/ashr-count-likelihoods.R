# Canonical syntax and stored data only; all derived-parameter checks run in Rust.
.rsdeseq2_native_ash_count <- function(lik) {
  kind <- lik[["name"]]
  if (!kind %in% c("pois", "binom") || !identical(lik[["const"]], TRUE)) return(NULL)
  fields <- c("lcdfFUN", "lpdfFUN", "etruncFUN", "e2truncFUN")
  functions <- lik[fields]
  if (!all(vapply(functions, function(f) typeof(f) == "closure", logical(1)))) return(NULL)
  env <- environment(functions[[1L]])
  if (!all(vapply(functions, function(f) identical(environment(f), env), logical(1)))) return(NULL)
  ns <- parent.env(env)
  parameters <- if (kind == "pois") c("y", "scale") else c("y", "n")
  if (!isNamespace(ns) || !identical(unname(getNamespaceName(ns)), "ashr") ||
      !all(ls(env, all.names = TRUE) %in% c(parameters, "link", "y1", "n1"))) return(NULL)
  constructor <- .rsdeseq2_frozen_binding(if (kind == "pois") "lik_pois" else "lik_binom", ns,
    context = "count likelihood constructor")
  if (typeof(constructor) != "closure") return(NULL)
  templates <- function(expression) {
    if (!is.call(expression)) return(list())
    pieces <- as.list(expression)
    if (identical(pieces[[1L]], as.name("list")) && all(fields %in% names(pieces))) return(list(pieces))
    unlist(lapply(pieces[-1L], templates), recursive = FALSE)
  }
  candidates <- templates(body(constructor))
  matches <- vapply(candidates, function(template) all(vapply(fields, function(field) {
    expected <- template[[field]]
    is.call(expected) && identical(expected[[1L]], as.name("function")) &&
      identical(formals(lik[[field]]), expected[[2L]]) && identical(body(lik[[field]]), expected[[3L]])
  }, logical(1))), logical(1))
  link <- .rsdeseq2_frozen_binding("link", env, context = "count likelihood link")
  choices <- if (kind == "pois") c("identity", "log") else c("identity", "logit")
  if (length(candidates) != 2L || !identical(which(matches), match(link, choices))) return(NULL)
  read <- function(name) {
    value <- .rsdeseq2_frozen_binding(name, env, context = "count likelihood parameter")
    if (is.object(value) || !typeof(value) %in% c("integer", "double") || !is.null(dim(value))) {
      stop("count likelihood parameters require stored primitive numeric vectors", call. = FALSE)
    }
    value
  }
  values <- stats::setNames(lapply(parameters, read), parameters)
  data <- lik[["data"]]
  if (typeof(data) != "list" || is.object(data) || !identical(names(data), c(parameters, "link")) ||
      !identical(data$link, link) || !all(vapply(parameters, function(n) identical(data[[n]], values[[n]]), logical(1)))) {
    stop("count likelihood metadata differs from stored parameters", call. = FALSE)
  }
  result <- list(type = kind, y = values$y, link = link)
  if (kind == "pois") result$scale <- values$scale else result$trials <- values$n
  if (link != "identity") {
    result$capture_y <- read("y1")
    if (kind == "binom") result$capture_trials <- read("n1")
  }
  result
}
