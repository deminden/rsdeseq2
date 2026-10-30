# Inspect language objects only. Never execute a formula to decide whether it is
# supported; even a nested or masked function must fail before model evaluation.
.rsdeseq2_validate_formula_calls <- function(formula) {
  builtins <- c("(", "I", "offset", "+", "-", "*", "/", "^", "%%", "%/%",
    "!", "&", "|", "<", "<=", ">", ">=", "==", "!=", "log", "log2", "log10",
    "log1p", "exp", "expm1", "sqrt", "abs", "sign", "floor", "ceiling", "trunc",
    "as.numeric", "as.double", "as.integer", "as.logical", "as.character", "is.na",
    "is.nan", "is.finite", "is.infinite", "c", "list", "cbind", "matrix", "rep",
    "seq_len", "length", "scale", "poly", "factor", "ordered", "as.factor",
    "as.ordered", "droplevels", "relevel", "mean", "sum", "sd", "ifelse",
    "~", ":", "%in%")
  env <- environment(formula)
  if (is.null(env)) env <- baseenv()
  trusted <- list(baseenv(), asNamespace("base"), asNamespace("stats"))
  for (name in intersect(search(), c("package:base", "package:stats"))) {
    trusted[[length(trusted)+1L]] <- as.environment(name)
  }
  reject <- function(name) stop(sprintf(
    "custom or unsupported formula function '%s'; use a supported built-in transform or supply a numeric design matrix",
    name), call.=FALSE)
  walk <- function(x, depth=0L) {
    if (depth > 128L) stop("formula expression nesting exceeds 128", call.=FALSE)
    if (!is.call(x)) return(invisible(NULL))
    head <- x[[1L]]
    namespace <- NULL
    if (is.symbol(head)) {
      name <- as.character(head)
    } else if (is.call(head) && length(head)==3L &&
        identical(head[[1L]], as.name("::")) &&
        is.symbol(head[[2L]]) && is.symbol(head[[3L]])) {
      namespace <- as.character(head[[2L]])
      name <- as.character(head[[3L]])
      if (!namespace %in% c("base", "stats")) reject(paste0(namespace, "::", name))
    } else reject("computed call")
    if (!name %in% builtins) reject(name)
    source <- if (exists(name, baseenv(), mode="function", inherits=FALSE)) baseenv() else asNamespace("stats")
    if (!exists(name, source, mode="function", inherits=FALSE)) reject(name)
    canonical <- get(name, source, mode="function", inherits=FALSE)
    if (is.null(namespace)) {
      scope <- env
      resolved <- NULL
      while (!identical(scope, emptyenv())) {
        if (exists(name, scope, inherits=FALSE)) {
          # Do not force a user binding to inspect it: even function lookup can
          # execute an active binding or a delayed definition. Aliases and
          # overrides must instead use an explicitly qualified built-in.
          if (!any(vapply(trusted, identical, logical(1), y=scope))) reject(name)
          resolved <- get(name, scope, mode="function", inherits=FALSE)
          break
        }
        scope <- parent.env(scope)
      }
    } else {
      resolved <- tryCatch(getExportedValue(namespace, name), error=function(e) NULL)
    }
    if (!identical(resolved, canonical)) reject(name)
    for (i in seq_len(length(x)-1L)) walk(x[[i+1L]], depth+1L)
    invisible(NULL)
  }
  walk(formula)
}
