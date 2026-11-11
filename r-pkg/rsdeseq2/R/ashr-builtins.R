# Read canonical constructor syntax and stored parameters; never execute its
# functions or reproduce its statistical calculations in the wrapper.
.rsdeseq2_native_ash_builtin <- function(lik) {
  kind <- lik[["name"]]
  if (identical(kind, "pois") || identical(kind, "binom")) return(.rsdeseq2_native_ash_count(lik))
  if (!identical(kind, "normalmix") && !identical(kind, "logF")) return(NULL)
  parameters <- if (kind == "normalmix") c("pilik", "sdlik") else c("df1", "df2")
  fields <- if (kind == "normalmix") c("lcdfFUN", "lpdfFUN", "etruncFUN", "e2truncFUN") else c("lcdfFUN", "lpdfFUN", "etruncFUN")
  functions <- lik[fields]
  if (!all(vapply(functions, function(f) typeof(f) == "closure", logical(1)))) return(NULL)
  env <- environment(functions[[1L]])
  if (!all(vapply(functions, function(f) identical(environment(f), env), logical(1)))) return(NULL)
  ns <- parent.env(env)
  if (!isNamespace(ns) || !identical(unname(getNamespaceName(ns)), "ashr") ||
      !all(ls(env, all.names = TRUE) %in% parameters)) return(NULL)
  constructor <- .rsdeseq2_frozen_binding(if (kind == "normalmix") "lik_normalmix" else "lik_logF", ns,
    context = "built-in likelihood constructor")
  if (typeof(constructor) != "closure") return(NULL)
  syntax <- body(constructor)
  if (is.call(syntax) && identical(syntax[[1L]], as.name("{")) && length(syntax) == 2L) syntax <- syntax[[2L]]
  if (!is.call(syntax) || !identical(syntax[[1L]], as.name("list"))) return(NULL)
  template <- as.list(syntax)
  for (field in fields) {
    expected <- template[[field]]
    if (!is.call(expected) || !identical(expected[[1L]], as.name("function")) ||
        !identical(names(formals(lik[[field]])), names(expected[[2L]])) ||
        !identical(body(lik[[field]]), expected[[3L]])) return(NULL)
  }
  if (kind == "logF" && !is.null(lik[["e2truncFUN"]])) return(NULL)
  values <- lapply(parameters, function(name) {
    value <- .rsdeseq2_frozen_binding(name, env, context = "built-in likelihood parameter")
    if (is.object(value) || !typeof(value) %in% c("integer", "double") ||
        (!is.null(dim(value)) && !(kind == "normalmix" && name == "sdlik" && length(dim(value)) == 2L))) {
      stop("built-in likelihood parameters require primitive numeric vectors or an SD matrix", call. = FALSE)
    }
    value
  })
  if (kind == "normalmix") list(type = "normalmix", weights = values[[1L]], scales = values[[2L]])
  else list(type = "logF", df1 = values[[1L]], df2 = values[[2L]])
}
