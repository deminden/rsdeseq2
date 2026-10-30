# Language/data translation only. Numerical expression/frame/offset work shares
# the standalone Rust evaluator; no eval(), model.frame() or numerical fallback.
.rsdeseq2_formula_request <- function(request) {
  .Call("rsdeseq2_formula_request", request, PACKAGE="rsdeseq2")
}
.rsdeseq2_formula_collation <- function() {
  locale <- base::icuGetCollate("actual")
  if (locale %in% c("unknown","ICU not in use","ICU not available","C","POSIX")) return("C")
  if (identical(locale,"root")) return("root")
  gsub("_","-",locale,fixed=TRUE)
}
.rsdeseq2_formula_ast <- function(expression) {
  if (is.symbol(expression)) return(list(type="symbol",name=as.character(expression)))
  if (!is.call(expression)) return(list(type="literal",value=expression))
  head <- expression[[1L]]
  name <- if (is.symbol(head)) as.character(head) else paste0(
    as.character(head[[2L]]),"::",as.character(head[[3L]]))
  inputs <- as.list(expression)[-1L]
  labels <- names(inputs)
  list(type="call",function_name=name,arguments=lapply(seq_along(inputs),function(i) {
    list(name=if(is.null(labels)||labels[[i]]=="")NULL else labels[[i]],
      expression=.rsdeseq2_formula_ast(inputs[[i]]))
  }))
}
.rsdeseq2_formula_input <- function(expressions,data,env) {
  names <- unique(unlist(lapply(expressions,all.vars),use.names=FALSE))
  if (!length(names)) names <- character()
  bindings <- stats::setNames(lapply(names,function(name) {
    if (name %in% names(data)) data[[name]] else get(name,envir=env,inherits=TRUE)
  }),names)
  list(operation="frame_evaluate",
    expressions=lapply(expressions,.rsdeseq2_formula_ast),variables=bindings,
    collation=.rsdeseq2_formula_collation())
}
.rsdeseq2_formula_predict_call <- function(value,call) {
  if (!is.call(call)) return(call)
  head <- call[[1L]]
  name <- if(is.symbol(head))as.character(head) else if(is.call(head)&&identical(head[[1L]],as.name("::")))paste0(as.character(head[[2L]]),"::",as.character(head[[3L]])) else ""
  if (name=="scale") {
    center <- attr(value,"scaled:center");scale <- attr(value,"scaled:scale")
    if (!is.null(center)||!is.null(scale)) {
      if (!is.null(center)) call$center <- center
      if (!is.null(scale)) call$scale <- scale
    }
  } else if (name %in% c("poly","stats::poly") && inherits(value,"poly")) {
    coefficients <- attr(value,"coefs")
    call$coefs <- coefficients
  }
  call
}
.rsdeseq2_model_frame <- function(formula,data) {
  .rsdeseq2_validate_formula_calls(formula)
  terms <- stats::terms(formula,data=data)
  expressions <- as.list(attr(terms,"variables"))[-1L]
  labels <- vapply(expressions,deparse1,"",width.cutoff=500L)
  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()
  action <- attr(data,"na.action")
  if (is.null(action)||is.numeric(action)) action <- getOption("na.action")
  if (is.null(action)) action <- stats::na.fail
  if (is.character(action)) {
    if(length(action)!=1L||!action %in% c("na.omit","na.exclude","na.pass","na.fail")) stop("custom formula missing-value functions are unsupported",call.=FALSE)
    action <- get(action,asNamespace("stats"))
  }
  policy <- if(identical(action,stats::na.pass))"pass" else if(identical(action,stats::na.fail))"fail" else if(identical(action,stats::na.omit)||identical(action,stats::na.exclude))"omit" else stop("custom formula missing-value functions are unsupported",call.=FALSE)
  rows <- nrow(data)
  if (is.null(rows)) rows <- 0L
  request <- .rsdeseq2_formula_input(expressions,data,env)
  request$data_rows <- rows;request$na_action <- policy
  selected <- .rsdeseq2_formula_request(request)
  frame <- stats::setNames(selected$values,labels)
  originalNames <- attr(data,"row.names")
  if (length(originalNames)!=selected$original_rows) originalNames <- seq_len(selected$original_rows)
  rownames <- originalNames[selected$kept]
  attr(frame,"row.names") <- rownames
  class(frame) <- "data.frame"
  predictions <- attr(terms,"predvars")
  if (is.null(predictions)) {
    predictions <- attr(terms,"variables")
    for (i in seq_along(expressions)) predictions[[i+1L]] <- .rsdeseq2_formula_predict_call(frame[[i]],expressions[[i]])
    attr(terms,"predvars") <- predictions
  }
  attr(terms,"dataClasses") <- vapply(frame,function(x) {
    if(is.logical(x))"logical" else if(is.ordered(x))"ordered" else if(is.factor(x))"factor" else if(is.character(x))"character" else if(is.matrix(x)&&is.numeric(x))paste0("nmatrix.",ncol(x)) else if(is.numeric(x))"numeric" else "other"
  },"")
  attr(frame,"terms") <- terms
  if (length(selected$omitted)) {
    omitted <- stats::setNames(selected$omitted,originalNames[selected$omitted])
    class(omitted) <- if(identical(action,stats::na.exclude))"exclude" else "omit"
    attr(frame,"na.action") <- omitted
  }
  frame
}
.rsdeseq2_model_offsets <- function(frame) {
  indices <- attr(attr(frame,"terms"),"offset")
  if (is.null(indices)) indices <- integer()
  if (!length(indices)&&is.null(frame[["(offset)"]])) return(NULL)
  .rsdeseq2_formula_request(list(operation="offsets",values=unname(unclass(frame)),
    indices=indices,base=frame[["(offset)"]]))
}
.rsdeseq2_model_values <- function(frame,response) {
  .rsdeseq2_formula_request(list(operation="prepare",values=unname(unclass(frame)),
    response=response,collation=.rsdeseq2_formula_collation()))
}
