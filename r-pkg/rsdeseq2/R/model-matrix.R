# Formula values, factor coding and products share standalone Rust kernels.
.rsdeseq2_model_matrix <- function(formula, data) {
  .rsdeseq2_validate_formula_calls(formula)
  frame <- if (is.null(attr(data, "terms"))) .rsdeseq2_model_frame(formula, data) else data
  terms <- attr(frame, "terms")
  codes <- attr(terms, "factors")
  if (!length(codes)) codes <- matrix(integer(), nrow=length(frame), ncol=0L)
  labels <- if (is.null(rownames(codes))) names(frame) else rownames(codes)
  contrastMetadata <- list()
  defaults <- as.character(getOption("contrasts"))
  response <- attr(terms,"response")
  if (is.null(response)) response <- 0L
  prepared <- .rsdeseq2_model_values(frame,response)
  variables <- lapply(seq_along(frame), function(i) {
    value <- prepared[[i]]
    name <- names(frame)[[i]]
    if (is.factor(value)) {
      coding <- attr(value, "contrasts")
      if (is.null(coding)) coding <- defaults[[1L + is.ordered(value)]]
      contrastMetadata[[name]] <<- coding
      if (is.matrix(coding)) {
        storage.mode(coding) <- "double"
        index <- 5L
      } else {
        index <- match(coding, c("contr.treatment", "contr.sum", "contr.helmert", "contr.poly", "contr.SAS")) - 1L
        if (length(index)!=1L || is.na(index)) stop("factor contrast functions must use a supported Rust coding or an explicit contrast matrix", call.=FALSE)
      }
      list(name=labels[[i]], kind=1L, values=as.integer(value), levels=levels(value),
        coding=index, contrasts=if(is.matrix(coding)) coding else NULL,
        columnNames=if(is.matrix(coding)) colnames(coding) else NULL)
    } else {
      if (!is.numeric(value) && !is.double(value) && !is.integer(value)) stop("model variables must be numeric or categorical", call.=FALSE)
      if (!is.matrix(value)) value <- matrix(value, ncol=1L)
      storage.mode(value) <- "double"
      list(name=labels[[i]], kind=0L, values=value, columnNames=colnames(value))
    }
  })
  result <- .Call("rsdeseq2_model_matrix", as.integer(nrow(frame)), variables,
    codes, as.logical(attr(terms,"intercept")), as.integer(response), PACKAGE="rsdeseq2")
  rownames(result) <- rownames(frame)
  if (length(contrastMetadata)) attr(result,"contrasts") <- contrastMetadata
  result
}
