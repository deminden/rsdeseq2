# Bounded result requests with one fitted-dataset validation per public call.
resultsList <- function(object, requests, ..., batchSize = 32L) {
  if(!is.list(requests)) stop("requests must be a list of named argument lists",call.=FALSE)
  if(!is.numeric(batchSize)||length(batchSize)!=1L||!is.finite(batchSize)||
      batchSize<1||batchSize>64||batchSize!=floor(batchSize)) {
    stop("batchSize must be an integer from 1 to 64",call.=FALSE)
  }
  batchSize <- as.integer(batchSize)
  common <- list(...)
  defaults <- .rsdeseq2_result_parameters(common,"resultsList")
  args <- lapply(seq_along(requests),function(i) {
    request <- requests[[i]]
    # Validate names and unknown arguments before merging; preserve explicit NULL.
    .rsdeseq2_result_parameters(request,sprintf("request %d",i))
    values <- defaults
    for(name in names(request)) values[name] <- request[name]
    values$altHypothesis <- match.arg(values$altHypothesis,
      c("greaterAbs","lessAbs","greater","less","greaterAbs2014","greaterAbsUPSHOT"))
    values
  })
  if(inherits(object,"DESeqDataSet")) {
    object <- .rsdeseq2_fitted_view(object)
    if(is.null(object)) stop("DESeqDataSet must be fitted with rsdeseq2::DESeq() before calling resultsList()",call.=FALSE)
  }
  fitted <- inherits(object,"rsdeseq2RustFit")
  if(!fitted && !inherits(object,c("rsdeseq2PrimitiveWaldFit","rsdeseq2PrimitiveLrtFit"))) {
    if(.rsdeseq2_is_list_like_lrt_input(object)) object <- nbinomLRTRust(object=object) else
      if(.rsdeseq2_is_list_like_wald_input(object)) object <- nbinomWaldTestRust(object=object) else
        .rsdeseq2_not_available("resultsList")
  }
  for(i in seq_along(args)) {
    cutoff <- args[[i]]$cooksCutoff
    if(!identical(cutoff,TRUE)) {
      validated <- .rsdeseq2_validate_cooks_cutoff(cutoff)
      if(!is.null(validated)&&validated<0) stop("cooksCutoff must be non-negative",call.=FALSE)
    }
    explicit <- "cooksCutoff" %in% c(names(common),names(requests[[i]]))
    if(!fitted && explicit && !identical(cutoff,FALSE) && !is.null(cutoff)) {
      stop("Cook's filtering requires a fitted dataset with retained Cook's context",call.=FALSE)
    }
    args[[i]]$alpha <- .rsdeseq2_validate_independent_filter_alpha(args[[i]]$alpha)
  }
  output <- vector("list",length(requests)); names(output) <- names(requests)
  if(!length(requests)) return(output)
  cached <- vapply(args,function(x)fitted&&.rsdeseq2_cached_result_request(x),logical(1))
  selected <- which(!cached)
  if(length(selected)) {
    lrt <- if(fitted)identical(object$test,"LRT") else inherits(object,"rsdeseq2PrimitiveLrtFit")
    primitive <- if(!fitted)object else if(lrt)nbinomLRTRust(object=object) else nbinomWaldTestRust(object=object)
    prepared <- lapply(selected,function(i).rsdeseq2_prepare_result_request(primitive,args[[i]],lrt))
    for(start in seq.int(1L,length(selected),by=batchSize)) {
      chunk <- seq.int(start,min(start+batchSize-1L,length(selected)))
      native <- lapply(prepared[chunk],function(p).rsdeseq2_native_result_request(primitive,p))
      raw <- .rsdeseq2_evaluate_result_batch(primitive,native)
      for(position in seq_along(chunk)) {
        index <- chunk[[position]]; request <- selected[[index]]
        value <- .rsdeseq2_format_primitive_result(primitive,prepared[[index]],raw[[position]],adjust=!fitted)
        output[[request]] <- if(fitted).rsdeseq2_finish_fitted_result(object,args[[request]],value) else value
      }
    }
  }
  for(i in which(cached)) output[[i]] <- .rsdeseq2_finish_fitted_result(object,args[[i]],object$results,cached=TRUE)
  output
}
