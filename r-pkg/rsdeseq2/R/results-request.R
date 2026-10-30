# Shared selector/metadata adapters for single and batch requests.
.rsdeseq2_result_parameters <- function(values, context = "results") {
  if (!is.list(values) || (length(values) && (is.null(names(values)) ||
      anyNA(names(values)) || any(!nzchar(names(values))) || anyDuplicated(names(values))))) {
    stop(sprintf("%s arguments must be a uniquely named list",context),call.=FALSE)
  }
  defaults <- list(contrast=NULL,name=NULL,listValues=c(1,-1),reference=NULL,
    lfcThreshold=0,altHypothesis="greaterAbs",useT=FALSE,degreesOfFreedom=NULL,
    alpha=.1,independentFiltering=TRUE,cooksCutoff=TRUE)
  unknown <- setdiff(names(values),names(defaults))
  if(length(unknown)) stop(sprintf("%s received unsupported arguments: %s",context,paste(unknown,collapse=", ")),call.=FALSE)
  for(name in names(values)) defaults[name] <- values[name]
  defaults$altHypothesis <- match.arg(defaults$altHypothesis,
    c("greaterAbs","lessAbs","greater","less","greaterAbs2014","greaterAbsUPSHOT"))
  defaults
}

.rsdeseq2_resolve_result_selector <- function(object,args,lrt) {
  names <- object$resultsNames
  reference <- args$reference %||% .rsdeseq2_reference_for_character_contrast(object,args$contrast)
  contrast <- .rsdeseq2_canonicalize_character_contrast(object,args$contrast)
  reference <- .rsdeseq2_canonicalize_character_contrast_reference(object,contrast,reference)
  if(is.null(contrast)) {
    name <- args$name
    if(is.null(name)) {
      if(length(names)<2L) stop("default results require at least two result names",call.=FALSE)
      name <- names[[length(names)]]
    }
    coefficient <- .rsdeseq2_resolve_results_name_index(name,names)
    weights <- numeric(length(names));weights[[coefficient]] <- 1;names(weights) <- names
    selector <- list(numeric=weights,resultName=names[[coefficient]],
      comparison=sprintf("coefficient %s",names[[coefficient]]),allZero=list(type="none"))
  } else {
    selector <- resolveResultsContrastRust(contrast,names,args$listValues,reference)
  }
  if(lrt) selector$comparison <- sprintf("%s; LRT full vs reduced",selector$comparison)
  selector
}

.rsdeseq2_prepare_result_request <- function(object,args,lrt) {
  if(lrt && (!is.numeric(args$lfcThreshold) || length(args$lfcThreshold)!=1L ||
      !is.finite(args$lfcThreshold) || args$lfcThreshold!=0 ||
      !identical(args$altHypothesis,"greaterAbs") || isTRUE(args$useT) || !is.null(args$degreesOfFreedom))) {
    stop("lfcThreshold, altHypothesis, useT, and degreesOfFreedom are only supported for Wald results",call.=FALSE)
  }
  selector <- .rsdeseq2_resolve_result_selector(object,args,lrt)
  options <- if(lrt)NULL else .rsdeseq2_validate_wald_result_options(
    args$lfcThreshold,args$altHypothesis,args$useT,args$degreesOfFreedom,nrow(object$beta))
  list(selector=selector,wald=options,lrt=lrt,args=args)
}

.rsdeseq2_native_result_request <- function(object,prepared) {
  options <- prepared$wald
  list(contrast=as.numeric(prepared$selector$numeric),
    wald=if(is.null(options))NULL else list(threshold=options$lfcThreshold,
      alternative=.rsdeseq2_wald_alternative_code(options$altHypothesis),useT=options$useT,df=options$degreesOfFreedom),
    zero=.rsdeseq2_contrast_zero_context(object,prepared$selector$allZero))
}

.rsdeseq2_format_primitive_result <- function(object,prepared,raw,adjust=TRUE) {
  selector <- prepared$selector; options <- prepared$wald; args <- prepared$args
  out <- resultsTableRust(object$baseMean,raw$log2FoldChange,raw$lfcSE,
    if(prepared$lrt)object$lrtStat else raw$stat,
    if(prepared$lrt)object$lrtPvalue else raw$pvalue,
    padj=rep(NA_real_,length(object$baseMean)),contrast=selector$numeric,rowNames=object$rowNames)
  attr(out,"resultName") <- selector$resultName
  attr(out,"comparison") <- selector$comparison
  attr(out,"contrastAllZero") <- selector$allZero$type
  if(prepared$lrt) {
    if(!is.null(object$lrtDf)) attr(out,"lrtDf") <- object$lrtDf
  } else {
    attr(out,"lfcThreshold") <- options$lfcThreshold
    attr(out,"altHypothesis") <- options$altHypothesis
    if(!is.null(options$degreesOfFreedom)) attr(out,"degreesOfFreedom") <- options$degreesOfFreedom
  }
  if(adjust) applyIndependentFilteringRust(out,alpha=args$alpha,enabled=args$independentFiltering) else out
}

.rsdeseq2_cached_result_request <- function(args) {
  is.null(args$contrast) && is.null(args$name) && identical(as.numeric(args$lfcThreshold),0) &&
    identical(args$altHypothesis,"greaterAbs") && !isTRUE(args$useT) && is.null(args$degreesOfFreedom)
}

.rsdeseq2_finish_fitted_result <- function(fit,args,out, cached=FALSE) {
  out <- .rsdeseq2_fit_cooks_policy(out,fit,args$cooksCutoff,cached=cached)
  out <- applyIndependentFilteringRust(out,alpha=args$alpha,enabled=args$independentFiltering)
  if(cached) {
    if(!isTRUE(args$independentFiltering)) out$filtered <- NULL
    resultName <- utils::tail(fit$resultsNames,1L)
    comparison <- sprintf("coefficient %s",resultName)
    if(identical(fit$test,"LRT")) comparison <- sprintf("%s; LRT full vs reduced",comparison)
    metadata <- list(resultName=resultName,comparison=comparison,test=fit$test,
      contrastAllZero="none",lfcThreshold=0,altHypothesis="greaterAbs",contrast=NULL,
      alpha=args$alpha,lrtDf=fit$lrtDf,independentFiltering=attr(out,"independentFiltering"),
      cooksReplacement=.rsdeseq2_cooks_replacement_metadata(fit))
  } else {
    metadata <- list(resultName=attr(out,"resultName"),comparison=attr(out,"comparison"),
      test=fit$test,contrastAllZero=attr(out,"contrastAllZero"),contrast=attr(out,"contrast"),
      lfcThreshold=attr(out,"lfcThreshold"),altHypothesis=attr(out,"altHypothesis"),
      degreesOfFreedom=attr(out,"degreesOfFreedom"),lrtDf=attr(out,"lrtDf"),alpha=args$alpha,
      independentFiltering=attr(out,"independentFiltering"))
    metadata$cooksReplacement <- .rsdeseq2_cooks_replacement_metadata(fit)
  }
  .rsdeseq2_finalize_result_metadata(out,metadata,fit)
}
