#!/usr/bin/env Rscript
# Complete public formula fits and later queries, not R-prepared matrix input.
args <- commandArgs(TRUE)
if (length(args) != 3L) stop("usage: STANDALONE_DRIVER R_LIBRARY FRESH_OUTPUT_DIR")
driver <- normalizePath(args[[1L]])
library(rsdeseq2, lib.loc=normalizePath(args[[2L]]))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~", " ", script, fixed=TRUE)), "..", ".."))
source(file.path(root, "r-pkg/rsdeseq2/tests/testthat/helper-fit-data.R"))
source(file.path(root, "r-pkg/rsdeseq2/tests/testthat/helper-formula.R"))
source(file.path(root, "scripts/validation/native_posterior_check_helpers.R"))
source(file.path(root, "scripts/validation/fresh_fit_check_helpers.R"))
out <- args[[3L]]
if (dir.exists(out)) stop("Use a fresh output directory")
dir.create(out, recursive=TRUE)
counts <- .rsdeseq2_test_counts(); counts[40L, ] <- 0
replaced <- counts; replaced[1L,8L] <- 2000; replaced[2L,] <- 0; replaced[2L,8L] <- 1000
size <- exp(seq(log(.8),log(1.2),length.out=8L))
factors <- outer(seq(.9,1.1,length.out=40L),size)
profiles <- list(
  factor=list(full="~condition", columns=2L, character=c("condition","B","A")),
  numeric_first=list(full="~exposure+condition", columns=3L, character=c("condition","B","A")),
  factor_first=list(full="~condition+exposure", columns=3L, character=c("condition","B","A")),
  interaction=list(full="~condition*exposure", columns=4L),
  numeric_first_interaction=list(full="~exposure*condition", columns=4L),
  relevel=list(full="~relevel(condition,ref='B')", columns=2L),
  no_intercept=list(full="~0+condition", columns=2L),
  polynomial=list(full="~condition+poly(exposure,2)", columns=4L),
  three_levels=list(full="~condition", columns=3L, three=TRUE, character=c("condition","C","B")),
  ordered=list(full="~ordered(condition)", columns=3L, three=TRUE),
  polynomial_first_interaction=list(full="~poly(exposure,2):condition", columns=7L, three=TRUE),
  factor_first_polynomial_interaction=list(full="~condition:poly(exposure,2)", columns=7L, three=TRUE),
  polynomial_main_and_interaction=list(full="~poly(exposure,2)*condition", columns=9L, three=TRUE, four=TRUE, character=c("condition","C","B")),
  repeated_polynomial=list(full="~poly(exposure,2)+poly(exposure, 2)", columns=3L),
  aliased_polynomial=list(full="~poly(exposure,2)+poly(exposure,degree=2)", columns=5L, error="rank"),
  repeated_identity=list(full="~I(exposure)+I( exposure)", columns=2L),
  removed_identity_alias=list(full="~I(exposure)+I(+exposure)-I(exposure)", columns=2L),
  identity_alias_interaction=list(full="~I(exposure)+condition:I(+exposure)", columns=4L, error="rank"),
  removed_interaction_alias=list(full="~condition:I(exposure)+condition:I(+exposure)-condition:I(exposure)", columns=3L),
  removed_polynomial_alias=list(full="~poly(exposure,2)+poly(exposure,degree=2)-poly(exposure,2)", columns=3L),
  aliased_identity=list(full="~I(exposure)+I(+exposure)", columns=3L, error="rank"),
  offset_estimated=list(full="~condition+offset(log(exposure))", columns=2L),
  offset_supplied=list(full="~condition+offset(log(exposure))", columns=2L, size=TRUE),
  offset_factors=list(full="~condition+offset(log(exposure))", columns=2L, factors=TRUE),
  factors_and_size=list(full="~condition", columns=2L, factors=TRUE, size=TRUE),
  offset_factors_and_size=list(full="~condition+offset(log(exposure))", columns=2L, factors=TRUE, size=TRUE),
  numeric_first_lrt=list(full="~exposure+condition", reduced="~exposure", columns=3L),
  lrt=list(full="~condition", reduced="~1", columns=2L, character=c("condition","B","A")),
  lrt_offset=list(full="~condition+offset(log(exposure))", reduced="~1", columns=2L, size=TRUE),
  lrt_matching_offset=list(full="~condition+offset(log(exposure))", reduced="~offset(log(exposure))", columns=2L, size=TRUE),
  lrt_matching_offset_replacement=list(full="~condition+offset(log(exposure))", reduced="~offset(log(exposure))", columns=2L, size=TRUE, replace=TRUE),
  replacement=list(full="~condition", columns=2L, replace=TRUE, character=c("condition","B","A")),
  lrt_mismatched_offset=list(full="~condition+offset(log(exposure))", reduced="~offset(exposure)", columns=2L, error="offset")
)
data_for <- function(p) {
  n <- if(isTRUE(p$four))12L else if (isTRUE(p$three)) 9L else 8L
  data.frame(condition=factor(if (n==12L)rep(c("A","B","C"),each=4L)else if(n==9L) rep(c("A","B","C"),each=3L) else rep(c("A","B"),each=4L)),
    exposure=if(n==12L)rep(1:4,3L)else if (n==9L) rep(1:3,3L) else rep(1:4,2L))
}
queries_for <- function(p) {
  contrast <- c(rep(0,p$columns-1L),1)
  q <- function(name, numeric=NULL, character=NULL, cooks="auto", threshold=0)
    list(name=name,contrast=if(is.null(numeric))list()else bits(numeric),character=if(is.null(character))NULL else as.list(character),
      threshold=bits(threshold)[[1L]],alternative="greaterAbs",df=NULL,cooks=cooks,
      filtering=TRUE,alpha=bits(.1)[[1L]],missing=bits(NA_real_)[[1L]])
  result <- list(q("numeric",contrast),q("reverse_unmasked",-contrast,cooks="off"))
  if (is.null(p$reduced)) result <- c(result,list(q("threshold",contrast,threshold=.5)))
  if (!is.null(p$character)) result <- c(result,list(q("character",character=p$character),
    q("character_reverse",character=p$character[c(1L,3L,2L)],cooks="0")))
  result
}
cases <- lapply(names(profiles),function(name) {
  p <- profiles[[name]]; d <- data_for(p)
  y <- if(isTRUE(p$replace))replaced else counts
  if(isTRUE(p$four))y<-cbind(y,y[,5:8])else if (isTRUE(p$three)) y <- cbind(y,y[,6L])
  list(name=name,counts=matrix_descriptor(y),design=NULL,reduced=NULL,
    formula=list(full=p$full,reduced=p$reduced,
      factors=list(list(name="condition",sample_levels=as.list(as.character(d$condition)),levels=as.list(levels(d$condition)),reference=levels(d$condition)[[1L]])),
      numeric=list(list(name="exposure",values=bits(d$exposure)))),
    size_factors=if(isTRUE(p$size))bits(size)else NULL,
    normalization_factors=if(isTRUE(p$factors))matrix_descriptor(factors)else NULL,weights=NULL,
    fit_type="mean",test=if(is.null(p$reduced))"Wald"else"LRT",threads=1L,
    min_replicates=if(isTRUE(p$replace))3L else 0L,
    replacement_cutoff=bits(rsdeseq2:::.rsdeseq2_default_cooks_cutoff(matrix(0,nrow=nrow(d),ncol=p$columns)))[[1L]],
    queries=queries_for(p))
})
input <- file.path(out,"input.json")
jsonlite::write_json(list(fresh_fit_cases=cases),input,auto_unbox=TRUE,null="null",digits=NA)
output <- system2(driver,shQuote(normalizePath(input)),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(output,"status")))stop(paste(output,collapse="\n"))
writeLines(output,file.path(out,"standalone.json"))
native <- jsonlite::fromJSON(paste(output,collapse="\n"),simplifyVector=FALSE)$cases
failures <- character()
records <- lapply(seq_along(cases),function(i) {
  c <- cases[[i]]; p <- profiles[[c$name]]; d <- data_for(p)
  y <- matrix_input(c$counts)
  rownames(y) <- paste0("gene",seq_len(nrow(y))); colnames(y) <- paste0("sample",seq_len(ncol(y)))
  rownames(d) <- colnames(y)
  dds <- tryCatch(DESeq2::DESeqDataSetFromMatrix(y,d,as.formula(p$full)),error=identity)
  fitted <- if(inherits(dds,"error"))dds else tryCatch(blocked(rsdeseq2::DESeq(dds,fitType="mean",test=c$test,
    reduced=if(is.null(p$reduced))NULL else as.formula(p$reduced),
    sizeFactors=if(is.null(c$size_factors))NULL else decode(c$size_factors),
    normalizationFactors=matrix_input(c$normalization_factors),
    minReplicatesForReplace=if(c$min_replicates==0L)Inf else c$min_replicates)),error=identity)
  expected <- native[[i]]
  if (!is.null(p$error)) {
    valid <- inherits(fitted,"error") && !is.null(expected$error) &&
      grepl(p$error,conditionMessage(fitted),ignore.case=TRUE) && grepl(p$error,expected$error,ignore.case=TRUE)
    if(!valid)failures <<- c(failures,paste(c$name,"required rejection differs"))
    return(list(name=c$name,rejected=valid))
  }
  if(inherits(fitted,"error") || !is.null(expected$error)) {
    failures <<- c(failures,paste(c$name,"unexpected failure:",if(inherits(fitted,"error"))conditionMessage(fitted)else expected$error))
    return(list(name=c$name,error=TRUE))
  }
  fit <- blocked(rsdeseq2:::.rsdeseq2_fitted_view(fitted))
  stopifnot(identical(rownames(fit$beta),rownames(y)),nrow(fit$results)==nrow(y))
  actual <- encode(fit); actual$queries <- query_values(fitted,c$queries)
  actual <- jsonlite::fromJSON(jsonlite::toJSON(actual,auto_unbox=TRUE,null="null"),simplifyVector=FALSE)
  jsonlite::write_json(actual,file.path(out,paste0(c$name,"-wrapped.json")),auto_unbox=TRUE,null="null",pretty=TRUE)
  saveRDS(fitted,file.path(out,paste0(c$name,".rds")))
  if(!identical(canonical(actual),canonical(expected$value)))
    failures <<- c(failures,paste(c$name,paste(all.equal(canonical(actual),canonical(expected$value),tolerance=0),collapse="; ")))
  list(name=c$name,value=actual)
})
jsonlite::write_json(list(cases=records),file.path(out,"wrapped.json"),auto_unbox=TRUE,null="null",pretty=TRUE)
jsonlite::write_json(list(cases=length(cases),successful_fits=sum(vapply(records,function(r)!is.null(r$value),logical(1))),rejections=sum(vapply(records,function(r)isTRUE(r$rejected),logical(1))),queries=sum(vapply(records,function(r)length(r$value$queries),integer(1))),failures=as.list(failures),R_calculations_blocked=TRUE),file.path(out,"validation.json"),auto_unbox=TRUE,pretty=TRUE)
if(length(failures))stop(paste(failures,collapse="\n"))
cat(length(cases),"complete formula fits and their later queries match every numerical bit\n")
