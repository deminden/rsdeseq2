#!/usr/bin/env Rscript
# Exact comparison; numerical corrections require an independent original oracle.
args <- commandArgs(TRUE)
if(!length(args)%in%c(5L,6L)) stop("usage: BASELINE CANDIDATE WORKFLOW STRICT|ORIGIN_V2|ASHR_SQP REPORT [CANDIDATE_LIBRARY]")
a <- readRDS(args[[1L]]); b <- readRDS(args[[2L]])
workflow <- args[[3L]]; mode <- args[[4L]]
if(!mode %in% c("STRICT","ORIGIN_V2","ASHR_SQP")) stop("Unsupported comparison mode")
expectedFormula <- if(workflow %in% c("rlog","vst","replacement","two_group_results")) ~condition else ~batch+condition
expectedOrigin <- list(kind="formula",variables=all.vars(expectedFormula))
migrations <- character()
sqp <- list()
compare <- function(a,b,path="payload") {
  if(inherits(a,"rsdeseq2RustFit")||inherits(b,"rsdeseq2RustFit")) {
    stopifnot(inherits(a,"rsdeseq2RustFit"),inherits(b,"rsdeseq2RustFit"))
    if(mode=="ORIGIN_V2"&&identical(a$schemaVersion,1L)&&is.null(a$designOrigin)&&identical(b$schemaVersion,2L)) {
      # Verify the added metadata against the actual registered workflow and
      # its complete fitted model matrix before normalizing just two fields.
      stopifnot(identical(b$designOrigin,expectedOrigin))
      matrix <- model.matrix(expectedFormula,b$colData)
      stopifnot(identical(dimnames(matrix),dimnames(b$design)),
        identical(dim(matrix),dim(b$design)),identical(as.vector(matrix),as.vector(b$design)))
      migrations <<- c(migrations,path)
      b$schemaVersion <- 1L;b$designOrigin <- NULL
    }
    stopifnot(isTRUE(all.equal(a,b,tolerance=0,check.attributes=TRUE)))
    return(invisible(NULL))
  }
  # Only descend through ordinary payload/metadata lists. Numerical arrays,
  # tables, attributes, nulls and other classes retain exact comparison.
  if(is.list(a)&&is.list(b)&&!is.object(a)&&!is.object(b)) {
    stopifnot(identical(attributes(a),attributes(b)),length(a)==length(b))
    for(i in seq_along(a)) compare(a[[i]],b[[i]],paste0(path,"/",if(is.null(names(a)))i else names(a)[i]))
  } else stopifnot(isTRUE(all.equal(a,b,tolerance=0,check.attributes=TRUE)))
}
if(mode=="ASHR_SQP" && workflow=="shrink_ashr") {
  if(length(args)!=6L) stop("ASHR_SQP requires the candidate library")
  # Establish unchanged inputs/fit before using its unshrunk estimates. Never
  # normalize fitted-state, metadata, p-values or unrelated numerical columns.
  compare(a$counts,b$counts,"payload/counts")
  compare(a$fitted,b$fitted,"payload/fitted")
  stopifnot(is.data.frame(a$output),is.data.frame(b$output),
    identical(attributes(a$output),attributes(b$output)))
  suppressPackageStartupMessages(library(rsdeseq2,lib.loc=normalizePath(args[[6L]])))
  stopifnot(getRversion()=="4.6.2",as.character(packageVersion("ashr"))=="2.2.63")
  raw <- rsdeseq2::results(b$fitted,name=tail(b$fitted$resultsNames,1L))
  expected <- ashr::ash(raw$log2FoldChange,raw$lfcSE,mixcompdist="normal",method="shrink")$result
  for(field in c("log2FoldChange","lfcSE")) {
    target <- expected[[if(field=="log2FoldChange")"PosteriorMean"else"PosteriorSD"]]
    before <- a$output[[field]]; after <- b$output[[field]]
    stopifnot(is.numeric(before),is.numeric(after),length(after)==length(target),
      length(before)==length(target),identical(is.na(after),is.na(target)),
      !any(is.infinite(after)),!any(is.infinite(target)))
    finite <- is.finite(target)
    stopifnot(any(finite))
    error <- max(abs(after[finite]-target[finite]))
    # Reuse the existing standalone SQP oracle bound; no broad all.equal
    # tolerance is introduced for the workflow or its fitted arrays.
    stopifnot(is.finite(error),error<=2e-6)
    sqp[[length(sqp)+1L]] <- data.frame(path=paste0("payload/output/",field),
      migration="original_ashr_SQP_checked",bound=2e-6,candidate_original_error=error,
      baseline_original_error=if(identical(is.na(before),is.na(target))&&!any(is.infinite(before)))max(abs(before[finite]-target[finite]))else NA_real_,
      candidate_baseline_difference=max(abs(after[finite]-before[finite])))
    b$output[[field]] <- before
  }
}
compare(a,b)
report <- if(length(sqp))do.call(rbind,sqp)else
  data.frame(path=migrations,migration=rep("schema1_to_2_validated_formula_origin",length(migrations)))
write.table(report,
  args[[5L]],sep="\t",quote=FALSE,row.names=FALSE)
cat("Payload comparison passed; validated origin migrations:",length(migrations),
  "; original-ashr posterior columns checked:",length(sqp),"\n")
