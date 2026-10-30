#!/usr/bin/env Rscript
# Exact workflow comparison, with one explicitly enabled provenance migration.
args <- commandArgs(TRUE)
if(length(args)!=5L) stop("usage: BASELINE CANDIDATE WORKFLOW STRICT|ORIGIN_V2 REPORT")
a <- readRDS(args[[1L]]); b <- readRDS(args[[2L]])
workflow <- args[[3L]]; mode <- args[[4L]]
if(!mode %in% c("STRICT","ORIGIN_V2")) stop("Unsupported comparison mode")
expectedFormula <- if(workflow %in% c("rlog","vst","replacement","two_group_results")) ~condition else ~batch+condition
expectedOrigin <- list(kind="formula",variables=all.vars(expectedFormula))
migrations <- character()
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
compare(a,b)
write.table(data.frame(path=migrations,migration=rep("schema1_to_2_validated_formula_origin",length(migrations))),
  args[[5L]],sep="\t",quote=FALSE,row.names=FALSE)
cat("Exact payload comparison passed; validated origin migrations:",length(migrations),"\n")
