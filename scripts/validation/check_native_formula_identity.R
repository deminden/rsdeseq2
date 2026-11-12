#!/usr/bin/env Rscript
# Compare the real native string parser with the installed Rust-backed R path.
args <- commandArgs(TRUE)
if(length(args)!=3L)stop("usage: STANDALONE_DRIVER R_LIBRARY FRESH_OUTPUT_DIR")
driver <- normalizePath(args[[1L]])
library(rsdeseq2,lib.loc=normalizePath(args[[2L]]))
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~"," ",script,fixed=TRUE)),"..", ".."))
source(file.path(root,"r-pkg/rsdeseq2/tests/testthat/helper-reference-values.R"))
out <- args[[3L]]
if(dir.exists(out))stop("Use a fresh identity output directory")
dir.create(out,recursive=TRUE)
reference <- jsonlite::fromJSON(.formula_reference_path(file.path(root,"crates/rsdeseq2/tests/data"),"native_formula_r462.json"),simplifyVector=FALSE)
path <- file.path(out,"input.json")
jsonlite::write_json(reference,path,auto_unbox=TRUE,null="null",digits=NA)
output <- system2(driver,shQuote(normalizePath(path)),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(output,"status")))stop(paste(output,collapse="\n"))
writeLines(output,file.path(out,"standalone.json"))
standalone <- jsonlite::fromJSON(paste(output,collapse="\n"),simplifyVector=FALSE)
decode <- function(hex)readBin(as.raw(strtoi(substring(hex,seq(1,15,2),seq(2,16,2)),16L)),"double",n=1L,size=8,endian="big")
hexes <- function(v)lapply(unname(v),function(x)paste(format(writeBin(as.double(x),raw(),size=8,endian="big")),collapse=""))
d <- as.data.frame(stats::setNames(lapply(reference$numeric,function(v)vapply(v$values,decode,0.0)),vapply(reference$numeric,`[[`,"","name")))
for(v in reference$factors)d[[v$name]] <- factor(unlist(v$sample_levels),levels=unlist(v$levels))
options(na.action="na.pass")
wrapped <- lapply(reference$native_formula_cases,function(case) {
  result <- tryCatch({
    formula <- as.formula(case$formula)
    frame <- rsdeseq2:::.rsdeseq2_model_frame(formula,d)
    matrix <- rsdeseq2:::.rsdeseq2_model_matrix(formula,frame)
    offsets <- rsdeseq2:::.rsdeseq2_model_offsets(frame)
    if(any(!is.finite(matrix))||any(!is.finite(offsets)))stop("non-finite fitting input")
    list(matrix=matrix[,unlist(case$r_columns),drop=FALSE],offsets=offsets)
  },error=identity)
  if(inherits(result,"error"))return(list(name=case$name,ok=FALSE))
  # Public compatibility labels/order are translated explicitly; the numerical
  # matrix and offsets are compared bit for bit, including signed zero.
  list(name=case$name,ok=TRUE,names=case$native_names,
    values=hexes(as.vector(t(result$matrix))),offsets=hexes(result$offsets))
})
jsonlite::write_json(list(cases=wrapped),file.path(out,"wrapped.json"),auto_unbox=TRUE,null="null",pretty=TRUE)
for(i in seq_along(wrapped)) {
  actual <- wrapped[[i]];native <- standalone$cases[[i]];referenceCase <- reference$native_formula_cases[[i]]
  stopifnot(identical(actual$name,native$name),identical(actual$ok,native$ok),identical(actual$ok,!is.null(referenceCase$expected)))
  if(actual$ok) {
    for(field in c("names","values","offsets")) {
      if(!identical(actual[[field]],native[[field]]))stop(actual$name,": ",field," cross-interface difference")
    }
    stopifnot(identical(actual$values,referenceCase$expected),identical(actual$offsets,referenceCase$offsets))
  }
}
cat(length(wrapped),"actual native-parser/R-wrapper cases match original R bits after explicit column mapping\n")
