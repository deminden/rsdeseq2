#!/usr/bin/env Rscript
# Compare actual standalone and wrapped factor paths with original-R fixtures.
args <- commandArgs(TRUE)
if(length(args)!=3L)stop("usage: STANDALONE_DRIVER R_LIBRARY FRESH_OUTPUT_DIR")
driver <- normalizePath(args[[1L]])
library(rsdeseq2,lib.loc=normalizePath(args[[2L]]))
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~"," ",script,fixed=TRUE)),".."))
source(file.path(root,"r-pkg/rsdeseq2/tests/testthat/helper-reference-values.R"))
source(file.path(root,"r-pkg/rsdeseq2/tests/testthat/helper-formula.R"))
canonical <- function(value) {
  if(!is.list(value))return(value)
  if(!is.null(names(value)))value<-value[sort(names(value))]
  lapply(value,canonical)
}
out <- args[[3L]]
if(dir.exists(out))stop("Use a fresh identity output directory")
dir.create(out,recursive=TRUE)
reference <- jsonlite::fromJSON(.formula_reference_path(file.path(root,"crates/rsdeseq2/tests/data"),"native_factor_r462.json"),simplifyVector=FALSE)
path <- file.path(out,"input.json")
jsonlite::write_json(reference,path,auto_unbox=TRUE,null="null",digits=NA)
output <- system2(driver,shQuote(normalizePath(path)),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(output,"status")))stop(paste(output,collapse="\n"))
writeLines(output,file.path(out,"standalone.json"))
standalone <- jsonlite::fromJSON(paste(output,collapse="\n"),simplifyVector=FALSE)
decode <- function(hex)readBin(as.raw(strtoi(substring(hex,seq(1,15,2),seq(2,16,2)),16L)),"double",1L,size=8,endian="big")
bits <- function(matrix)lapply(unname(as.vector(t(matrix))),hex)
d <- as.data.frame(stats::setNames(lapply(reference$numeric,function(v)vapply(v$values,decode,0.0)),vapply(reference$numeric,`[[`,"","name")))
for(v in reference$factors)d[[v$name]] <- factor(unlist(v$sample_levels),levels=unlist(v$levels))
options(na.action="na.pass")
wrapped <- lapply(seq_along(reference$native_formula_cases),function(i) {
  case <- reference$native_formula_cases[[i]]
  result <- tryCatch({
    formula <- as.formula(case$formula)
    frame <- rsdeseq2:::.rsdeseq2_model_frame(formula,d)
    matrix <- rsdeseq2:::.rsdeseq2_model_matrix(formula,frame)
    if(any(!is.finite(matrix)))stop("non-finite fitting matrix")
    list(formula=formula,frame=frame,matrix=matrix)
  },error=identity)
  native <- standalone$cases[[i]]
  valid <- !inherits(result,"error")
  stopifnot(identical(native$name,case$name),identical(native$ok,valid),identical(valid,!is.null(case$expected)))
  if(!valid)return(list(name=case$name,ok=FALSE))
  factors <- result$frame[vapply(result$frame,is.factor,TRUE)]
  frames <- lapply(unname(factors),function(value)list(levels=as.list(levels(value)),reference=levels(value)[[1L]],sample_levels=as.list(as.character(value)),ordered=is.ordered(value)))
  nativeFrames <- lapply(native$frames,function(value){value$name<-NULL;value})
  frames <- canonical(frames);nativeFrames <- canonical(nativeFrames)
  if(!identical(frames,nativeFrames))stop(case$name,": factor metadata: ",paste(all.equal(frames,nativeFrames),collapse="; "))
  for(name in names(factors))if(!is.ordered(result$frame[[name]])) {
    identity <- diag(nlevels(result$frame[[name]]))
    dimnames(identity) <- list(levels(result$frame[[name]]),levels(result$frame[[name]]))
    attr(result$frame[[name]],"contrasts") <- identity
  }
  expanded <- rsdeseq2:::.rsdeseq2_model_matrix(result$formula,result$frame)
  values <- bits(result$matrix[,unlist(case$r_columns),drop=FALSE])
  expandedValues <- bits(expanded[,unlist(case$r_expanded_columns),drop=FALSE])
  stopifnot(identical(values,native$values),identical(values,case$expected),
    identical(expandedValues,native$expanded_values),identical(expandedValues,case$expanded_expected),
    identical(frames,nativeFrames),identical(canonical(native$frames),canonical(case$frames)),
    identical(native$names,case$native_names),identical(native$expanded_names,case$expanded_names),
    identical(native$groups,case$groups))
  list(name=case$name,ok=TRUE,r_columns=case$r_columns,values=values,
    r_expanded_columns=case$r_expanded_columns,expanded_values=expandedValues,frames=frames)
})
expansions <- lapply(seq_along(reference$factor_expansions),function(i) {
  case <- reference$factor_expansions[[i]]
  result <- tryCatch(rsdeseq2:::.rsdeseq2_formula_request(list(operation="expanded_factors",values=list(.formula_decode(case$input)),null_level=case$null_level,collation="C"))[[1L]],error=identity)
  valid <- !inherits(result,"error")
  native <- standalone$expansions[[i]]
  stopifnot(identical(valid,native$ok),identical(valid,!is.null(case$expected)))
  if(!valid)return(list(name=case$name,ok=FALSE))
  value <- .formula_canonical(.formula_descriptor(result))
  stopifnot(identical(canonical(value),canonical(.formula_canonical(case$expected))),identical(canonical(value),canonical(native$value)))
  list(name=case$name,ok=TRUE,value=value)
})
jsonlite::write_json(list(cases=wrapped,expansions=expansions),file.path(out,"wrapped.json"),auto_unbox=TRUE,null="null",pretty=TRUE)
cat(length(wrapped),"factor matrix/frame cases and",length(expansions),"expanded-factor cases agree across original R, standalone Rust and wrapped Rust after explicit column mapping\n")
