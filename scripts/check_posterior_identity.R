#!/usr/bin/env Rscript
# Test-only bit comparison of standalone and native-bridge posterior primitives.
args <- commandArgs(TRUE)
if (length(args) != 3L) stop('usage: STANDALONE_DRIVER R_LIBRARY FRESH_OUTPUT_DIR')
driver <- normalizePath(args[[1L]])
library(rsdeseq2,lib.loc=normalizePath(args[[2L]]))
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
out <- args[[3L]]
if (dir.exists(out)) stop('Use a fresh identity output directory')
dir.create(out,recursive=TRUE)
reference <- jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/normal_mixture_posterior_r462.json'),simplifyVector=FALSE)
input <- list(posterior_cases=reference$cases,svalue_cases=reference$svalues)
path <- file.path(out,'input.json')
jsonlite::write_json(input,path,auto_unbox=TRUE,null='null',digits=NA)
output <- system2(driver,shQuote(normalizePath(path)),stdout=TRUE,stderr=TRUE)
if (!is.null(attr(output,'status'))) stop(paste(output,collapse='\n'))
standalone <- jsonlite::fromJSON(paste(output,collapse='\n'),simplifyVector=FALSE)
number <- function(x) {
  x <- unlist(x,use.names=FALSE)
  x[x=='NA'] <- NA_character_
  as.numeric(x)
}
blocked <- function(expression) {
  functions <- list(base=c('sqrt','exp','log','mean','cumsum','sort','rank'),stats=c('pnorm','dnorm'))
  on.exit(for (package in names(functions)) for (name in functions[[package]])
    suppressMessages(untrace(name,where=asNamespace(package))),add=TRUE)
  for (package in names(functions)) for (name in functions[[package]])
    suppressWarnings(suppressMessages(trace(name,where=asNamespace(package),print=FALSE,
      tracer=quote(stop('posterior bridge executed an R calculation',call.=FALSE)))))
  force(expression)
}
guard <- tryCatch(blocked(sqrt(2)),error=identity)
stopifnot(inherits(guard,'error'),grepl('executed an R calculation',conditionMessage(guard)))
cases <- lapply(reference$cases,function(case) {
  fields <- c('weights','means','sds','estimates','standard_errors','thresholds')
  request <- c(list(operation='posterior_normal_mixture'),lapply(case[fields],number))
  result <- blocked(rsdeseq2:::.rsdeseq2_formula_request(request))
  c(list(name=case$name),lapply(result,function(v)lapply(unname(v),hex)))
})
sv <- lapply(reference$svalues,function(case) {
  request <- list(operation='posterior_svalues',rates=number(case$rates))
  result <- blocked(rsdeseq2:::.rsdeseq2_formula_request(request))
  list(name=case$name,values=lapply(unname(result),function(v)if(is.na(v))NULL else hex(v)))
})
canonical <- function(x) {
  if(!is.list(x))return(x)
  if(!is.null(names(x)))x<-x[sort(names(x))]
  lapply(x,canonical)
}
wrapped <- .formula_canonical(list(cases=cases,svalues=sv))
stopifnot(identical(canonical(wrapped),canonical(standalone)))
jsonlite::write_json(wrapped,file.path(out,'wrapped.json'),auto_unbox=TRUE,null='null',pretty=TRUE)
writeLines(output,file.path(out,'standalone.json'))
cat('520 normal-mixture posterior observations, both tails at ten thresholds, and 8 s-value panels: every output bit identical; R calculations blocked\n')
