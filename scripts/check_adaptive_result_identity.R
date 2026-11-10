#!/usr/bin/env Rscript
# Actual native/R-free outputs, including the wrapper's independent queries.
args <- commandArgs(TRUE)
if(length(args)!=3L) stop('usage: STANDALONE_DRIVER R_LIBRARY FRESH_OUTPUT_DIR')
driver <- normalizePath(args[[1L]])
library(rsdeseq2,lib.loc=normalizePath(args[[2L]]))
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
out <- args[[3L]]
if(dir.exists(out)) stop('Use a fresh output directory')
dir.create(out,recursive=TRUE)
reference <- jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/adaptive_result_r462.json'),simplifyVector=FALSE)
fields <- c('name','estimates','standard_errors','query_estimates','query_standard_errors','prior','alpha','df','threshold','svalue')
cases <- lapply(reference$cases,function(c)c[fields])
input <- file.path(out,'input.json')
jsonlite::write_json(list(adaptive_cases=cases),input,auto_unbox=TRUE,null='null')
output <- system2(driver,shQuote(normalizePath(input)),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(output,'status')))stop(paste(output,collapse='\n'))
standalone <- jsonlite::fromJSON(paste(output,collapse='\n'),simplifyVector=FALSE)
decode <- function(v)vapply(v,function(h){h<-sub('^bits:','',h);readBin(as.raw(strtoi(substring(h,seq(1L,15L,2L),seq(2L,16L,2L)),16L)),'double',n=1L,size=8L,endian='big')},numeric(1L))
blocked <- function(expr) {
  functions <- list(base=c('sqrt','exp','expm1','log','log1p','log2','mean','sum','sort','order','rank','cumsum','^'),stats=c('pnorm','dnorm','pt','dt','integrate'))
  on.exit(for(package in names(functions))for(name in functions[[package]])suppressMessages(untrace(name,where=asNamespace(package))),add=TRUE)
  for(package in names(functions))for(name in functions[[package]])suppressWarnings(suppressMessages(trace(name,where=asNamespace(package),print=FALSE,tracer=quote(stop('adaptive bridge executed an R calculation',call.=FALSE)))))
  force(expr)
}
guard <- tryCatch(blocked(sqrt(2)),error=identity)
stopifnot(inherits(guard,'error'),grepl('executed an R calculation',conditionMessage(guard)))
encode <- function(v){if(is.numeric(v))return(lapply(unname(v),hex));if(is.logical(v)||is.character(v))return(as.list(unname(v)));if(is.list(v))return(lapply(v,encode));v}
wrapped <- lapply(cases,function(c){
  prior <- c$prior
  for(field in c('weights','means','sds','lower','upper'))if(!is.null(prior[[field]]))prior[[field]]<-decode(prior[[field]])
  arguments <- list(betahat=decode(c$estimates),sebetahat=decode(c$standard_errors),g=prior,fixg=TRUE,alpha=decode(c$alpha),method='shrink',mixcompdist=if(is.null(c$df))'normal'else'uniform')
  if(!is.null(c$df))arguments$df<-decode(c$df)
  result <- blocked(rsdeseq2:::.rsdeseq2_formula_request(list(operation='ashr_shrinkage',arguments=arguments,
    query_estimates=decode(c$query_estimates),query_standard_errors=decode(c$query_standard_errors),
    lfc_threshold=decode(c$threshold),svalue=c$svalue)))
  list(name=c$name,value=encode(result[c('posterior_mean','posterior_sd','svalue','warnings')]))
})
canonical <- function(x){if(!is.list(x))return(x);if(!is.null(names(x)))x<-x[sort(names(x))];lapply(x,canonical)}
wrapped <- .formula_canonical(list(cases=wrapped))
jsonlite::write_json(wrapped,file.path(out,'wrapped.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
writeLines(output,file.path(out,'standalone.json'))
stopifnot(identical(canonical(wrapped),canonical(standalone)))
cat(length(cases),'adaptive result batches: every output bit matches; R calculations blocked\n')
