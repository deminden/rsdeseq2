#!/usr/bin/env Rscript
# Untimed fit reuse; measure public result calls and registered-native segments.
args <- commandArgs(TRUE)
if(length(args)<3L) stop('usage: R_LIBRARY FRESH_OUTPUT_TSV FITTED_WORKFLOW_RDS...')
if(file.exists(args[[2L]])) stop('Use a fresh output path')
library(rsdeseq2,lib.loc=normalizePath(args[[1L]]))
stopifnot(getRversion()=='4.6.2')
rows <- list()
namespace <- asNamespace('rsdeseq2')
native <- Filter(function(name) {
  f <- get(name,envir=namespace,inherits=FALSE)
  typeof(f)=='closure' && '.Call'%in%all.names(body(f))
},ls(namespace,all.names=TRUE))
stopifnot(length(native)>0L)
for(path in args[-c(1L,2L)]) {
  payload <- readRDS(path)
  fit <- payload$fitted
  stopifnot(inherits(fit,'rsdeseq2RustFit'),length(fit$resultsNames)>1L)
  requests <- lapply(seq_len(10L),function(i)
    list(name=fit$resultsNames[[1L+i%%ncol(fit$beta)]]))
  workloads <- list(
    one=function()rsdeseq2::results(fit,name=tail(fit$resultsNames,1L)),
    ten_single=function()lapply(requests,function(r)do.call(rsdeseq2::results,c(list(object=fit),r))),
    ten_batch=function()rsdeseq2::resultsList(fit,requests))
  # These untraced values establish that profiling does not change outputs.
  reference <- lapply(workloads,function(f)f())
  for(name in names(workloads)) {
    .rs_boundary_profile <- new.env(parent=emptyenv())
    .rs_boundary_profile$elapsed <- 0
    .rs_boundary_profile$calls <- 0L
    .rs_boundary_profile$symbols <- character()
    .rs_boundary_profile$depth <- 0L
    for(symbol in native) suppressMessages(trace(symbol,where=namespace,print=FALSE,
      tracer=substitute({
        .rs_state <- get('.rs_boundary_profile',envir=.GlobalEnv)
        if(.rs_state$depth==0L) .rs_state$start <- proc.time()[['elapsed']]
        .rs_state$depth <- .rs_state$depth+1L
        .rs_state$calls <- .rs_state$calls+1L
        .rs_state$symbols <- c(.rs_state$symbols,SYMBOL)
      },list(SYMBOL=symbol)),
      exit=quote({
        .rs_state <- get('.rs_boundary_profile',envir=.GlobalEnv)
        .rs_state$depth <- .rs_state$depth-1L
        if(.rs_state$depth==0L)
          .rs_state$elapsed <- .rs_state$elapsed+(proc.time()[['elapsed']]-.rs_state$start)
      })))
    start <- proc.time()[['elapsed']]
    value <- tryCatch({
      result <- workloads[[name]]()
      elapsed <- proc.time()[['elapsed']]-start
      result
    },finally=for(symbol in native)suppressMessages(untrace(symbol,where=namespace)))
    stopifnot(identical(value,reference[[name]]),.rs_boundary_profile$calls>0L,
      .rs_boundary_profile$depth==0L,
      elapsed>=.rs_boundary_profile$elapsed-0.003)
    rows[[length(rows)+1L]] <- data.frame(input=basename(path),genes=nrow(fit$beta),
      coefficients=ncol(fit$beta),workload=name,queries=if(name=='one')1L else 10L,
      total_s=elapsed,native_boundary_s=.rs_boundary_profile$elapsed,
      outside_native_boundary_s=max(0,elapsed-.rs_boundary_profile$elapsed),
      native_calls=.rs_boundary_profile$calls,
      native_symbols=paste(unique(.rs_boundary_profile$symbols),collapse=','))
  }
}
write.table(do.call(rbind,rows),args[[2L]],sep='\t',quote=FALSE,row.names=FALSE)
cat('All profiled public outputs are identical to their untraced controls.\n',
  'Native segments cover outermost R adapters containing .Call, including their validation/marshalling and Rust computation.\n',
  'Remaining time includes R presentation/validation and tracing overhead; nested adapters are not double-counted.\n',
  'Elapsed timers have millisecond resolution; no isolated speed improvement is inferred.\n',sep='')
