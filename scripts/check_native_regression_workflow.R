#!/usr/bin/env Rscript
# Focused complete bridge/R-free checks and public general/nbinomR comparisons.
args <- commandArgs(TRUE)
if (length(args) != 4L) stop('usage: DRIVER R_LIBRARY FRESH_OUTPUT_DIR PAIRED_FITS_RDS')
driver <- normalizePath(args[[1L]])
library(rsdeseq2, lib.loc=normalizePath(args[[2L]]))
requireNamespace('DESeq2'); requireNamespace('apeglm')
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub('~+~', ' ', script, fixed=TRUE)), '..'))
source(file.path(root, 'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
source(file.path(root, 'r-pkg/rsdeseq2/tests/testthat/helper-fit-data.R'))
out <- args[[3L]]
if (dir.exists(out)) stop('Use a fresh output directory')
dir.create(out, recursive=TRUE)
decode <- function(v) vapply(v, function(h) {
  h <- sub('^bits:', '', h)
  readBin(as.raw(strtoi(substring(h,seq(1L,15L,2L),seq(2L,16L,2L)),16L)), 'double', n=1L, size=8L, endian='big')
}, numeric(1L))
encode <- function(v) lapply(as.numeric(v), hex)
encode_optional <- function(v) lapply(as.numeric(v), function(x) if (is.na(x)) NULL else hex(x))
matrix_input <- function(m) if (is.null(m)) NULL else matrix(decode(m$values), nrow=m$rows, ncol=m$columns)
compare <- function(a,b,bound) {
  stopifnot(length(a)==length(b), identical(unname(is.na(a)),unname(is.na(b))))
  finite <- is.finite(a)&is.finite(b)
  stopifnot(all(finite|is.na(a)), all(finite|is.na(b)))
  error <- if (any(finite)) max(abs(a[finite]-b[finite])) else 0
  stopifnot(error<=bound)
  error
}
blocked <- function(expr) {
  functions <- list(base=c('sqrt','exp','expm1','log','log1p','log2','mean','rank','cumsum','^','%*%','lgamma'),
    stats=c('pnorm','dnorm','dnbinom','dbinom','dt','integrate','quantile','optim','optimHess','uniroot'),
    apeglm=c('apeglm','priorVar','svalue'))
  on.exit(for (package in names(functions)) for (name in functions[[package]])
    suppressMessages(untrace(name, where=asNamespace(package))), add=TRUE)
  for (package in names(functions)) for (name in functions[[package]])
    suppressWarnings(suppressMessages(trace(name,where=asNamespace(package),print=FALSE,
      tracer=substitute(stop(MESSAGE,call.=FALSE),list(MESSAGE=paste('native regression executed R math:',package,name))))))
  force(expr)
}
guard <- tryCatch(blocked(stats::dnorm(0)),error=identity)
stopifnot(inherits(guard,'error'),grepl('executed R math',conditionMessage(guard)))
reference <- jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/regression_posterior_r462.json'), simplifyVector=FALSE)
input <- file.path(out,'input.json')
fields <- c('name','observations','design','parameters','weights','offsets','mle','coefficient','nuisance',
  'scale','nuisance_scale','prior_df','prior_mean','nuisance_mean','method','optimizer','bounds','log_link','threshold','likelihood')
jsonlite::write_json(list(regression_cases=lapply(reference$cases,function(c)c[fields]),density_cases=reference$leaves),input,auto_unbox=TRUE,null='null')
output <- system2(driver,shQuote(normalizePath(input)),stdout=TRUE,stderr=TRUE)
if (!is.null(attr(output,'status'))) stop(paste(output,collapse='\n'))
writeLines(output,file.path(out,'standalone.json'))
standalone <- jsonlite::fromJSON(paste(output,collapse='\n'),simplifyVector=FALSE)
wrapped <- lapply(seq_along(reference$cases),function(i) {
  c <- reference$cases[[i]]
  model <- c$likelihood
  if (!is.null(model)) {
    model$arguments <- unlist(model$arguments)
    model$expression <- .formula_ast_decode(model$expression)
    model$constants <- lapply(model$constants,.formula_decode)
  }
  arguments <- list(Y=matrix_input(c$observations),x=matrix_input(c$design),param=matrix_input(c$parameters),
    coef=c$coefficient+1L,weights=matrix_input(c$weights),offset=matrix_input(c$offsets),mle=matrix_input(c$mle),
    prior.control=list(no.shrink=as.numeric(unlist(c$nuisance))+1L,prior.mean=decode(c$prior_mean),
      prior.scale=decode(c$scale),prior.df=decode(c$prior_df),prior.no.shrink.mean=decode(c$nuisance_mean),
      prior.no.shrink.scale=decode(c$nuisance_scale)),threshold=decode(c$threshold),
    method=c$method,optim.method=c$optimizer,bounds=decode(c$bounds),log.link=c$log_link,
    log.lik=if(is.null(model))list(type='negative_binomial')else list(type='restricted',model=model))
  r <- tryCatch(blocked(rsdeseq2:::.rsdeseq2_formula_request(list(operation='apeglm_shrinkage',arguments=arguments))), error=identity)
  s <- standalone$cases[[i]]
  if (!is.null(c$expected_error)) {
    stopifnot(inherits(r,'error'),!is.null(s$error))
    return(list(name=c$name,rejected=TRUE))
  }
  if (inherits(r,'error')) stop(c$name,': ',conditionMessage(r))
  stopifnot(is.null(s$error))
  s <- s$value
  expected <- list(mode=s$summary_mode,sd=s$summary_sd,svalue=s$sign_svalues,
    threshold_svalue=s$threshold_svalues,prior_scale=s$prior_scale,prior_variance=s$prior_variance,
    nuisance=lapply(s$nuisance,function(v)v+1L),convergence=lapply(s$diagnostics,`[[`,'convergence'))
  actual <- list(mode=encode(r$posterior_mode),sd=encode_optional(r$posterior_sd),svalue=encode_optional(r$svalue),
    threshold_svalue=encode_optional(r$threshold_svalue),prior_scale=encode(r$prior_control$prior.scale),
    prior_variance=if(is.null(r$prior_control$prior.var))NULL else encode(r$prior_control$prior.var),
    nuisance=as.list(as.integer(r$prior_control$no.shrink)),
    convergence=lapply(r$convergence,function(v)if(is.na(v))NULL else as.integer(v)))
  stopifnot(identical(.formula_canonical(actual),.formula_canonical(expected)))
  list(name=c$name,value=actual)
})
jsonlite::write_json(wrapped,file.path(out,'wrapped.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
# Reuse the accepted count fit. Only the newly changed posterior paths are fitted.
fixture <- .rsdeseq2_large_transform_fixture()
cached <- readRDS(args[[4L]])
rust <- cached$rust
stopifnot(isTRUE(all.equal(as.matrix(SummarizedExperiment::assay(rust,'counts')),fixture$counts,
  check.attributes=TRUE,tolerance=0)))
settings <- list(general=list(apeMethod='general'),nbinomR=list(apeMethod='nbinomR'))
if (all(names(settings) %in% names(cached$reference))) {
  original <- cached$reference[names(settings)]
} else {
  original_fit <- DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(fixture$counts,fixture$colData,~condition),fitType='mean',quiet=TRUE)
  coefficient <- tail(DESeq2::resultsNames(original_fit),1L)
  original <- lapply(settings,function(options) do.call(DESeq2::lfcShrink,
    c(list(dds=original_fit,coef=coefficient,type='apeglm',lfcThreshold=.5,quiet=TRUE),options)))
}
# The custom density expression goes through the public closure importer, then
# runs in Rust. Its body is never evaluated by R during a candidate call.
custom <- function(y,x,beta,param,offset)stats::dnbinom(y,mu=exp(x%*%beta+offset),size=1/param,log=TRUE)
settings$custom <- list(apeMethod='general',log.lik=custom)
actual <- lapply(settings,function(options)blocked(do.call(rsdeseq2::lfcShrink,
  c(list(dds=rust,coef='conditionB',type='apeglm',lfcThreshold=.5,quiet=TRUE),options))))
errors <- lapply(names(settings),function(name) {
  a <- actual[[name]]; b <- original[[if(name=='custom')'general'else name]]
  stopifnot(identical(DESeq2::priorInfo(a)$engine,'Rust'))
  list(profile=name,mode=compare(a$log2FoldChange,b$log2FoldChange,.01),
    sd=compare(a$lfcSE,b$lfcSE,.01),threshold_svalue=compare(a$svalue,b$svalue,.05))
})
saveRDS(list(rust=rust,reference=original,actual=actual),file.path(out,'paired_fits.rds'))
jsonlite::write_json(list(R=as.character(getRversion()),DESeq2=as.character(packageVersion('DESeq2')),
  apeglm=as.character(packageVersion('apeglm')),genes=nrow(fixture$counts),complete_native_standalone_batches=length(wrapped),
  R_calculations_blocked=TRUE,existing_bounds_unchanged=TRUE,public_errors=errors),
  file.path(out,'validation.json'),auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
cat(length(wrapped),'general/analytic NB batches match standalone bits or reject together; public general/nbinomR/custom checks pass with R math blocked\n')
