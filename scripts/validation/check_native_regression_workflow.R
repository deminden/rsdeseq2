#!/usr/bin/env Rscript
# Focused complete bridge/R-free checks and public general/nbinomR comparisons.
args <- commandArgs(TRUE)
if (!length(args) %in% c(4L,5L,6L)) stop('usage: DRIVER R_LIBRARY FRESH_OUTPUT_DIR PAIRED_FITS_RDS [REFERENCE_JSON [identity-only]]')
driver <- normalizePath(args[[1L]])
library(rsdeseq2, lib.loc=normalizePath(args[[2L]]))
requireNamespace('DESeq2'); requireNamespace('apeglm')
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub('~+~', ' ', script, fixed=TRUE)), '..', '..'))
source(file.path(root, 'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
source(file.path(root, 'r-pkg/rsdeseq2/tests/testthat/helper-fit-data.R'))
out <- args[[3L]]
if (dir.exists(out)) stop('Use a fresh output directory')
dir.create(out, recursive=TRUE)
source(file.path(root,'scripts/validation/native_posterior_check_helpers.R'))
reference_path <- if(length(args)>=5L)args[[5L]]else file.path(root,'crates/rsdeseq2/tests/data/regression_posterior_r462.json')
reference <- jsonlite::fromJSON(reference_path,simplifyVector=FALSE)
random_mode <- !is.null(reference$rng_states)
decode_state <- function(v)vapply(v,function(h)readBin(as.raw(strtoi(substring(h,seq(1L,7L,2L),seq(2L,8L,2L)),16L)), 'integer',n=1L,size=4L,endian='big'),integer(1L))
encode_state <- function(v)lapply(v,function(x)paste(format(writeBin(as.integer(x),raw(),size=4L,endian='big')),collapse=''))
if(random_mode) {
  expand <- function(c,nb=FALSE) {
    r <- reference$rng_states[[c$rng_key]]
    c$rng_state <- r$before
    c$expected_random_state <- r$after[[as.character(c$design$columns)]]
    if(nb) {
      c$observations <- c$counts
      c$parameters <- list(rows=c$counts$rows,columns=1L,values=c$dispersions)
      c$prior_df <- 'bits:3ff0000000000000';c$prior_mean <- c$nuisance_mean <- 'bits:0000000000000000'
      c$method <- 'nbinomC*';c$optimizer <- 'BFGS';c$bounds <- list('bits:fff0000000000000','bits:7ff0000000000000');c$log_link <- TRUE
    }
    c
  }
  reference$cases <- c(lapply(reference$heavy_tail_cases,expand,nb=TRUE),lapply(reference$regression_cases,expand))
}
beta_binomial <- identical(reference$cases[[1L]]$method,'betabinR')
input <- file.path(out,'input.json')
fields <- c('name','observations','design','parameters','weights','offsets','mle','coefficient','nuisance',
  'scale','nuisance_scale','prior_df','prior_mean','nuisance_mean','method','optimizer','bounds','log_link','threshold','likelihood','rng_state','flip_sign')
jsonlite::write_json(list(regression_cases=lapply(reference$cases,function(c)c[fields]),density_cases=if(is.null(reference$leaves))list()else reference$leaves),input,auto_unbox=TRUE,null='null')
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
    flip.sign=if(is.null(c$flip_sign))TRUE else c$flip_sign,
    log.lik=if(is.null(model))list(type='negative_binomial')else list(type='restricted',model=model))
  r <- tryCatch(blocked(rsdeseq2:::.rsdeseq2_formula_request(list(operation='apeglm_shrinkage',arguments=arguments,rng_state=if(is.null(c$rng_state))NULL else decode_state(c$rng_state)))), error=identity)
  s <- standalone$cases[[i]]
  if (!is.null(c$expected_error)) {
    stopifnot(inherits(r,'error'),!is.null(s$error))
    return(list(name=c$name,rejected=TRUE))
  }
  if (inherits(r,'error')) stop(c$name,': ',conditionMessage(r))
  stopifnot(is.null(s$error))
  s <- s$value
  expected <- list(random_state=s$random_state,warnings=s$warnings,mode=s$summary_mode,sd=s$summary_sd,svalue=s$sign_svalues,
    threshold_svalue=s$threshold_svalues,prior_scale=s$prior_scale,prior_variance=s$prior_variance,
    nuisance=lapply(s$nuisance,function(v)v+1L),convergence=lapply(s$diagnostics,`[[`,'convergence'))
  actual <- list(random_state=if(is.null(r$random_state))NULL else encode_state(r$random_state),warnings=as.list(r$warnings),mode=encode(r$posterior_mode),sd=encode_optional(r$posterior_sd),svalue=encode_optional(r$svalue),
    threshold_svalue=encode_optional(r$threshold_svalue),prior_scale=encode(r$prior_control$prior.scale),
    prior_variance=if(is.null(r$prior_control$prior.var))NULL else encode(r$prior_control$prior.var),
    nuisance=as.list(as.integer(r$prior_control$no.shrink)),
    convergence=lapply(r$convergence,function(v)if(is.na(v))NULL else as.integer(v)))
  stopifnot(identical(.formula_canonical(actual),.formula_canonical(expected)))
  if(!is.null(c$expected_random_state))stopifnot(identical(actual$random_state,c$expected_random_state))
  list(name=c$name,value=actual,original_warning_messages_identical=identical(as.list(r$warnings),c$warnings),
    original_warning_events_identical=identical(sub(',.*','',r$warnings),sub(',.*','',unlist(c$warnings))))
})
jsonlite::write_json(wrapped,file.path(out,'wrapped.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
if(length(args)==6L){
  stopifnot(identical(args[[6L]],'identity-only'))
  jsonlite::write_json(list(complete_native_standalone_batches=length(reference$cases),R_calculations_blocked=TRUE,
    native_wrapped_bits='exact',source_warning_messages=lapply(wrapped,`[[`,'original_warning_messages_identical')),
    file.path(out,'validation.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
  cat(length(reference$cases),'complete regression batches have exact wrapped/native bits; R calculations blocked\n')
  quit(status=0)
}
# Reuse the accepted count fit. Only the newly changed posterior paths are fitted.
fixture <- .rsdeseq2_large_transform_fixture()
cached <- readRDS(args[[4L]])
rust <- cached$rust
stopifnot(isTRUE(all.equal(as.matrix(SummarizedExperiment::assay(rust,'counts')),fixture$counts,
  check.attributes=TRUE,tolerance=0)))
if (random_mode) {
  analysis_dds <- rust[seq_len(256L),]
  fit <- rsdeseq2:::.rsdeseq2_fitted_view(analysis_dds)
  raw <- rsdeseq2::results(analysis_dds,name='conditionB')
  observations <- as.matrix(SummarizedExperiment::assay(analysis_dds,'counts'))
  parameters <- cbind(8,observations+20)
  settings <- list(); public_states <- list()
  for(kind in c('Mersenne-Twister',"L'Ecuyer-CMRG"))for(method in c('nbinomC*','betabinC*')) {
    name <- paste(kind,method,sep='_')
    settings[[name]] <- list(apeMethod=method)
    if(method=='betabinC*')settings[[name]]$param<-parameters
    public_states[[name]]<-decode_state(reference$rng_states[[paste(kind,'42',sep='_')]]$before)
  }
  original <- if(all(names(settings)%in%names(cached$reference)))cached$reference[names(settings)]else lapply(names(settings),function(name) {
    assign('.Random.seed',public_states[[name]],.GlobalEnv)
    method<-settings[[name]]$apeMethod
    offset<-if(is.null(fit$normalizationFactors))matrix(log(fit$sizeFactors),nrow(observations),ncol(observations),byrow=TRUE)else log(fit$normalizationFactors)
    warnings <- character()
    r<-withCallingHandlers(apeglm::apeglm(Y=observations,x=fit$design,
      param=if(method=='nbinomC*')fit$dispersions else parameters,offset=offset,coef=2L,
      mle=cbind(raw$log2FoldChange,raw$lfcSE)*log(2),threshold=.5*log(2),method=method),warning=function(w) {
        warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    list(log2FoldChange=r$map[,2L]/log(2),lfcSE=r$sd[,2L]/log(2),svalue=apeglm::svalue(r$thresh),
      warnings=warnings,random_state=.Random.seed)
  })
  names(original)<-names(settings)
} else if (beta_binomial) {
  analysis_dds <- rust[seq_len(256L),]
  fit <- rsdeseq2:::.rsdeseq2_fitted_view(analysis_dds)
  raw <- rsdeseq2::results(analysis_dds,name='conditionB')
  observations <- as.matrix(SummarizedExperiment::assay(analysis_dds,'counts'))
  parameters <- cbind(8,observations+20)
  settings <- list(betabinR=list(apeMethod='betabinR',param=parameters),
    betabinCR=list(apeMethod='betabinCR',param=parameters),betabinC=list(apeMethod='betabinC',param=parameters))
  original <- if(all(names(settings)%in%names(cached$reference)))cached$reference[names(settings)]else lapply(names(settings),function(method) {
    warnings <- character()
    r <- withCallingHandlers(apeglm::apeglm(Y=observations,x=fit$design,param=parameters,coef=2L,
      mle=cbind(raw$log2FoldChange,raw$lfcSE)*log(2),threshold=.5*log(2),method=method),warning=function(w) {
        warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    list(log2FoldChange=r$map[,2L]/log(2),lfcSE=r$sd[,2L]/log(2),svalue=apeglm::svalue(r$thresh),warnings=warnings)
  })
  names(original) <- names(settings)
} else {
  analysis_dds <- rust
  settings <- list(general=list(apeMethod='general'),nbinomR=list(apeMethod='nbinomR'))
  if (all(names(settings) %in% names(cached$reference))) {
    original <- cached$reference[names(settings)]
  } else {
    original_fit <- DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(fixture$counts,fixture$colData,~condition),fitType='mean',quiet=TRUE)
    coefficient <- tail(DESeq2::resultsNames(original_fit),1L)
    original <- lapply(settings,function(options) do.call(DESeq2::lfcShrink,
      c(list(dds=original_fit,coef=coefficient,type='apeglm',lfcThreshold=.5,quiet=TRUE),options)))
  }
  # This expression is imported as syntax. R never evaluates its body in a
  # candidate call, even though the function is accepted by the public API.
  custom <- function(y,x,beta,param,offset)stats::dnbinom(y,mu=exp(x%*%beta+offset),size=1/param,log=TRUE)
  settings$custom <- list(apeMethod='general',log.lik=custom)
}
public_warnings <- list()
actual <- lapply(names(settings),function(name) {
  warnings <- character()
  if(random_mode)assign('.Random.seed',public_states[[name]],.GlobalEnv)
  r <- withCallingHandlers(blocked(do.call(rsdeseq2::lfcShrink,
    c(list(dds=analysis_dds,coef='conditionB',type='apeglm',lfcThreshold=.5,quiet=TRUE),settings[[name]]))),warning=function(w) {
      warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  if(random_mode)stopifnot(identical(.Random.seed,original[[name]]$random_state))
  public_warnings[[name]]<<-warnings
  r
})
names(actual)<-names(settings)
errors <- lapply(names(settings),function(name) {
  a <- actual[[name]]; b <- original[[if(name=='custom')'general'else name]]
  stopifnot(identical(DESeq2::priorInfo(a)$engine,'Rust'))
  list(profile=name,mode=compare(a$log2FoldChange,b$log2FoldChange,.01),
    sd=compare(a$lfcSE,b$lfcSE,.01),threshold_svalue=compare(a$svalue,b$svalue,.05))
})
rng_error_checks <- NULL
if(random_mode) {
  name <- names(settings)[[1L]]
  before <- public_states[[name]]
  expected_after <- original[[name]]$random_state
  assign('.Random.seed',before,.GlobalEnv)
  overflow <- matrix(.Machine$double.xmax,nrow(analysis_dds),ncol(analysis_dds))
  late <- tryCatch(suppressWarnings(blocked(rsdeseq2::lfcShrink(analysis_dds,coef='conditionB',type='apeglm',
    apeMethod='nbinomC*',Y=overflow,quiet=TRUE))),error=identity)
  stopifnot(inherits(late,'error'),grepl('representable|dispersions',conditionMessage(late)),identical(.Random.seed,expected_after))
  unsupported <- before;unsupported[[1L]] <- 10203L
  assign('.Random.seed',unsupported,.GlobalEnv)
  early <- tryCatch(blocked(rsdeseq2::lfcShrink(analysis_dds,coef='conditionB',type='apeglm',apeMethod='nbinomC*',quiet=TRUE)),error=identity)
  stopifnot(inherits(early,'error'),grepl('Inversion-normal',conditionMessage(early)),identical(.Random.seed,unsupported))
  rm('.Random.seed',envir=.GlobalEnv)
  unseeded <- suppressWarnings(blocked(rsdeseq2::lfcShrink(analysis_dds,coef='conditionB',type='apeglm',apeMethod='nbinomC*',quiet=TRUE)))
  stopifnot(is.integer(.Random.seed),length(.Random.seed)==626L,.Random.seed[[1L]]==10403L,
    identical(DESeq2::priorInfo(unseeded)$engine,'Rust'))
  rng_error_checks <- list(late_error_retains_consumed_state=TRUE,unsupported_state_rejects_without_changes=TRUE,
    unseeded_request_initializes_in_Rust=TRUE,R_random_calls_blocked=TRUE)
}
saveRDS(list(rust=rust,reference=original,actual=actual),file.path(out,'paired_fits.rds'))
jsonlite::write_json(list(R=as.character(getRversion()),DESeq2=as.character(packageVersion('DESeq2')),
  apeglm=as.character(packageVersion('apeglm')),genes=nrow(analysis_dds),complete_native_standalone_batches=length(wrapped),
  warning_messages_identical=sum(vapply(wrapped,function(x)isTRUE(x$original_warning_messages_identical),logical(1))),
  warning_events_identical=sum(vapply(wrapped,function(x)isTRUE(x$original_warning_events_identical),logical(1))),
  public_warning_counts=lapply(public_warnings,length),
  original_public_warning_counts=if(beta_binomial||random_mode)lapply(original,function(v)length(v$warnings))else NULL,
  public_random_states_identical=if(random_mode)TRUE else NULL,rng_error_checks=rng_error_checks,
  R_calculations_blocked=TRUE,existing_bounds_unchanged=TRUE,public_errors=errors),
  file.path(out,'validation.json'),auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
cat(length(wrapped),'regression batches match standalone bits or reject together; public method checks pass with R math blocked\n')
