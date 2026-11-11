#!/usr/bin/env Rscript
# Complete native/wrapped comparison; original numerical oracles are tested
# separately by location_shrinkage.rs, including exact NA/NaN distinctions.
args <- commandArgs(TRUE)
if(!length(args)%in%c(3L,4L))stop('usage: STANDALONE_DRIVER R_LIBRARY FRESH_OUTPUT_DIR [REFERENCE_JSON]')
driver <- normalizePath(args[[1L]])
library(rsdeseq2,lib.loc=normalizePath(args[[2L]]))
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
out <- args[[3L]];if(dir.exists(out))stop('Use a fresh output directory');dir.create(out,recursive=TRUE)
references <- if(length(args)==4L)list(jsonlite::fromJSON(args[[4L]],simplifyVector=FALSE))else lapply(c('location_shrinkage_r462.json','custom_likelihood_r462.json'),function(file)jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data',file),simplifyVector=FALSE))
cases <- unlist(lapply(references,`[[`,'cases'),recursive=FALSE)
fields <- c('estimates','standard_errors','thresholds','family','point_mass','nullweight','scales','mode','range','alpha','prune_threshold','observation_weights','fixed_prior','supplied_prior','df')
inputs <- lapply(cases,function(c){wanted<-c('name',fields,intersect(c('optimizer','control','penalty'),names(c)),if(!is.null(c$custom_likelihood))'custom_likelihood',if(!is.null(c$builtin_likelihood))'builtin_likelihood');stats::setNames(lapply(wanted,function(n)c[[n]]),wanted)})
path <- file.path(out,'input.json');jsonlite::write_json(list(location_cases=inputs),path,auto_unbox=TRUE,null='null',digits=NA)
output <- system2(driver,shQuote(normalizePath(path)),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(output,'status')))stop(paste(output,collapse='\n'))
standalone <- jsonlite::fromJSON(paste(output,collapse='\n'),simplifyVector=FALSE)
decode <- function(v)vapply(v,function(h){h<-sub('^bits:','',h);readBin(as.raw(strtoi(substring(h,seq(1L,15L,2L),seq(2L,16L,2L)),16L)),'double',n=1L,size=8L,endian='big')},numeric(1L))
function_input <- function(f){
 if(is.null(f))return(NULL)
 list(arguments=unlist(f$arguments,use.names=FALSE),expression=.formula_ast_decode(f$expression),constants=lapply(f$constants,.formula_decode))
}
blocked <- function(expression){
 functions <- list(base=c('sqrt','exp','expm1','log','log1p','log2','mean','sum','sort','order','rank','cumsum','^','lgamma','digamma'),stats=c('pnorm','dnorm','pt','dt','integrate','pf','df','pgamma','dgamma','pbeta','dbeta'))
 on.exit(for(package in names(functions))for(name in functions[[package]])suppressMessages(untrace(name,where=asNamespace(package))),add=TRUE)
 for(package in names(functions))for(name in functions[[package]])suppressWarnings(suppressMessages(trace(name,where=asNamespace(package),print=FALSE,tracer=quote(stop('location bridge executed an R calculation',call.=FALSE)))))
 force(expression)
}
guard <- tryCatch(blocked(sqrt(2)),error=identity);stopifnot(inherits(guard,'error'),grepl('executed an R calculation',conditionMessage(guard)))
encode <- function(v){if(is.numeric(v))return(lapply(unname(v),hex));if(is.logical(v)||is.character(v))return(as.list(unname(v)));if(is.list(v))return(lapply(v,encode));v}
wrapped <- lapply(seq_along(inputs),function(i){
 c <- inputs[[i]];request <- c(list(operation='shrink_location_batch'),c[c(fields,intersect(c('optimizer','control','penalty'),names(c)))])
 for(field in c('estimates','standard_errors','thresholds','nullweight','scales','mode','range','alpha','prune_threshold','observation_weights','df'))if(!is.null(request[[field]]))request[[field]]<-decode(request[[field]])
 if(!is.null(request$penalty)&&!identical(request$penalty,'unit'))request$penalty<-decode(request$penalty)
 if(!is.null(request$control))request$control<-lapply(request$control,function(v)if(is.character(v)&&length(v)==1L&&v%in%c('Inf','-Inf','NaN'))as.numeric(v)else v)
 if(!is.null(request$supplied_prior))for(field in c('weights','means','sds','lower','upper'))if(!is.null(request$supplied_prior[[field]]))request$supplied_prior[[field]]<-decode(request$supplied_prior[[field]])
 if(!is.null(c$custom_likelihood))request$custom_likelihood <- c(list(name=c$custom_likelihood$name),lapply(c$custom_likelihood[c('log_cdf','log_density','first','second')],function_input))
 if(!is.null(c$builtin_likelihood)){
   model<-c$builtin_likelihood
   request$builtin_likelihood<-if(model$type=='normalmix')list(type='normalmix',weights=decode(model$weights),scales=.formula_decode(model$scales))else if(model$type=='logF')list(type='logF',df1=decode(model$df1),df2=decode(model$df2))else {v<-model;for(n in intersect(c('y','scale','trials','capture_y','capture_trials'),names(v)))v[[n]]<-decode(v[[n]]);v}
 }
 result <- tryCatch(blocked(rsdeseq2:::.rsdeseq2_formula_request(request)),error=identity)
 if(!identical(inherits(result,'error'),!is.null(cases[[i]]$expected_error)))stop(c$name,': ',if(inherits(result,'error'))conditionMessage(result)else 'unexpected success')
 if(inherits(result,'error'))return(list(name=c$name,ok=FALSE))
 list(name=c$name,ok=TRUE,value=encode(result))
})
canonical <- function(x){if(!is.list(x))return(x);if(!is.null(names(x)))x<-x[sort(names(x))];lapply(x,canonical)}
wrapped <- .formula_canonical(list(cases=wrapped))
jsonlite::write_json(wrapped,file.path(out,'wrapped.json'),auto_unbox=TRUE,pretty=TRUE,null='null');writeLines(output,file.path(out,'standalone.json'))
stopifnot(identical(canonical(wrapped),canonical(standalone)))
jsonlite::write_json(list(complete_native_standalone_batches=length(inputs),R_calculations_blocked=TRUE,reference=if(length(args)==4L)basename(args[[4L]])else 'accepted location/custom'),file.path(out,'validation.json'),auto_unbox=TRUE,pretty=TRUE)
cat(length(inputs),'complete location/custom likelihood batches: every output bit, warning and acceptance decision matches; R calculations blocked\n')
