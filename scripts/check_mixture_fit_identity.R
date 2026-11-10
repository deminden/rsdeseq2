#!/usr/bin/env Rscript
# Conformance check for complete native prior fits, separate from public lfcShrink.
args<-commandArgs(TRUE)
if(length(args)!=3L)stop('usage: STANDALONE_DRIVER R_LIBRARY FRESH_OUTPUT_DIR')
driver<-normalizePath(args[[1L]]);library(rsdeseq2,lib.loc=normalizePath(args[[2L]]))
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
out<-args[[3L]];if(dir.exists(out))stop('Use a fresh identity output directory');dir.create(out,recursive=TRUE)
reference<-jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/mixture_fit_r462.json'),simplifyVector=FALSE)
options<-list(max_iterations=1000L,gradient_tolerance='1e-8',em_iterations=20L,active_iterations=100L,zero_threshold='1e-8',likelihood_offset='1e-6')
cases<-lapply(reference$cases,function(c){c$options<-options;c['initial_weights']<-list(NULL);c})
# Transport the actual R input bits. R's decimal reader and Rust's correctly
# rounded parser can differ on 17-digit fixture strings; identical expressions
# in text are insufficient to establish identical numeric engine input.
inputCases<-lapply(cases,function(c){
 for(name in c('log_likelihoods','observation_weights','dirichlet'))
   c[[name]]<-lapply(as.numeric(unlist(c[[name]])),function(v)paste0('bits:',hex(v)))
 for(name in c('gradient_tolerance','zero_threshold','likelihood_offset'))
   c$options[[name]]<-paste0('bits:',hex(as.numeric(c$options[[name]])))
 c
})
input<-list(mixture_fit_cases=inputCases);path<-file.path(out,'input.json')
jsonlite::write_json(input,path,auto_unbox=TRUE,null='null',digits=NA)
output<-system2(driver,shQuote(normalizePath(path)),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(output,'status')))stop(paste(output,collapse='\n'))
standalone<-jsonlite::fromJSON(paste(output,collapse='\n'),simplifyVector=FALSE)
blocked<-function(expression){
 functions<-list(base=c('sqrt','exp','log','mean','sum'),stats=c('optim','dnorm'))
 on.exit(for(package in names(functions))for(name in functions[[package]])suppressMessages(untrace(name,where=asNamespace(package))),add=TRUE)
 for(package in names(functions))for(name in functions[[package]])suppressWarnings(suppressMessages(trace(name,where=asNamespace(package),print=FALSE,
   tracer=quote(stop('prior bridge executed an R calculation',call.=FALSE)))))
 force(expression)
}
guard<-tryCatch(blocked(exp(1)),error=identity);stopifnot(inherits(guard,'error'),grepl('executed an R calculation',conditionMessage(guard)))
wrapped<-lapply(cases,function(c){
 request<-c(list(operation='posterior_fit_mixture'),c[c('observations','components')],
   lapply(c[c('log_likelihoods','observation_weights','dirichlet')],function(v)as.numeric(unlist(v))))
 request['initial_weights']<-list(NULL)
 request$options<-options
 request$options[c('gradient_tolerance','zero_threshold','likelihood_offset')]<-lapply(options[c('gradient_tolerance','zero_threshold','likelihood_offset')],as.numeric)
 result<-blocked(rsdeseq2:::.rsdeseq2_formula_request(request))
 result$weights<-lapply(unname(result$weights),hex)
 for(name in c('log_likelihood','penalized_log_likelihood','dual_residual','optimizer_weight_sum','function_evaluations'))result[[name]]<-hex(result[[name]])
 c(list(name=c$name),result)
})
canonical<-function(x){if(!is.list(x))return(x);if(!is.null(names(x)))x<-x[sort(names(x))];lapply(x,canonical)}
wrapped<-.formula_canonical(list(cases=wrapped))
jsonlite::write_json(wrapped,file.path(out,'wrapped.json'),auto_unbox=TRUE,null='null',pretty=TRUE)
writeLines(output,file.path(out,'standalone.json'))
if(!identical(canonical(wrapped),canonical(standalone)))stop(paste(all.equal(canonical(wrapped),canonical(standalone),tolerance=0),collapse='; '))
cat('8 complete mixture prior fits: all weights, likelihoods and diagnostics bit-identical; R calculations blocked\n')
