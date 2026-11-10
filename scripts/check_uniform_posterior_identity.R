#!/usr/bin/env Rscript
# Independent accuracy plus actual native-bridge/R-free bit conformance.
args<-commandArgs(TRUE)
if(!length(args)%in%c(3L,4L))stop('usage: STANDALONE_DRIVER R_LIBRARY FRESH_OUTPUT_DIR [ORIGINAL_REFERENCE_JSON]')
driver<-normalizePath(args[[1L]]);library(rsdeseq2,lib.loc=normalizePath(args[[2L]]))
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
out<-args[[3L]];if(dir.exists(out))stop('Use a fresh identity output directory');dir.create(out,recursive=TRUE)
referencePath<-if(length(args)==4L)args[[4L]]else file.path(root,'crates/rsdeseq2/tests/data/uniform_posterior_r462.json')
reference<-jsonlite::fromJSON(referencePath,simplifyVector=FALSE)
stopifnot(identical(reference$source$input_encoding,'f64_bits'))
truncated<-!is.null(reference$cases[[1L]]$means)
fields<-c('weights',if(!is.null(reference$cases[[1L]]$df))'df',if(truncated)c('means','sds'),'lower','upper','estimates','standard_errors','thresholds')
cases<-lapply(reference$cases,function(c)c[c('name',fields)])
path<-file.path(out,'input.json');jsonlite::write_json(setNames(list(cases),if(truncated)'truncated_cases'else'uniform_cases'),path,auto_unbox=TRUE,null='null',digits=NA)
output<-system2(driver,shQuote(normalizePath(path)),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(output,'status')))stop(paste(output,collapse='\n'))
standalone<-jsonlite::fromJSON(paste(output,collapse='\n'),simplifyVector=FALSE)
decode<-function(v)vapply(v,function(h){
 h<-sub('^bits:','',h)
 readBin(as.raw(strtoi(substring(h,seq(1L,15L,2L),seq(2L,16L,2L)),16L)),'double',n=1L,size=8L,endian='big')
},numeric(1L))
blocked<-function(expression){
 functions<-list(base=c('sqrt','exp','log','mean','sum'),stats=c('pnorm','dnorm'))
 on.exit(for(package in names(functions))for(name in functions[[package]])suppressMessages(untrace(name,where=asNamespace(package))),add=TRUE)
 for(package in names(functions))for(name in functions[[package]])suppressWarnings(suppressMessages(trace(name,where=asNamespace(package),print=FALSE,
  tracer=quote(stop('uniform bridge executed an R calculation',call.=FALSE)))))
 force(expression)
}
guard<-tryCatch(blocked(exp(1)),error=identity);stopifnot(inherits(guard,'error'),grepl('executed an R calculation',conditionMessage(guard)))
close<-function(a,b,absolute,relative,context){
 b<-as.numeric(unlist(b))
 if(!identical(is.nan(a),is.nan(b))||!identical(is.infinite(a),is.infinite(b)))stop(context,': missingness differs')
 infinite<-is.infinite(b);if(any(a[infinite]!=b[infinite]))stop(context,': infinity signs differ')
 finite<-is.finite(b)
 if(any(abs(a[finite]-b[finite])>absolute+relative*abs(b[finite])))stop(context,': independent numerical bound exceeded')
}
wrapped<-lapply(seq_along(cases),function(i){
 c<-cases[[i]];request<-c(list(operation=if(truncated)'posterior_truncated_normal_mixture'else'posterior_uniform_mixture'),lapply(c[fields],decode))
 result<-blocked(rsdeseq2:::.rsdeseq2_formula_request(request))
 expected<-reference$cases[[i]]
 close(result$log_likelihoods,expected$log_likelihoods,1e-8,1e-12,paste(c$name,'likelihood'))
 close(result$posterior_mean,expected$posterior_mean,2e-6,0,paste(c$name,'mean'))
 close(result$posterior_sd,expected$posterior_sd,2e-6,0,paste(c$name,'SD'))
 close(result$posterior_weights,expected$posterior_weights,2e-8,0,paste(c$name,'weights'))
 close(result$cdf,expected$cdf,2e-8,0,paste(c$name,'CDF'))
 c(list(name=c$name),lapply(result,function(v)lapply(unname(v),hex)))
})
canonical<-function(x){if(!is.list(x))return(x);if(!is.null(names(x)))x<-x[sort(names(x))];lapply(x,canonical)}
wrapped<-.formula_canonical(list(cases=wrapped))
jsonlite::write_json(wrapped,file.path(out,'wrapped.json'),auto_unbox=TRUE,null='null',pretty=TRUE)
writeLines(output,file.path(out,'standalone.json'))
stopifnot(identical(canonical(wrapped),canonical(standalone)))
cat(length(cases), if(truncated)'truncated-normal'else'uniform', 'mixture profiles: independent bounds pass and every output bit matches; R calculations blocked\n')
