#!/usr/bin/env Rscript
# Generic actual-bridge/standalone kernel identity with R math blocked.
args<-commandArgs(TRUE)
if(length(args)!=4L)stop('usage: REFERENCE_JSON DRIVER R_LIBRARY FRESH_OUTPUT_DIR')
reference<-jsonlite::fromJSON(args[[1L]],simplifyVector=FALSE)
driver<-normalizePath(args[[2L]]);library(rsdeseq2,lib.loc=normalizePath(args[[3L]]))
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
source(file.path(root,'scripts/native_posterior_check_helpers.R'))
out<-args[[4L]];if(dir.exists(out))stop('Use a fresh output directory');dir.create(out,recursive=TRUE)
path<-file.path(out,'input.json');jsonlite::write_json(list(builtin_likelihood_cases=reference$kernel_cases),path,auto_unbox=TRUE,null='null',digits=NA)
output<-system2(driver,shQuote(normalizePath(path)),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(output,'status')))stop(paste(output,collapse='\n'))
writeLines(output,file.path(out,'standalone.json'));standalone<-jsonlite::fromJSON(paste(output,collapse='\n'),simplifyVector=FALSE)
wrapped<-lapply(reference$kernel_cases,function(c){
  m<-c$builtin_likelihood
  for(n in intersect(c('weights','df1','df2','y','scale','trials','capture_y','capture_trials'),names(m)))m[[n]]<-if(n=='weights'&&!is.null(m[[n]]$type)).formula_decode(m[[n]])else decode(m[[n]])
  if(!is.null(m$scales))m$scales<-.formula_decode(m$scales)
  r<-tryCatch(blocked(rsdeseq2:::.rsdeseq2_formula_request(list(operation='posterior_builtin_likelihood',likelihood=m,
    arguments=decode(c$arguments),lower=decode(c$lower),upper=decode(c$upper)))),error=identity)
  stopifnot(identical(inherits(r,'error'),!is.null(c$expected_error)))
  if(inherits(r,'error'))return(list(name=c$name,ok=FALSE))
  list(name=c$name,ok=TRUE,value=lapply(r,encode))
})
canonical<-function(v){if(!is.list(v))return(v);if(!is.null(names(v)))v<-v[order(names(v))];lapply(v,canonical)}
wrapped<-.formula_canonical(list(cases=wrapped));jsonlite::write_json(wrapped,file.path(out,'wrapped.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
stopifnot(identical(canonical(wrapped),canonical(standalone)))
jsonlite::write_json(list(kernel_batches=length(wrapped$cases),R_calculations_blocked=TRUE,reference=basename(args[[1L]]),native_wrapped_bits='exact'),file.path(out,'validation.json'),auto_unbox=TRUE,pretty=TRUE)
cat(length(wrapped$cases),'kernel batches have exact wrapped/native bits; R calculations blocked\n')
