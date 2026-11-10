#!/usr/bin/env Rscript
# Complete direct optimizer diagnostics plus actual public R workflow checks.
args<-commandArgs(TRUE)
if(!length(args)%in%c(3L,4L))stop('usage: DRIVER R_LIBRARY FRESH_OUTPUT_DIR [PAIRED_FITS_RDS]')
driver<-normalizePath(args[[1L]])
library(rsdeseq2,lib.loc=normalizePath(args[[2L]]))
requireNamespace('ashr');requireNamespace('apeglm')
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
source(file.path(root,'scripts/native_posterior_check_helpers.R'))
out<-args[[3L]];if(dir.exists(out))stop('Use a fresh output directory');dir.create(out,recursive=TRUE)
reference<-jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/ashr_optimizer_r462.json'),simplifyVector=FALSE)
input<-file.path(out,'input.json');jsonlite::write_json(list(em_fit_cases=reference$em_fit_cases),input,auto_unbox=TRUE,null='null',digits=NA)
output<-system2(driver,shQuote(normalizePath(input)),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(output,'status')))stop(paste(output,collapse='\n'))
writeLines(output,file.path(out,'standalone.json'))
standalone<-jsonlite::fromJSON(paste(output,collapse='\n'),simplifyVector=FALSE)
controls<-function(c)lapply(c,function(v)if(is.character(v)&&length(v)==1L&&v%in%c('Inf','-Inf','NaN'))as.numeric(v)else v)
canonical<-function(v){if(is.list(v)){if(!is.null(names(v)))v<-v[order(names(v))];return(lapply(v,canonical))};.formula_canonical(v)}
wrapped<-lapply(seq_along(reference$em_fit_cases),function(i){
 c<-reference$em_fit_cases[[i]]
 request<-c(list(operation='posterior_fit_em'),c[c('observations','components','log_likelihoods','observation_weights','dirichlet','initial_weights','optimizer','control')])
 for(field in c('log_likelihoods','observation_weights','dirichlet','initial_weights'))if(!is.null(request[[field]]))request[[field]]<-decode(request[[field]])
 request$control<-controls(request$control)
 r<-blocked(rsdeseq2:::.rsdeseq2_formula_request(request))
 value<-list(weights=encode(r$weights),log_likelihood=hex(r$log_likelihood),penalized_log_likelihood=hex(r$penalized_log_likelihood),converged=r$converged,
   iterations=if(is.null(r$iterations))NULL else as.integer(r$iterations),evaluations=as.integer(r$evaluations),residual=hex(r$residual),weight_sum=hex(r$weight_sum),
   posterior_concentrations=if(is.null(r$posterior_concentrations))NULL else encode(r$posterior_concentrations),optimizer_objective=if(is.null(r$optimizer_objective))NULL else hex(r$optimizer_objective),message=r$message)
 s<-standalone$cases[[i]];stopifnot(is.null(s$error))
 if(!identical(canonical(value),canonical(s$value)))stop(c$name,': ',paste(all.equal(canonical(value),canonical(s$value)),collapse='; '))
 list(name=c$name,value=value)
})
jsonlite::write_json(wrapped,file.path(out,'wrapped.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
public_errors<-NULL
if(length(args)==4L){
 cached<-readRDS(args[[4L]]);dds<-cached$rust[seq_len(256L),]
 raw<-rsdeseq2::results(dds,name='conditionB')
 settings<-list(mixEM=list(optmethod='mixEM'),w_mixEM=list(optmethod='w_mixEM',weights=rep(c(.5,1,2,1),64)),
   mixVBEM=list(optmethod='mixVBEM'),cxxMixSquarem=list(optmethod='cxxMixSquarem'))
 original<-lapply(settings,function(options){
   r<-suppressWarnings(do.call(ashr::ash,c(list(betahat=raw$log2FoldChange,sebetahat=raw$lfcSE,mixcompdist='normal',method='shrink'),options)))
   list(mean=ashr::get_pm(r),sd=ashr::get_psd(r),svalue=ashr::get_svalue(r))
 })
 actual<-lapply(settings,function(options)suppressWarnings(blocked(do.call(rsdeseq2::lfcShrink,c(list(dds=dds,coef='conditionB',type='ashr',svalue=TRUE,quiet=TRUE),options)))))
 public_errors<-lapply(names(settings),function(name){a<-actual[[name]];b<-original[[name]];stopifnot(identical(DESeq2::priorInfo(a)$engine,'Rust'))
   list(profile=name,mean=compare(a$log2FoldChange,b$mean,1e-3),sd=compare(a$lfcSE,b$sd,1e-3),svalue=compare(a$svalue,b$svalue,.02))})
 saveRDS(list(rust=cached$rust,reference=original,actual=actual),file.path(out,'paired_fits.rds'))
}
jsonlite::write_json(list(complete_optimizer_batches=length(wrapped),R_calculations_blocked=TRUE,public_errors=public_errors,existing_bounds_unchanged=TRUE),file.path(out,'validation.json'),auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
cat('Native EM optimizer diagnostics match standalone bits:',length(wrapped),'batches\n')
