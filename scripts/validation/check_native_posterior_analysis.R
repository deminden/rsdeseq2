#!/usr/bin/env Rscript
# Complete direct/wrapped checks; all candidate statistical R calls are blocked.
args<-commandArgs(TRUE)
if(!length(args)%in%c(3L,4L))stop('usage: DRIVER R_LIBRARY FRESH_OUTPUT_DIR [PAIRED_FITS_RDS]')
driver<-normalizePath(args[[1L]])
library(rsdeseq2,lib.loc=normalizePath(args[[2L]]))
requireNamespace('apeglm')
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..', '..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
source(file.path(root,'scripts/validation/native_posterior_check_helpers.R'))
out<-args[[3L]]
if(dir.exists(out))stop('Use a fresh output directory')
dir.create(out,recursive=TRUE)
reference<-jsonlite::fromJSON(file.path(root,'crates/rsdeseq2/tests/data/posterior_analysis_r462.json'),simplifyVector=FALSE)
input<-file.path(out,'input.json')
jsonlite::write_json(list(regression_cases=reference$cases),input,auto_unbox=TRUE,null='null')
output<-system2(driver,shQuote(normalizePath(input)),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(output,'status')))stop(paste(output,collapse='\n'))
writeLines(output,file.path(out,'standalone.json'))
standalone<-jsonlite::fromJSON(paste(output,collapse='\n'),simplifyVector=FALSE)
.posterior_canonical<-function(x) {
 if(is.list(x)) {if(!is.null(names(x)))x<-x[order(names(x))];return(lapply(x,.posterior_canonical))}
 .formula_canonical(x)
}
wrapped<-lapply(seq_along(reference$cases),function(i){
 c<-reference$cases[[i]];a<-c$analysis;model<-c$likelihood
 if(!is.null(model)){
   model$arguments<-unlist(model$arguments);model$expression<-.formula_ast_decode(model$expression)
   model$constants<-lapply(model$constants,.formula_decode)
 }
 arguments<-list(Y=matrix_input(c$observations),x=matrix_input(c$design),param=matrix_input(c$parameters),
   coef=c$coefficient+1L,weights=matrix_input(c$weights),offset=matrix_input(c$offsets),
   prior.control=list(no.shrink=as.numeric(unlist(c$nuisance))+1L,prior.mean=decode(c$prior_mean),prior.scale=decode(c$scale),prior.df=decode(c$prior_df),prior.no.shrink.mean=decode(c$nuisance_mean),prior.no.shrink.scale=decode(c$nuisance_scale)),
   method=c$method,optim.method=c$optimizer,bounds=decode(c$bounds),log.link=c$log_link,
   interval.type=a$kind,interval.level=decode(a$level),ngrid=a$grid_intervals,nsd=decode(a$grid_standard_deviations),ngrid.nuis=a$nuisance_points,nsd.nuis=decode(a$nuisance_standard_deviations),
   contrasts=matrix_input(a$contrasts),param.sd=if(is.null(a$parameter_standard_errors))NULL else decode(a$parameter_standard_errors),
   log.lik=if(is.null(model))NULL else list(type='restricted',model=model))
 if(c$method=='general'&&is.null(model))arguments$log.lik<-list(type='negative_binomial')
 r<-blocked(rsdeseq2:::.rsdeseq2_formula_request(list(operation='apeglm_shrinkage',arguments=arguments)))
 s<-standalone$cases[[i]];stopifnot(is.null(s$error));s<-s$value
 grid<-lapply(r$analysis$grid,function(d)if(is.null(d))NULL else list(monotone=d$monotone,left=if(is.na(d$left))NULL else hex(d$left),right=if(is.na(d$right))NULL else hex(d$right),evaluations=as.integer(d$evaluations)))
 rowbits<-function(m)if(is.null(m))NULL else encode(t(m))
 actual<-list(intervals=rowbits(r$analysis$intervals),contrast_mode=rowbits(r$analysis$contrast_mode),contrast_sd=rowbits(r$analysis$contrast_sd),grid=grid)
 if(!identical(.posterior_canonical(actual),.posterior_canonical(s$analysis))) {
   saveRDS(list(name=c$name,actual=actual,expected=s$analysis),file.path(out,"mismatch.rds"))
   stop(c$name, ": ", paste(all.equal(.posterior_canonical(actual),.posterior_canonical(s$analysis)),collapse="; "))
 }
 stopifnot(
   identical(encode(r$posterior_mode),s$summary_mode),identical(encode_optional(r$posterior_sd),s$summary_sd),
   identical(encode_optional(r$svalue),s$sign_svalues))
 list(name=c$name,analysis=actual,mode=encode(r$posterior_mode),sd=encode_optional(r$posterior_sd),svalue=encode_optional(r$svalue))
})
jsonlite::write_json(wrapped,file.path(out,'wrapped.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
public_errors<-NULL
if(length(args)==4L) {
 cached<-readRDS(args[[4L]])
 dds<-cached$rust[seq_len(128L),]
 raw<-rsdeseq2::results(dds,name='conditionB')
 fit<-rsdeseq2:::.rsdeseq2_fitted_view(dds)
 observations<-as.matrix(SummarizedExperiment::assay(dds,'counts'))
 offset<-if(is.null(fit$normalizationFactors))matrix(log(fit$sizeFactors),nrow(observations),ncol(observations),byrow=TRUE)else log(fit$normalizationFactors)
 custom<-function(y,x,beta,param,offset).9*stats::dnbinom(y,mu=exp(x%*%beta+offset),size=1/param,log=TRUE)
 settings<-list(HPD=list(interval.type='HPD'),credible=list(interval.type='credible',contrasts=cbind(c(0,1),c(-1,1))),
   uncertain=list(interval.type='credible',param.sd=.05*fit$dispersions),
   custom_grid=list(interval.type='credible',log.lik=custom))
 original<-lapply(settings,function(options) {
   arguments<-utils::modifyList(list(Y=observations,x=fit$design,param=fit$dispersions,offset=offset,coef=2L,
     mle=cbind(raw$log2FoldChange,raw$lfcSE)*log(2),method='nbinomCR',log.lik=apeglm::logLikNB),options)
   r<-suppressWarnings(do.call(apeglm::apeglm,arguments))
   list(mode=r$map[,2L]/log(2),sd=r$sd[,2L]/log(2),svalue=r$svalue)
 })
 actual<-lapply(settings,function(options) suppressWarnings(blocked(do.call(rsdeseq2::lfcShrink,
   c(list(dds=dds,coef='conditionB',type='apeglm',apeMethod='nbinomCR',lfcThreshold=0,svalue=TRUE,quiet=TRUE),options)))))
 public_errors<-lapply(names(settings),function(name) {
   a<-actual[[name]];b<-original[[name]]
   stopifnot(identical(DESeq2::priorInfo(a)$engine,'Rust'))
   list(profile=name,mode=compare(a$log2FoldChange,b$mode,.01),sd=compare(a$lfcSE,b$sd,.01),svalue=compare(a$svalue,as.numeric(b$svalue),.02))
 })
 saveRDS(list(rust=cached$rust,reference=original,actual=actual),file.path(out,'paired_fits.rds'))
}
jsonlite::write_json(list(complete_native_standalone_batches=length(wrapped),R_calculations_blocked=TRUE,original_reference_profiles=length(reference$cases),public_errors=public_errors,fixed_grid_leaves=length(reference$grid_leaves),existing_posterior_bounds_unchanged=TRUE),file.path(out,'validation.json'),auto_unbox=TRUE,pretty=TRUE)
cat('Posterior analysis complete-output identity passed:',length(wrapped),'batches\n')
