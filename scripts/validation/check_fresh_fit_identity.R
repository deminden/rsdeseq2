#!/usr/bin/env Rscript
# Complete fresh numerical fit state: actual public wrapper versus an R-free process.
args<-commandArgs(TRUE)
if(length(args)!=3L)stop('usage: STANDALONE_DRIVER R_LIBRARY FRESH_OUTPUT_DIR')
driver<-normalizePath(args[[1L]]);library(rsdeseq2,lib.loc=normalizePath(args[[2L]]))
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..', '..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-fit-data.R'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
source(file.path(root,'scripts/validation/native_posterior_check_helpers.R'))
source(file.path(root,'scripts/validation/fresh_fit_check_helpers.R'))
out<-args[[3L]];if(dir.exists(out))stop('Use a fresh output directory');dir.create(out,recursive=TRUE)
counts<-.rsdeseq2_test_counts();counts[nrow(counts),]<-0
design<-cbind(Intercept=1,conditionB=rep(0:1,each=4));reduced<-matrix(1,8,1)
size<-exp(seq(log(.8),log(1.2),length.out=8))
factors<-outer(seq(.9,1.1,length.out=nrow(counts)),size)
weights<-matrix(rep(c(.3,1,.7,1,1,.4,1,.8),each=nrow(counts)),nrow(counts),8)
failed<-weights;failed[39,]<-c(1,rep(0,7))
replaced<-counts;replaced[1,8]<-2000;replaced[2,]<-0;replaced[2,8]<-1000
newzero<-counts;newzero[2,]<-0;newzero[2,8]<-1000
large<-.rsdeseq2_large_transform_fixture()
profiles<-list(
 mean_estimated=list(),lrt_supplied=list(test='LRT',size=size),
 lrt_two_df=list(test='LRT',design=cbind(design,dose=rep(1:4,2)),factors=factors,weights=weights),
 parametric=list(fit_type='parametric',counts=large$counts,design=cbind(Intercept=1,conditionB=rep(0:1,each=6))),
 local=list(fit_type='local'),
 weighted_failed=list(weights=failed),weighted_local_lrt=list(test='LRT',fit_type='local',weights=weights),
 weighted_failed_lrt=list(test='LRT',weights=failed,factors=factors),
 weighted_failed_no_refit_wald=list(weights=failed,size=size,min_replicates=3L),
 weighted_failed_no_refit_lrt=list(test='LRT',weights=failed,factors=factors,min_replicates=3L),
 factor_precedence=list(size=size,factors=factors),multiworker=list(threads=4L),
 enabled_no_refit_wald=list(min_replicates=3L),enabled_no_refit_lrt=list(test='LRT',min_replicates=3L),
 weighted_enabled_no_refit=list(min_replicates=3L,weights=weights),
 replacement_wald=list(counts=replaced,min_replicates=3L),replacement_lrt=list(test='LRT',counts=replaced,min_replicates=3L),
 replacement_weighted_wald=list(counts=replaced,min_replicates=3L,weights=weights),
 replacement_weighted_lrt=list(test='LRT',counts=replaced,min_replicates=3L,weights=weights,factors=factors),
 replacement_weighted_failed_wald=list(counts=replaced,min_replicates=3L,weights=failed,size=size),
 replacement_weighted_failed_lrt=list(test='LRT',counts=replaced,min_replicates=3L,weights=failed,factors=factors),
 newzero_wald=list(counts=newzero,min_replicates=3L),newzero_lrt=list(test='LRT',counts=newzero,min_replicates=3L),
 newzero_weighted_failed_wald=list(counts=newzero,min_replicates=3L,weights=failed,size=size),
 newzero_weighted_failed_lrt=list(test='LRT',counts=newzero,min_replicates=3L,weights=failed,factors=factors),
 unreplaceable=list(counts=replaced,min_replicates=7L),
 weighted_unreplaceable=list(counts=replaced,min_replicates=7L,weights=weights),
 bad_reduced=list(test='LRT',reduced=matrix(1:8,8,1),error='nested'),
 bad_rank=list(design=matrix(1,8,2),error='rank'),
 bad_counts=list(counts=replace(counts,1,.5),error='count'))
profiles$parametric_boundary<-list(fit_type='parametric',allow_error='dispersion')
query_cases<-function(p){
 n<-if(is.null(p$design))ncol(design)else ncol(p$design)
 contrast<-c(rep(0,n-1L),1)
 query<-function(name,contrast,cooks='auto',filtering=TRUE,alpha=.1,threshold=0,alternative='greaterAbs',df=NULL)
  list(name=name,contrast=bits(contrast),cooks=cooks,filtering=filtering,alpha=bits(alpha)[[1L]],
   threshold=bits(threshold)[[1L]],alternative=alternative,df=if(is.null(df))NULL else bits(df)[[1L]],missing=bits(NA_real_)[[1L]])
 q<-list(query('forward',contrast),query('reverse_unmasked',-contrast,cooks='off',filtering=FALSE),
  query('fractional_zero_cutoff',contrast*.5,cooks='0',alpha=.05))
 if(!identical(p$test,'LRT'))q<-c(q,list(query('threshold',contrast,threshold=1),
  query('threshold_t',contrast,threshold=.5,alternative='less',df=5)))
 q
}
cases<-lapply(names(profiles),function(name){
 p<-profiles[[name]]
 list(name=name,counts=matrix_descriptor(if(is.null(p$counts))counts else p$counts),
  design=matrix_descriptor(if(is.null(p$design))design else p$design),
  reduced=if(identical(p$test,'LRT'))matrix_descriptor(if(is.null(p$reduced))reduced else p$reduced)else NULL,
  size_factors=if(is.null(p$size))NULL else bits(p$size),normalization_factors=matrix_descriptor(p$factors),weights=matrix_descriptor(p$weights),
  fit_type=if(is.null(p$fit_type))'mean'else p$fit_type,test=if(is.null(p$test))'Wald'else p$test,
  threads=if(is.null(p$threads))1L else p$threads,min_replicates=if(is.null(p$min_replicates))0L else p$min_replicates,
  replacement_cutoff=bits(rsdeseq2:::.rsdeseq2_default_cooks_cutoff(if(is.null(p$design))design else p$design))[[1L]],
  queries=query_cases(p))
})
input<-file.path(out,'input.json');jsonlite::write_json(list(fresh_fit_cases=cases),input,auto_unbox=TRUE,null='null',digits=NA)
output<-system2(driver,shQuote(normalizePath(input)),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(output,'status')))stop(paste(output,collapse='\n'))
writeLines(output,file.path(out,'standalone.json'));standalone<-jsonlite::fromJSON(paste(output,collapse='\n'),simplifyVector=FALSE)$cases
records<-lapply(seq_along(cases),function(i){
 c<-cases[[i]];p<-profiles[[c$name]]
 fit<-tryCatch(blocked(rsdeseq2::DESeq(list(counts=matrix_input(c$counts),design=matrix_input(c$design)),
  fitType=c$fit_type,test=c$test,reduced=matrix_input(c$reduced),sizeFactors=if(is.null(c$size_factors))NULL else decode(c$size_factors),
  normalizationFactors=matrix_input(c$normalization_factors),weights=matrix_input(c$weights),
  nproc=c$threads,minReplicatesForReplace=if(c$min_replicates==0L)Inf else c$min_replicates)),error=identity)
 expected<-standalone[[i]]
 if(!is.null(p$error)||inherits(fit,'error')||!is.null(expected$error)){
  category<-if(!is.null(p$error))p$error else p$allow_error
  if(is.null(category))stop(c$name,': unexpected fitting failure')
  stopifnot(inherits(fit,'error'),!is.null(expected$error),grepl(category,conditionMessage(fit),ignore.case=TRUE),grepl(category,expected$error,ignore.case=TRUE))
  return(list(name=c$name,rejected=TRUE,native_error=expected$error,wrapper_error=conditionMessage(fit)))
 }
 if(inherits(fit,'error'))stop(c$name,': ',conditionMessage(fit))
 if(!is.null(expected$error))stop(c$name,': ',expected$error)
 stopifnot(nrow(fit$results)==c$counts$rows,is.null(rownames(fit$beta)),
  is.null(dimnames(fit$betaCovariance)[[1L]]),
  identical(rownames(fit$results),as.character(seq_len(c$counts$rows))))
 actual<-encode(fit);actual$queries<-query_values(fit,c$queries);wanted<-expected$value
 # Plain JSON canonicalization normalizes integer/logical/list containers;
 # floating-point values are all raw 64-bit encodings, including missing values.
 actual<-jsonlite::fromJSON(jsonlite::toJSON(actual,auto_unbox=TRUE,null='null'),simplifyVector=FALSE)
 jsonlite::write_json(actual,file.path(out,paste0(c$name,'-wrapped.json')),auto_unbox=TRUE,null='null',pretty=TRUE)
 if(!identical(canonical(actual),canonical(wanted)))
  stop(c$name,': ',paste(all.equal(canonical(actual),canonical(wanted),tolerance=0),collapse='\n'))
 saveRDS(fit,file.path(out,paste0(c$name,'.rds')))
 if(startsWith(c$name,'replacement_'))stopifnot(fit$nRefit>0L,any(fit$refitRows),any(fit$newAllZeroRows))
 if(startsWith(c$name,'newzero_'))stopifnot(fit$nRefit>0L,!any(fit$refitRows),any(fit$newAllZeroRows))
 if(grepl('no_refit',c$name,fixed=TRUE))stopifnot(fit$nRefit==0L)
 if(grepl('unreplaceable',c$name,fixed=TRUE))stopifnot(fit$nRefit==0L,!any(fit$replaceable),all(is.na(fit$replace)))
 if(grepl('weighted_failed',c$name,fixed=TRUE))stopifnot(fit$allZero[[39L]] == (fit$nRefit == 0L),
  is.na(fit$rawStat[[39L]]),is.na(fit$rawPvalue[[39L]]),is.na(fit$results$converged[[39L]]))
 list(name=c$name,value=actual)
})
names(records)<-names(profiles)
stopifnot(identical(canonical(records$mean_estimated$value),canonical(records$multiworker$value)))
names(records)<-NULL
jsonlite::write_json(list(cases=records),file.path(out,'wrapped.json'),auto_unbox=TRUE,null='null',pretty=TRUE)
rejections<-sum(vapply(records,function(c)isTRUE(c$rejected),logical(1)))
later_queries<-sum(vapply(records,function(r)if(is.null(r$value))0L else length(r$value$queries),integer(1)))
jsonlite::write_json(list(complete_fresh_fit_cases=length(cases),valid=length(cases)-rejections,
 later_queries=later_queries,
 controlled_rejections=rejections,R_calculations_blocked=TRUE,
 floating_point_bits_including_missing_values='exact'),file.path(out,'validation.json'),auto_unbox=TRUE,pretty=TRUE)
cat(length(cases),'complete fresh fits and',later_queries,'later queries: every numerical output bit and acceptance decision matches; R calculations blocked\n')
