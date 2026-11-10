#!/usr/bin/env Rscript
# Independent original R tail reference; never loads rsdeseq2.
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: generate_normal_log_reference.R OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2')
numbers<-function(x)as.list(sprintf('%.17g',as.numeric(x)))
extra<-c(0,1e-15,1e-8,8-1e-14,8,8+1e-14,37.5,37.5193,38,40,100,1e4,1e6,1e154,1.4e154,1e155,Inf)
x<-sort(unique(c(seq(-40,40,by=.125),extra,-extra)))
x<-c(x,NaN)
jsonlite::write_json(list(source=list(R=as.character(getRversion()),function_name='stats::pnorm(log.p=TRUE)'),
 arguments=numbers(x),log_cdf=numbers(stats::pnorm(x,log.p=TRUE)),log_upper=numbers(stats::pnorm(x,lower.tail=FALSE,log.p=TRUE))),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
