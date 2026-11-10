#!/usr/bin/env Rscript
# Independent normalization leaves used to localize conditioned moment errors.
args <- commandArgs(TRUE)
if(length(args)!=1L)stop('usage: generate_student_component_leaves.R FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',!file.exists(args[[1L]]))
bits <- function(x)lapply(as.numeric(x),function(v)paste0('bits:',paste(format(writeBin(v,raw(),size=8,endian='big')),collapse='')))
numbers <- function(x)as.list(sprintf('%.17g',as.numeric(x)))
gamma_arguments <- c(.25,.5,.75,1,1.5,4,4.5,5,5.5,10,10.5,49,49.5,50,50.5,51,100)
cases <- list()
for(df in c(10,100))for(q in c(-150,-2.5e-7,2.5e-7,2.5))for(order in 0:2){
 a <- .5*(df-order);b <- .5*(order+1);p <- df/(df+q*q)
 cases[[length(cases)+1L]] <- list(df=bits(df)[[1L]],q=bits(q)[[1L]],order=order,a=bits(a)[[1L]],b=bits(b)[[1L]],p=bits(p)[[1L]],
  beta=numbers(beta(a,b))[[1L]],regularized=numbers(pbeta(p,a,b))[[1L]],
  G=numbers(etrunct:::G(order,df,q))[[1L]])
}
jsonlite::write_json(list(source=list(R=as.character(getRversion()),etrunct=as.character(packageVersion('etrunct')),input_encoding='f64_bits'),
 gamma_arguments=bits(gamma_arguments),gamma_values=numbers(gamma(gamma_arguments)),cases=cases),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null='null')
