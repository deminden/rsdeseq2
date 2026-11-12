#!/usr/bin/env Rscript
# Original interpreted vector expressions. Byte-code scalar intrinsics may have
# different rounding; model-frame expressions use the vector operation here.
args <- commandArgs(TRUE)
if (length(args)!=1L) stop('usage: generate_elementary_expression_reference.R FRESH_OUTPUT_JSON')
stopifnot(as.character(getRversion())=='4.6.2',!file.exists(args[[1L]]))
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub('~+~',' ',script,fixed=TRUE)),'..', '..'))
source(file.path(root,'r-pkg/rsdeseq2/tests/testthat/helper-formula.R'))
eps <- .Machine$double.eps
x <- c(-0,0,eps*c(-2,-1.5,-1,-.5,.5,1,1.5,2),
       sqrt(eps)*c(-2,-1-2*eps,-1,-1+eps,-.5,.5,1-eps,1,1+2*eps,2),
       eps*c(-eps,eps),-.Machine$double.xmin,.Machine$double.xmin)
cases <- lapply(c('exp','log1p','expm1'),function(f) {
  e <- call(f,as.name('x'))
  list(name=paste0(f,'_boundaries'),
       expression=list(type='call',function_name=f,arguments=list(list(name=NULL,expression=list(type='symbol',name='x')))),
       expected=.formula_descriptor(eval(e,list(x=x))))
})
jsonlite::write_json(list(metadata=list(R=as.character(getRversion()),
  platform=R.version$platform,collation=icuGetCollate('actual')),
  variables=list(x=.formula_descriptor(x)),cases=cases),args[[1L]],
  auto_unbox=TRUE,pretty=TRUE,null='null',digits=NA)
