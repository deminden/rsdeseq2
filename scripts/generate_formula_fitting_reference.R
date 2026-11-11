#!/usr/bin/env Rscript
# Original R model.frame/model.matrix only; never imports candidate calculations.
args <- commandArgs(TRUE)
if(length(args)!=1L)stop("usage: FRESH_OUTPUT_JSON")
stopifnot(as.character(getRversion())=="4.6.2",!file.exists(args[[1L]]))
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub("~+~"," ",script,fixed=TRUE)),".."))
source(file.path(root,"r-pkg/rsdeseq2/tests/testthat/helper-formula.R"))
bits<-function(v)lapply(as.numeric(v),function(x)paste0("bits:",hex(x)))
d<-data.frame(condition=factor(rep(c("A","B","C"),each=4L)),
  batch=factor(rep(c("N","N","Y","Y"),3L)),exposure=rep(1:4,3L),
  u=c(.2,.6,.1,.8,.7,.3,.9,.4,.5,.1,.6,.9))
formulas<-c(
 removed_alias="~I(exposure)+I(+exposure)-I(exposure)",
 alias_interaction="~I(exposure)+condition:I(+exposure)",
 removed_polynomial_alias="~poly(exposure,2)+poly(exposure,degree=2)-poly(exposure,2)",
 numeric_star="~exposure*condition", mixed_star="~exposure*condition*batch",
 mixed_power="~(exposure+condition+u+batch)^3",
 numeric_nested="~exposure/condition/batch",mixed_nested="~condition/exposure/batch",
 numeric_in="~exposure%in%condition",mixed_group="~(exposure+condition)*batch",
 interaction_subsets="~condition:exposure:batch+u:exposure+batch:u",
 factor_readded="~condition*batch-condition+condition",
 numeric_first="~exposure+condition+batch",numeric_between="~condition+exposure+batch",
 removed_readded="~condition+exposure-condition+condition",
 interaction_order="~condition:exposure+batch:exposure+condition+batch",
 numeric_interaction_first="~exposure:condition+condition:batch",
 no_intercept_numeric_first="~0+exposure+condition+batch",
 polynomial_main="~poly(exposure,2)+condition",
 polynomial_first_interaction="~poly(exposure,2):condition",
 factor_first_polynomial="~condition:poly(exposure,2)",
 polynomial_main_interaction="~poly(exposure,2)*condition",
 polynomial_power="~(poly(exposure,2)+condition)^2",
 grouped_numeric_interactions="~(poly(exposure,2)+u)*condition",
 ordered_numeric_first="~exposure+ordered(condition)",
 dot="~.",dot_power="~.^2",subtracted_dot="~.-batch",
 repeated_polynomial="~poly(exposure,2)+poly(exposure, 2)",
 repeated_identity="~I(exposure)+I( exposure)",
 offsets="~exposure+condition+offset(log(exposure))",
 two_offsets="~condition+offset(log(exposure))+offset(u)")
cases<-lapply(names(formulas),function(name){
 f<-as.formula(formulas[[name]])
 frame<-stats::model.frame(f,d,na.action=stats::na.pass)
 matrix<-stats::model.matrix(f,frame)
 offsets<-stats::model.offset(frame)
 list(name=name,formula=formulas[[name]],columns=ncol(matrix),names=as.list(colnames(matrix)),
  values=bits(as.vector(t(matrix))),offsets=bits(offsets))
})
factors<-lapply(c("condition","batch"),function(name)list(name=name,sample_levels=as.list(as.character(d[[name]])),
 levels=as.list(levels(d[[name]])),reference=levels(d[[name]])[[1L]]))
numeric<-lapply(c("exposure","u"),function(name)list(name=name,values=bits(d[[name]])))
jsonlite::write_json(list(source=list(R=as.character(getRversion())),samples=nrow(d),
 factors=factors,numeric=numeric,cases=cases),args[[1L]],auto_unbox=TRUE,pretty=TRUE,null="null",digits=NA)
