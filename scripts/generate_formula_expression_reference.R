#!/usr/bin/env Rscript
# Independent built-in expression values and metadata; never loads rsdeseq2.
stopifnot(getRversion()=="4.6.2")
args<-commandArgs(TRUE)
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub("~+~"," ",script,fixed=TRUE)),".."))
paths<-if(length(args))args else file.path(root,"crates/rsdeseq2/tests/data/formula_expression_r462.json")
if(any(file.exists(paths)))stop("Use fresh paths; references are never overwritten")
hex<-function(v)paste(format(writeBin(as.double(v),raw(),size=8,endian="big")),collapse="")
strings<-function(v)lapply(v,function(x)if(is.na(x))NULL else x)
value<-function(v) {
  type<-if(is.null(v))"null"else if(typeof(v)%in%c("double","integer"))"number"else if(is.logical(v))"logical"else if(is.character(v))"text"else if(is.list(v))"list"else stop("Unsupported reference value")
  flat<-if(is.matrix(v))as.vector(t(v))else v
  data<-switch(type,number=lapply(unname(flat),hex),logical=lapply(unname(flat),function(x)if(is.na(x))NULL else x),text=strings(unname(flat)),list=lapply(unname(flat),value),NULL)
  keep<-c("class","levels","contrasts","scaled:center","scaled:scale","coefs","degree","units","tzone")
  attrs<-attributes(v);attrs<-attrs[names(attrs)%in%keep]
  list(type=type,integer=typeof(v)=="integer",values=data,dimensions=if(is.matrix(v))as.list(dim(v))else NULL,
    names=if(is.null(names(v)))NULL else as.list(names(v)),dimnames=if(is.matrix(v))lapply(dimnames(v),function(x)if(is.null(x))NULL else as.list(x))else list(NULL,NULL),attributes=lapply(attrs,value))
}
expression<-function(e) {
  if(is.symbol(e))return(list(type="symbol",name=as.character(e)))
  if(!is.call(e))return(list(type="literal",value=value(e)))
  call<-as.list(e);head<-call[[1L]]
  name<-if(is.symbol(head))as.character(head)else if(is.call(head)&&identical(head[[1L]],as.name("::")))paste0(as.character(head[[2L]]),"::",as.character(head[[3L]]))else stop("custom call head")
  inputs<-call[-1L];n<-names(inputs)
  list(type="call",function_name=name,arguments=lapply(seq_along(inputs),function(i)list(name=if(is.null(n)||n[[i]]=="")NULL else n[[i]],expression=expression(inputs[[i]]))))
}
x<-c(-2,-1,-0,0,.25,.5,1,2,3,5,7,9);y<-seq(.25,3,length.out=12)
f<-factor(rep(c("C","a","B"),4),levels=c("unused","B","C","a"))
ch<-c("a","A","á","ä","b","B","10","2","a-b","a_b","é","e")
M<-cbind(first=x,second=y);rownames(M)<-paste0("s",seq_len(nrow(M)))
z<-x;z[c(3L,5L)]<-c(NA_real_,NaN)
vars<-list(x=x,y=y,f=f,ch=ch,M=M,z=z,l=rep(c(TRUE,FALSE,NA),4),"TRUE"=TRUE,"FALSE"=FALSE,
  "NA"=NA,"NA_real_"=NA_real_,"NA_integer_"=NA_integer_,"NA_character_"=NA_character_,"NaN"=NaN,"Inf"=Inf,"NULL"=NULL,pi=pi,
  coefficients=attr(stats::poly(y,3),"coefs"))
specs<-alist(identity=x,negate=-x,plus=x+1,minus=1-x,product=x*y,quotient=x/y,power=I(x^3),
  fractional_power=y^.5,modulo=x%%1.5,integer_division=x%/%1.5,comparison=x>=y,
  logical_and=l&(x>0),logical_or=l|(x>0),logical_not=!l,
  logarithm=log(y),base_two=log(y,2),base_ten=log(y,base=10),base_three=log(y,base=3),
  log2=log2(y),log10=log10(y),log1p=log1p(x),sqrt=sqrt(y),exp=exp(x),expm1=expm1(y),
  abs=abs(x),sign=sign(x),floor=floor(y),ceiling=ceiling(y),trunc=trunc(x),
  numeric_factor=as.numeric(f),integer=as.integer(x),integer_missing=as.integer(c(-1.1,NA_real_,NaN,1e10)),
  numeric_matrix=as.numeric(M),text_numeric=as.character(c(0,1.25,1000,1e5,1e-4,1e-7)),text_factor=as.character(f),
  logical_numeric=as.logical(z),missing=is.na(z),nan=is.nan(z),finite=is.finite(z),infinite=is.infinite(z),
  text_missing=is.na(ch),conditional=ifelse(is.na(z),0,z),logical_conditional=ifelse(l,x,y),
  integer_conditional=ifelse(c(TRUE,FALSE),1L,2L),
  concatenate=c(x,y),list=list(alpha=c(1,2),norm2=c(1,3,4)),bind=cbind(x,y),bind_named=cbind(dose=x,y),
  bind_scalar=cbind(x,2),bind_matrix=cbind(M,x),matrix_default=matrix(x,nrow=4),matrix_byrow=matrix(x,4,3,TRUE),
  matrix_names=matrix(x,4,3,dimnames=list(c("a","b","c","d"),c("one","two","three"))),
  replicate=rep(c("A","B"),each=3),rep_times=rep(c(1,2),times=c(2,3)),rep_length=rep(x,length.out=17),
  seq=seq_len(4),length=length(x),scale=scale(x),scale_nested=scale(log(y)),scale_matrix=scale(M),
  scale_explicit=scale(M,center=c(1,2),scale=c(-2,4)),scale_positional=scale(x,FALSE,TRUE),
  scale_mixed=scale(x=x,FALSE),poly=poly(y,3),poly_raw=poly(x,3,raw=TRUE),poly_simple=poly(x,2,simple=TRUE),
  poly_reuse=poly(y,3,coefs=coefficients),
  factor_existing=factor(f),factor_numeric=factor(x),factor_text=factor(ch),factor_logical=factor(l),
  factor_named=factor(f,levels=c("a","C","B","unused"),labels=c("early","middle","late","absent")),
  factor_collapse=factor(f,labels=c("zero","same","same","other")),factor_prefix=factor(f,labels="group"),
  as_factor=as.factor(f),ordered=ordered(f),as_ordered=as.ordered(f),drop=droplevels(f),relevel=relevel(f,"C"),
  nested_factor=relevel(droplevels(factor(f)),ref="C"),
  I_matrix=I(cbind(x,y)*2),offset=offset(log(y)),sum=sum(x),mean=mean(x),sd=sd(y),
  mean_missing=mean(z),mean_remove_missing=mean(z,na.rm=TRUE),sum_missing=sum(z),sum_remove_missing=sum(z,na.rm=TRUE),
  sd_missing=sd(z),sd_remove_missing=sd(z,na.rm=TRUE),
  logical_text=as.logical(c("T","false","False","TRUE","oops",NA_character_)),
  missing_text=is.na(c("a",NA_character_)),nan_text=is.nan(ch),finite_text=is.finite(ch),
  infinity_text=is.infinite(ch),text_missing_numbers=as.character(z),
  integer_overflow=c(2147483647L,NA_integer_)+1L,
  matrix_predicate=is.na(M),logical_matrix=as.logical(M),
  scale_missing=scale(z),scale_constant=scale(rep(1,12)),
  scale_repair=ifelse(is.na(scale(rep(1,12))),0,scale(rep(1,12))),
  factor_missing=factor(z),factor_na_level=factor(c("a",NA_character_),exclude=NULL),
  conditional_unused_integer=ifelse(c(TRUE,TRUE),1L,2.5),conditional_all_missing=ifelse(c(NA,NA),1,2),
  conditional_text=ifelse(x>0,"positive","other"),
  invalid_custom=custom_function(x),invalid_nested=log(custom_function(x)),
  invalid_namespace=base::poly(x,2),invalid_duplicate=scale(x,center=TRUE,center=FALSE),
  invalid_ambiguous=poly(x,c=coefficients),invalid_matrix=cbind(M,matrix(x,4,3)),
  invalid_scale_length=scale(M,center=1),invalid_factor_levels=factor(ch,levels=c("a","a")),
  invalid_relevel=relevel(f,"absent"),invalid_poly_degree=poly(x,20),invalid_poly_missing=poly(z,2),
  conditional_matrix=ifelse(M>1,seq_len(5),seq_len(7)),
  conditional_matrix_branch=ifelse(M>1,M,M+2),
  conditional_unused_error=ifelse(TRUE,1,poly(z,2)),
  conditional_scalar_missing=ifelse(NA,poly(z,2),poly(z,2)),
  conditional_only_false=ifelse(c(FALSE,FALSE),1.5,2L),
  conditional_logical=ifelse(l,TRUE,FALSE),conditional_empty_branch=ifelse(c(TRUE,FALSE),NULL,2),
  invalid_matrix_names=matrix(x,4,3,dimnames=list(c("a"),c("x"))),
  invalid_nested_matrix_names=cbind(matrix(x,4,3,dimnames=list(c("a"),c("x"))),x))
cases<-lapply(names(specs),function(name) {
  expr<-specs[[name]]
  result<-tryCatch(withCallingHandlers(eval(expr,vars,enclos=asNamespace("stats")),warning=function(w)invokeRestart("muffleWarning")),error=identity)
  list(name=name,expression=expression(expr),expected=if(inherits(result,"error"))NULL else tryCatch(value(result),error=function(e)stop(name,": ",typeof(result),": ",conditionMessage(e))))
})
reference<-list(metadata=list(R=as.character(getRversion()),platform=R.version$platform,BLAS=unname(extSoftVersion()["BLAS"]),blas_core=Sys.getenv("OPENBLAS_CORETYPE"),ICU=unname(extSoftVersion()["ICU"]),collation=icuGetCollate("actual")),variables=lapply(vars,value),cases=cases)
for(path in paths){dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);jsonlite::write_json(reference,path,pretty=TRUE,auto_unbox=TRUE,null="null",digits=NA);cat("Wrote",length(cases),"cases to",path,"\n")}
