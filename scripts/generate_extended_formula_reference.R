#!/usr/bin/env Rscript
# Independent R calculations; the helper supplies primitive encoding only.
stopifnot(getRversion()=="4.6.2",packageVersion("DESeq2")=="1.52.0")
args<-commandArgs(TRUE)
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub("~+~"," ",script,fixed=TRUE)),".."))
source(file.path(root,"r-pkg/rsdeseq2/tests/testthat/helper-formula.R"))
paths<-if(length(args))args else file.path(root,"crates/rsdeseq2/tests/data/formula_expression_extended_r462.json")
if(any(file.exists(paths)))stop("Use fresh paths; references are never overwritten")
ast<-function(e){
 if(is.symbol(e))return(list(type="symbol",name=as.character(e)))
 if(!is.call(e))return(list(type="literal",value=.formula_descriptor(e)))
 h<-e[[1L]];name<-if(is.symbol(h))as.character(h)else paste0(as.character(h[[2L]]),"::",as.character(h[[3L]]))
 a<-as.list(e)[-1L];n<-names(a)
 list(type="call",function_name=name,arguments=lapply(seq_along(a),function(i)list(name=if(is.null(n)||n[[i]]=="")NULL else n[[i]],expression=ast(a[[i]]))))
}
x<-c(-2,-1,-0,0,.25,.5,1,2,3,5,7,9);y<-seq(.25,3,length.out=12);xi<-seq_len(12)
z<-x;z[c(3,5,7,9)]<-c(NA_real_,NaN,Inf,-Inf)
g<-factor(rep(c("A","C","B"),4),levels=c("C","B","A"));orderedg<-ordered(g)
ch<-c("A","a","ä","b","x","X","a-b","a_b","2","10",NA_character_,"é")
M<-cbind(x=x,y=y);rownames(M)<-paste0("s",seq_len(nrow(M)))
NM<-matrix(seq_len(6),2,3);names(NM)<-letters[1:6]
coefficients<-attr(stats::poly(x,2),"coefs");multiCoefficients<-attr(stats::poly(M,degree=2),"coefs")
vars<-list(x=x,y=y,xi=xi,z=z,g=g,orderedg=orderedg,ch=ch,M=M,NM=NM,coefficients=coefficients,multiCoefficients=multiCoefficients)
specs<-alist(sum_integer=sum(xi),sum_logical=sum(xi>4),sum_empty=sum(),sum_null=sum(NULL),sum_empty_double=sum(as.numeric(NULL)),
 sum_variadic=sum(x,y),sum_named=sum(first=x,second=y),sum_missing=sum(z),sum_clean=sum(z,na.rm=TRUE),sum_overflow=sum(c(2147483647L,1L)),
 mean_integer=mean(xi),mean_trim=mean(x,trim=.2),mean_trim_positional=mean(x,.2),mean_median=mean(x,trim=.5),mean_large_trim=mean(x,trim=1),
 mean_trim_missing=mean(c(x,NA_real_),trim=.1),mean_trim_clean=mean(c(x,NA_real_),trim=.1,na.rm=TRUE),mean_named_flag=mean(x,na.rm=TRUE),
 integer_center=I(xi-sum(xi)),character_equal=ch=="a",character_less=ch<"b",mixed_character=xi=="2",factor_equal=g=="A",factor_not_equal=g!="C",
 factor_compare=g==factor(g,levels=c("A","B","C")),ordered_less=orderedg<"A",conditional_factor=ifelse(g=="A",x,y),
 raw_missing=poly(z,2,raw=TRUE),raw_repair=ifelse(is.na(poly(z,2,raw=TRUE)),0,poly(z,2,raw=TRUE)),recurrence_missing=poly(z,2,coefs=coefficients),
 polynomial_matrix=poly(M,degree=2),polynomial_matrix_positional=poly(M,2),polynomial_dots=poly(x,y,degree=2),
 polynomial_matrix_raw=poly(M,degree=2,raw=TRUE),polynomial_dots_raw=poly(x,y,degree=3,raw=TRUE),
 polynomial_matrix_simple=poly(M,degree=2,simple=TRUE),polynomial_recurrence=poly(M,degree=2,coefs=multiCoefficients),
 scale_named=scale(c(first=1,second=2,third=4)),scale_nested=scale(scale(x)),scale_asis=scale(I(M)),
 sum_nan_na=sum(c(NaN,NA_real_)),sum_na_nan=sum(c(NA_real_,NaN)),mean_nan_na=mean(c(NaN,NA_real_)),sd_nonfinite=sd(z),
 missing_condition=ifelse(is.na(z),0,z),invalid_degree=poly(M,degree=20),invalid_sum=sum(ch),invalid_trim=mean(x,trim=NA_real_),nested_asis=I(I(x)),
 named_matrix=NM,named_matrix_concat=c(NM),named_matrix_scale=scale(NM),named_matrix_add=NM+1,
 integer_neg_missing=-c(0L,1L,NA_integer_),integer_zero_product=c(0L,1L,NA_integer_)*(-1L),
 mean_null=mean(NULL),mean_text=mean(ch),concatenate_null=c(NULL),mean_concatenate_null=mean(c()))
cases<-lapply(names(specs),function(name){
 e<-specs[[name]];result<-tryCatch(withCallingHandlers(eval(e,vars,enclos=asNamespace("stats")),warning=function(w)invokeRestart("muffleWarning")),error=identity)
 list(name=name,expression=ast(e),expected=if(inherits(result,"error"))NULL else .formula_descriptor(result))
})
reference<-list(metadata=list(R=as.character(getRversion()),DESeq2=as.character(packageVersion("DESeq2")),collation=icuGetCollate("actual")),variables=lapply(vars,.formula_descriptor),cases=cases)
for(path in paths){dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);jsonlite::write_json(reference,path,pretty=TRUE,auto_unbox=TRUE,null="null",digits=NA);cat("Wrote",length(cases),"independent extended cases to",path,"\n")}
