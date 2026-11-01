#!/usr/bin/env Rscript
# Independent model.matrix references; never loads rsdeseq2.
stopifnot(getRversion()=="4.6.2")
args <- commandArgs(TRUE)
allLevels <- "--all-levels" %in% args
args <- args[args!="--all-levels"]
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~"," ",script,fixed=TRUE)),".."))
paths <- if(length(args))args else file.path(root,c("crates/rsdeseq2/tests/data/model_matrix_r462.json","r-pkg/rsdeseq2/tests/testthat/data/model_matrix_r462.json"))
if(any(file.exists(paths)))stop("Use fresh paths; references are never overwritten")
hex <- function(v)paste(format(writeBin(as.double(v),raw(),size=8,endian="big")),collapse="")
hexes <- function(v)lapply(unname(v),hex)
d <- data.frame(a=factor(rep(c("A","B","C"),4),levels=c("C","A","B","unused")),
 b=factor(rep(c("x","y"),each=6)), c=factor(rep(c("i","j"),6)),
 x=c(-2,-1,-0,0,.5,1,2,3,4,5,6,7), y=seq(1,3,length.out=12),
 flag=rep(c(TRUE,FALSE),6), text=rep(c("z","p","q"),4), check.names=FALSE)
d[["cell type"]]<-d$a;d[["dose value"]]<-d$x
d$date<-as.Date("2020-01-01")+seq_len(nrow(d))
d$single<-factor(rep("only",nrow(d)))
specs<-list(intercept=~1,empty=~0,additive=~a+b,reverse=~x+a,
 factor_interaction=~a*b,numeric_interaction=~x*y,factor_numeric=~x*a,
 interaction_only=~a:b,no_intercept=~0+a+b,no_intercept_interaction=~0+a*b,
 nested=~a/b,third_order=~a*b*c,third_numeric=~a*b*x,
 without_mains=~x:a,grouped=~(a+b+x)^2,removed=~a*b-a,
 logical=~flag*a,character=~text+x,raw_poly=~poly(x,3,raw=TRUE),
 orthogonal_poly=~poly(x,3),numeric_matrix=~I(cbind(x,y)),
 matrix_interaction=~I(cbind(x,y)):a,single_matrix=~I(cbind(x)),
 unnamed_matrix=~I(cbind(x,y)*2),transform=~log1p(y)+I(x^2),
 offset=~a+offset(log(y)),backticks=~`cell type`*`dose value`,
 sum=~a*b,helmert=~a*x,sas=~a+b,custom=~a*x,ordered=~a*x)
specs$date<-~date+x
specs$response<-y~a
specs$response_factor<-single~x
specs$response_rhs<-y~y*x
levelCounts<-if(allLevels)2:95 else c(2L,4L,8L,16L,24L,32L,53L,64L,68L,85L,91L,95L)
for(n in levelCounts)specs[[paste0("ordered_",n)]]<-~a
cases<-lapply(names(specs),function(name) {
 data<-d
 if(startsWith(name,"ordered_"))data<-data.frame(a=ordered(seq_len(as.integer(sub("ordered_","",name,fixed=TRUE)))))
 if(name=="sum")attr(data$a,"contrasts")<-"contr.sum"
 if(name=="helmert")attr(data$a,"contrasts")<-"contr.helmert"
 if(name=="sas")attr(data$a,"contrasts")<-"contr.SAS"
 if(name=="custom")attr(data$a,"contrasts")<-matrix(c(1,.5,-.25,0,-.5,1,0,.25),4,2,dimnames=list(levels(data$a),c("early","late")))
 if(name=="ordered")data$a<-ordered(data$a)
 f<-specs[[name]];frame<-model.frame(f,data);t<-terms(frame);codes<-attr(t,"factors")
 if(!length(codes))codes<-matrix(integer(),length(frame),0)
 variables<-lapply(seq_along(frame),function(i) {
  v<-frame[[i]];label<-if(is.null(rownames(codes)))names(frame)[[i]]else rownames(codes)[[i]]
  if(is.character(v))v<-factor(v)
  if(i==attr(t,"response"))v<-as.double(v)
  if(is.logical(v))v<-factor(v,levels=c(FALSE,TRUE))
  if(is.factor(v)) {
   coding<-attr(v,"contrasts")
   if(is.null(coding))coding<-if(is.ordered(v))"contr.poly" else "contr.treatment"
   list(name=label,raw_name=names(frame)[[i]],type="factor",codes=as.list(as.integer(v)),levels=as.list(levels(v)),
    coding=if(is.matrix(coding))"matrix"else coding,contrast_columns=if(is.matrix(coding))ncol(coding)else NULL,
    contrast_values=if(is.matrix(coding))hexes(as.vector(t(coding)))else NULL,
    contrast_names=if(is.matrix(coding)&&!is.null(colnames(coding)))as.list(colnames(coding))else NULL)
  } else {
   if(!is.matrix(v))v<-matrix(v,ncol=1)
   list(name=label,raw_name=names(frame)[[i]],type="numeric",columns=ncol(v),values=hexes(as.vector(t(v))),
    column_names=if(is.null(colnames(v)))NULL else as.list(colnames(v)))
  }
 })
 expected<-suppressWarnings(model.matrix(t,frame))
 list(name=name,samples=nrow(frame),variables=variables,terms=lapply(seq_len(ncol(codes)),function(j)as.list(unname(codes[,j]))),
  intercept=attr(t,"intercept")==1L,response=attr(t,"response"),expected_names=as.list(colnames(expected)),assignments=as.list(attr(expected,"assign")),
  expected=hexes(as.vector(t(expected))),tolerance=0)
})
averageSpecs<-list(two_vs_one=list(p=c(0L,2L),n=5L),thirds_vs_pair=list(p=c(0L,2L,3L),n=c(4L,6L)),
 positive_only=list(p=c(1L,3L,6L),n=integer()),negative_only=list(p=integer(),n=c(0L,2L,4L)),
 coefficient=list(p=5L,n=integer()),pair=list(p=4L,n=1L),empty=list(p=integer(),n=integer()),
 duplicate=list(p=c(0L,0L),n=5L),overlap=list(p=0L,n=0L),outside=list(p=8L,n=integer()))
averages<-lapply(names(averageSpecs),function(name){x<-averageSpecs[[name]];selected<-c(x$p,x$n)
 valid<-length(selected)>0&&!anyDuplicated(selected)&&all(selected>=0L&selected<7L)
 value<-numeric(7)
 if(valid){if(length(x$p))value[x$p+1L]<-1/length(x$p);if(length(x$n))value[x$n+1L]<- -1/length(x$n)}
 list(name=name,coefficients=7L,positive=as.list(x$p),negative=as.list(x$n),expected=if(valid)hexes(value)else NULL)
})
ref<-list(metadata=list(R=as.character(getRversion()),platform=R.version$platform,BLAS=unname(extSoftVersion()["BLAS"]),blas_core=Sys.getenv("OPENBLAS_CORETYPE"),
 purpose="Independent stats::model.matrix assembly and expanded list-contrast means; expression evaluation tested separately"),cases=cases,average_contrasts=averages)
for(path in paths){dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);jsonlite::write_json(ref,path,pretty=TRUE,auto_unbox=TRUE,null="null",digits=NA);cat("Wrote",path,"\n")}
