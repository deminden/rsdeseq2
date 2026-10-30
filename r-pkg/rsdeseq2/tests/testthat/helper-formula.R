# Test-only exact primitive metadata encoder/decoder. Never used by production.
hex<-function(v)paste(format(writeBin(as.double(v),raw(),size=8,endian="big")),collapse="")
strings<-function(v)lapply(v,function(x)if(is.na(x))NULL else x)
.formula_descriptor<-function(v) {
  type<-if(is.null(v))"null"else if(typeof(v)%in%c("double","integer"))"number"else if(is.logical(v))"logical"else if(is.character(v))"text"else if(is.list(v))"list"else stop("Unsupported reference value")
  flat<-if(is.matrix(v))as.vector(t(v))else v
  elementNames<-names(v)
  if(is.matrix(v)&&!is.null(elementNames)){dim(elementNames)<-dim(v);elementNames<-as.vector(t(elementNames))}
  data<-switch(type,number=lapply(unname(flat),hex),logical=lapply(unname(flat),function(x)if(is.na(x))NULL else x),text=strings(unname(flat)),list=lapply(unname(flat),.formula_descriptor),NULL)
  keep<-c("class","levels","contrasts","scaled:center","scaled:scale","coefs","degree","units","tzone")
  attrs<-attributes(v);attrs<-attrs[names(attrs)%in%keep]; if(length(attrs))attrs<-attrs[sort(names(attrs))]
  list(type=type,integer=typeof(v)=="integer",values=data,dimensions=if(is.matrix(v))as.list(dim(v))else NULL,
    names=if(is.null(elementNames))NULL else as.list(elementNames),dimnames=if(is.matrix(v))lapply(dimnames(v),function(x)if(is.null(x))NULL else as.list(unname(x)))else list(NULL,NULL),dimname_names=if(is.matrix(v)&&!is.null(dimnames(v)))lapply(dimnames(v),function(x)if(is.null(names(x)))NULL else as.list(names(x)))else list(NULL,NULL),attributes=lapply(attrs,.formula_descriptor))
}

.formula_decode <- function(v) {
 rawDouble <- function(hex)readBin(as.raw(strtoi(substring(hex,seq(1,15,2),seq(2,16,2)),16L)),"double",n=1L,size=8,endian="big")
 x <- switch(v$type,number=vapply(v$values,rawDouble,0.0),logical=vapply(v$values,function(x)if(is.null(x))NA else x,TRUE),text=vapply(v$values,function(x)if(is.null(x))NA_character_ else x,""),list=lapply(v$values,.formula_decode),null=NULL)
 if(isTRUE(v$integer))x<-as.integer(x)
 if(length(v$dimensions)){d<-as.integer(unlist(v$dimensions));dim(x)<-rev(d);x<-t(x);dimnames(x)<-lapply(v$dimnames,unlist)}
 if(length(v$dimname_names)){
   dn<-dimnames(x)
   for(i in seq_along(v$dimname_names))if(!is.null(v$dimname_names[[i]]))names(dn[[i]])<-unlist(v$dimname_names[[i]])
   if(!is.null(dn))dimnames(x)<-dn
 }
 if(!is.null(v$names)){
  nm<-unlist(v$names)
  if(length(v$dimensions)){dim(nm)<-rev(as.integer(unlist(v$dimensions)));nm<-as.vector(t(nm))}
  names(x)<-nm
 }
 if(length(v$attributes))for(name in names(v$attributes))attr(x,name)<-.formula_decode(v$attributes[[name]])
 x
}
.formula_ast_decode <- function(v) {
 if(v$type=="literal")v["value"]<-list(.formula_decode(v$value))
 if(v$type=="call")v$arguments<-lapply(v$arguments,function(a){a$expression<-.formula_ast_decode(a$expression);a})
 v
}
.formula_canonical <- function(x)jsonlite::fromJSON(jsonlite::toJSON(x,auto_unbox=TRUE,null="null",digits=NA),simplifyVector=FALSE)
