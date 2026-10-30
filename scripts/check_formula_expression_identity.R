#!/usr/bin/env Rscript
# Independent reference descriptors, direct standalone output, and actual bridge.
args<-commandArgs(TRUE)
if(length(args)!=3L)stop("usage: STANDALONE_DRIVER R_LIBRARY FRESH_OUTPUT_DIR")
driver<-normalizePath(args[[1L]]);library(rsdeseq2,lib.loc=normalizePath(args[[2L]]))
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[[1L]])
root<-normalizePath(file.path(dirname(gsub("~+~"," ",script,fixed=TRUE)),".."))
source(file.path(root,"r-pkg/rsdeseq2/tests/testthat/helper-formula.R"))
out<-args[[3L]]
if(dir.exists(out))stop("Use a fresh identity output directory")
dir.create(out,recursive=TRUE)
canonical<-function(x){
 if(!is.list(x))return(x)
 if(!is.null(names(x)))x<-x[sort(names(x))]
 lapply(x,canonical)
}
total<-0L
for(file in c("formula_expression_r462.json","formula_expression_extended_r462.json")) {
 panelOut<-file.path(out,sub(".json$","",file));dir.create(panelOut)
reference<-jsonlite::fromJSON(file.path(root,"crates/rsdeseq2/tests/data",file),simplifyVector=FALSE)
input<-list(formula_cases=reference$cases,variables=reference$variables,collation=reference$metadata$collation)
path<-file.path(panelOut,"input.json");jsonlite::write_json(input,path,auto_unbox=TRUE,null="null",digits=NA)
output<-system2(driver,shQuote(normalizePath(path)),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(output,"status")))stop(paste(output,collapse="\n"))
standalone<-jsonlite::fromJSON(paste(output,collapse="\n"),simplifyVector=FALSE)
variables<-lapply(reference$variables,.formula_decode)
wrapped<-lapply(reference$cases,function(case){
 request<-list(operation="evaluate",expressions=list(.formula_ast_decode(case$expression)),variables=variables,collation=reference$metadata$collation)
 result<-tryCatch(rsdeseq2:::.rsdeseq2_formula_request(request),error=identity)
 if(inherits(result,"error"))return(list(name=case$name,ok=FALSE))
 list(name=case$name,ok=TRUE,value=.formula_descriptor(result[[1L]]))
})
wrapped<-.formula_canonical(list(cases=wrapped))
jsonlite::write_json(wrapped,file.path(panelOut,"wrapped.json"),auto_unbox=TRUE,null="null",pretty=TRUE)
writeLines(output,file.path(panelOut,"standalone.json"))
if(!identical(canonical(wrapped),canonical(standalone))) {
 for(i in seq_along(wrapped$cases)) if(!identical(canonical(wrapped$cases[[i]]),canonical(standalone$cases[[i]]))) {
   stop(wrapped$cases[[i]]$name,": ",paste(all.equal(canonical(wrapped$cases[[i]]),canonical(standalone$cases[[i]]),tolerance=0),collapse="; "))
 }
 stop("formula panel envelope differs")
}
for(i in seq_along(wrapped$cases)){
 case<-wrapped$cases[[i]];expected<-reference$cases[[i]]$expected
 stopifnot(identical(case$ok,!is.null(expected)))
 if(case$ok)stopifnot(identical(canonical(case$value),canonical(.formula_canonical(.formula_descriptor(.formula_decode(expected))))))
}
cat(length(reference$cases),"direct formula expression cases: every bit and primitive metadata field identical\n")
 total<-total+length(reference$cases)
}
cat(total,"total direct expression cases identical\n")
