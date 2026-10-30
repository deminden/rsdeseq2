#!/usr/bin/env Rscript
# Independent original-R matrices, metadata and expanded-column maps.
stopifnot(getRversion()=="4.6.2",packageVersion("DESeq2")=="1.52.0")
args <- commandArgs(TRUE)
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~"," ",script,fixed=TRUE)),".."))
paths <- if(length(args))args else file.path(root,"crates/rsdeseq2/tests/data/native_factor_r462.json")
if(any(file.exists(paths)))stop("Use fresh paths; references are never overwritten")
hex <- function(v)paste(format(writeBin(as.double(v),raw(),size=8,endian="big")),collapse="")
hexes <- function(v)lapply(unname(v),hex)
d <- data.frame(x=c(-3.1,-.7,-0,0,.13,.7,1.1,2.3,4.7,-1.3,.001,-.001),
  g=factor(rep(c("A","B","C"),4),levels=c("unused","A","B","C")),
  h=factor(rep(c("X","X","Y","Y"),3),levels=c("X","Y")),
  q=factor(rep(c("C","A","B"),4),levels=c("C","A","B")))
numeric <- list(list(name="x",values=hexes(d$x)))
d$y <- c(.7,1.3,-2.1,.17,-.13,3.7,-7.3,1e-7,11.1,.03,2.9,-.9)
d$z <- c(1.1,-.3,2.7,1e-9,3.3,-1.7,.19,.003,7.7,-.23,1.9,-4.1)
numeric <- c(numeric,list(list(name="y",values=hexes(d$y)),list(name="z",values=hexes(d$z))))
d$r <- factor(rep(c("B","A"),6),levels=c("Z","B","A"))
contrasts(d$r) <- stats::contr.treatment(levels(d$r),base=3L)
factors <- lapply(c("g","h","q","r"),function(name)list(name=name,sample_levels=as.list(as.character(d[[name]])),levels=as.list(levels(d[[name]])),reference=if(name=="r")"A"else levels(d[[name]])[[1L]]))
cases <- list()
record <- function(name,formula,r_columns=NULL,r_expanded_columns=NULL,native_label=NULL,native_columns=NULL,native_expanded_columns=NULL) {
  f <- as.formula(formula)
  result <- tryCatch(suppressWarnings({
    frame <- stats::model.frame(f,d,na.action=stats::na.pass)
    standard <- stats::model.matrix(f,frame)
    names <- names(frame)[vapply(frame,is.factor,TRUE)]
    unordered <- names[!vapply(frame[names],is.ordered,TRUE)]
    treatment <- stats::setNames(lapply(frame[unordered],function(v){
      m<-stats::contr.treatment(levels(v));colnames(m)<-paste0("_",levels(v)[-1L],"_vs_",levels(v)[[1L]]);m
    }),unordered)
    identity <- stats::setNames(lapply(frame[unordered],function(v){
      m<-diag(nlevels(v));dimnames(m)<-list(levels(v),levels(v));m
    }),unordered)
    named <- stats::model.matrix(f,frame,contrasts.arg=if(length(treatment))treatment else NULL)
    stopifnot(identical(hexes(as.vector(t(named))),hexes(as.vector(t(standard)))))
    expanded <- stats::model.matrix(f,frame,contrasts.arg=if(length(identity))identity else NULL)
    if(is.null(r_columns))r_columns<-colnames(standard)
    if(is.null(r_expanded_columns))r_expanded_columns<-colnames(expanded)
    rename <- function(labels) {
      labels[labels=="(Intercept)"]<-"Intercept"
      if(!is.null(native_label))labels<-gsub(names(frame)[[1L]],native_label,labels,fixed=TRUE)
      labels
    }
    frameNames<-names
    if(!is.null(native_label))frameNames[frameNames==names(frame)[[1L]]]<-native_label
    frames<-lapply(seq_along(names),function(i){v<-frame[[names[[i]]]];list(name=frameNames[[i]],levels=as.list(levels(v)),reference=levels(v)[[1L]],sample_levels=as.list(as.character(v)),ordered=is.ordered(v))})
    indices<-match(r_columns,colnames(standard));expandedIndices<-match(r_expanded_columns,colnames(expanded))
    groups<-match(r_columns,r_expanded_columns)-1L
    stopifnot(!anyNA(indices),!anyNA(expandedIndices),!anyNA(groups))
    list(names=as.list(if(is.null(native_columns))rename(colnames(named)[indices])else native_columns),values=hexes(as.vector(t(standard[,indices,drop=FALSE]))),
      expanded_names=as.list(if(is.null(native_expanded_columns))rename(colnames(expanded)[expandedIndices])else native_expanded_columns),expanded_values=hexes(as.vector(t(expanded[,expandedIndices,drop=FALSE]))),
      groups=lapply(groups,function(i)as.list(i)),frames=frames,frame_names=as.list(frameNames))
  }),error=identity)
  valid<-!inherits(result,"error")&&all(vapply(result$values,function(v)is.finite(readBin(as.raw(strtoi(substring(v,seq(1,15,2),seq(2,16,2)),16L)),"double",1L,size=8,endian="big")),TRUE))
  cases[[length(cases)+1L]] <<- list(name=name,formula=formula,details=TRUE,
    frame_names=if(valid)result$frame_names else list(),r_columns=if(valid)as.list(r_columns)else NULL,
    r_expanded_columns=if(valid)as.list(r_expanded_columns)else NULL,
    native_names=if(valid)result$names else NULL,expected=if(valid)result$values else NULL,
    expanded_names=if(valid)result$expanded_names else NULL,expanded_expected=if(valid)result$expanded_values else NULL,
    groups=if(valid)result$groups else NULL,frames=if(valid)result$frames else NULL,
    reference_error=if(inherits(result,"error"))conditionMessage(result)else NULL)
}
for(expr in c("factor(g)","as.factor(g)","ordered(g)","as.ordered(g)","droplevels(g)",
  "factor(q)","ordered(q)","as.factor(ordered(g))","factor(as.ordered(g))","droplevels(as.ordered(g))",
  "ordered(relevel(g,ref='B'))","relevel(droplevels(g),ref='B')",
  "factor(g,levels=c('C','B','A'))","ordered(g,levels=c('C','B','A'))",
  "factor(g,labels=c('a','b','c'))","ordered(g,labels=c('a','b','c'))",
  "factor(g,labels=c('AB','AB','C'))","ordered(g,labels=c('AB','AB','C'))",
  "factor(g,labels=c('P'))","ordered(g,labels=c('P'))",
  "factor(g,levels=c('unused','A','B','C'))","ordered(g,levels=c('unused','A','B','C'))")) {
  record(expr,paste0("~",expr))
}
record("scalar_label","~factor(g,labels='P')",native_label='factor(g, labels = c("P"))')
for(expr in c("ordered(g)","as.factor(g)","factor(g)","ordered(q)"))record(paste0("no_intercept_",expr),paste0("~0+",expr))
for(f in c("~ordered(g)+x","~ordered(g)*x","~ordered(g):x","~ordered(g)+h+ordered(g):h",
  "~ordered(g):h","~0+ordered(g):h","~ordered(g)+ordered(h)","~ordered(g)*ordered(h)",
  "~ordered(g)+h:x","~ordered(g):h:x","~ordered(g)+h+x+ordered(g):h:x",
  "~0+ordered(g)+h","~0+h+ordered(g)"))record(f,f)
record("numeric_first","~x+ordered(g)",r_columns=c("(Intercept)","ordered(g).L","ordered(g).Q","x"),r_expanded_columns=c("(Intercept)","ordered(g).L","ordered(g).Q","x"))
record("numeric_first_product","~x:ordered(g):ordered(h)")
for(f in c("~relevel(ordered(g),ref='B')","~relevel(as.factor(ordered(g)),ref='B')",
  "~factor(g,levels=c('A','B'))","~factor(g,labels=c('all','all','all'))",
  "~factor(g,labels=c('a','b'))"))record(f,f)
for(suffix in c("x","x:y")) {
  products<-unlist(lapply(levels(d$g),function(g)paste0("g",g,":h",levels(d$h),":",suffix)))
  main<-c("(Intercept)",paste0("g",levels(d$g)[-1L]),"hY","x",if(suffix=="x:y")"y")
  expandedMain<-c("(Intercept)",paste0("g",levels(d$g)),paste0("h",levels(d$h)),"x",if(suffix=="x:y")"y")
  record(paste0("unordered_higher_",suffix),paste0("~g+h+",if(suffix=="x")"x"else"x+y","+g:h:",suffix),r_columns=c(main,products),r_expanded_columns=c(expandedMain,products))
}
record("numeric_product_frame_order","~z+x+y+x:y:z",native_columns=c("Intercept","z","x","y","x:y:z"),native_expanded_columns=c("Intercept","z","x","y","x:y:z"))
record("droplevels_resets_contrasts","~droplevels(r)")
source(file.path(root,"r-pkg/rsdeseq2/tests/testthat/helper-formula.R"))
named <- d$q;names(named)<-paste0("sample",seq_along(named))
custom <- d$g;contrasts(custom)<-stats::contr.sum(nlevels(custom))
inputs <- list(unused=d$g,reordered=stats::relevel(d$q,"B"),named=named,custom_contrasts=custom,
  one_level=factor(c("A","A")),ordered=ordered(d$g),collision=factor(c("_null_level_","A")),empty=factor(character()),numeric=c(1,2))
expansions <- lapply(names(inputs),function(name) {
  v<-inputs[[name]]
  if(length(v))v<-v[c(seq_along(v),length(v))]
  expected<-tryCatch({
    if(!is.factor(v)||!length(v)||"_null_level_"%in%levels(v))stop("invalid expanded-factor input")
    levels(v)<-c(levels(v),"_null_level_")
    v<-stats::relevel(v,"_null_level_")
    v[length(v)]<-"_null_level_"
    .formula_descriptor(v)
  },error=function(e)NULL)
  input<-inputs[[name]];if(length(input))input<-input[c(seq_along(input),length(input))]
  list(name=name,input=.formula_descriptor(input),null_level="_null_level_",expected=expected)
})
output<-list(metadata=list(R=as.character(getRversion()),DESeq2=as.character(packageVersion("DESeq2")),expanded_policy="Original model.matrix with identity contrasts for unordered factors; ordered polynomial coordinates retained; not DESeq2 expanded ordered-formula betaPrior"),numeric=numeric,factors=factors,native_formula_cases=cases,factor_expansions=expansions)
for(path in paths)jsonlite::write_json(output,path,auto_unbox=TRUE,null="null",pretty=TRUE,digits=NA)
cat(length(cases),"independent native factor cases generated\n")
