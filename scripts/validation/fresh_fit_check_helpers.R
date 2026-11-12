# Test-only exact fit/query encoders shared by matrix and formula checks.
bits<-function(v)lapply(as.numeric(v),function(x)paste0('bits:',hex(x)))
matrix_descriptor<-function(m)if(is.null(m))NULL else list(rows=nrow(m),columns=ncol(m),values=bits(m))
plain<-function(v)if(is.null(v))NULL else as.list(unname(v))
floats<-function(v)if(is.null(v))NULL else lapply(as.numeric(v),hex)
row_bits<-function(m)if(is.null(m))NULL else floats(as.vector(t(m)))
query_values<-function(f,queries)lapply(queries,function(q){
 cutoff<-switch(q$cooks,auto=TRUE,off=FALSE,as.numeric(q$cooks))
 r<-blocked(rsdeseq2::results(f,contrast=if(is.null(q$character))decode(q$contrast)else unlist(q$character),cooksCutoff=cutoff,
  alpha=decode(q$alpha),independentFiltering=q$filtering,lfcThreshold=decode(q$threshold),
  altHypothesis=q$alternative,useT=!is.null(q$df),degreesOfFreedom=if(is.null(q$df))NULL else decode(q$df)))
 v<-lapply(r[c('baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')],floats)
 v$filtered<-plain(r$filtered)
 m<-attr(r,'independentFiltering');scalar<-function(x)if(is.na(x))NULL else hex(x)
 v$filtering<-list(enabled=m$enabled,theta=floats(m$theta),numRej=plain(m$numRej),
  selected=if(is.na(m$selected))NULL else m$selected,filterTheta=scalar(m$filterTheta),filterThreshold=scalar(m$filterThreshold),
  lowessX=if(is.null(m$lo.fit))NULL else floats(m$lo.fit$x),lowessY=if(is.null(m$lo.fit))NULL else floats(m$lo.fit$y),alpha=hex(m$alpha))
 list(name=q$name,value=v)
})
encode<-function(f){
 result<-lapply(f$results,function(v)if(is.logical(v))plain(v)else floats(v))
 list(results=result,rawStat=floats(f$rawStat),rawPvalue=floats(f$rawPvalue),baseMean=floats(f$baseMean),baseVar=floats(f$baseVar),allZero=plain(f$allZero),
  sizeFactors=floats(f$sizeFactors),dispersions=floats(f$dispersions),dispFit=floats(f$dispFit),
  beta=row_bits(f$beta),betaSE=row_bits(f$betaSE),betaCovariance=floats(as.vector(aperm(f$betaCovariance,c(3,2,1)))),
  normalizationFactors=row_bits(f$normalizationFactors),design=row_bits(f$design),reduced=row_bits(f$reduced),
  mu=row_bits(f$mu),cooks=row_bits(f$cooks),lrtDf=f$lrtDf,nRefit=f$nRefit,replaceCounts=row_bits(f$replaceCounts),
  replace=plain(f$replace),replaceable=plain(f$replaceable),refitRows=plain(f$refitRows),newAllZeroRows=plain(f$newAllZeroRows))
}
canonical<-function(x){
 if(!is.list(x))return(x)
 if(!is.null(names(x)))x<-x[sort(names(x))]
 lapply(x,canonical)
}
