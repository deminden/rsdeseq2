test_that("factor matrices and frames match independent original R references", {
  reference <- jsonlite::fromJSON(test_path("data","native_factor_r462.json"),simplifyVector=FALSE)
  decode <- function(hex)readBin(as.raw(strtoi(substring(hex,seq(1,15,2),seq(2,16,2)),16L)),"double",1L,size=8,endian="big")
  d <- as.data.frame(stats::setNames(lapply(reference$numeric,function(v)vapply(v$values,decode,0.0)),vapply(reference$numeric,`[[`,"","name")))
  for(v in reference$factors)d[[v$name]] <- factor(unlist(v$sample_levels),levels=unlist(v$levels))
  old <- options(na.action="na.pass");on.exit(options(old),add=TRUE)
  for(case in reference$native_formula_cases) {
    run <- function() {
      formula <- as.formula(case$formula)
      frame <- rsdeseq2:::.rsdeseq2_model_frame(formula,d)
      matrix <- rsdeseq2:::.rsdeseq2_model_matrix(formula,frame)
      if(any(!is.finite(matrix)))stop("non-finite fitting matrix")
      list(formula=formula,frame=frame,matrix=matrix)
    }
    if(is.null(case$expected)) { expect_error(run(),info=case$name);next }
    result <- run()
    bits <- function(matrix)lapply(unname(as.vector(t(matrix))),hex)
    expect_identical(bits(result$matrix[,unlist(case$r_columns),drop=FALSE]),case$expected,info=case$name)
    frames <- result$frame[vapply(result$frame,is.factor,TRUE)]
    expectedFrames <- lapply(case$frames,function(value){value$name<-NULL;value})
    actualFrames <- lapply(unname(frames),function(value)list(levels=as.list(levels(value)),reference=levels(value)[[1L]],sample_levels=as.list(as.character(value)),ordered=is.ordered(value)))
    expect_identical(actualFrames,expectedFrames,info=case$name)
    for(name in names(frames))if(!is.ordered(result$frame[[name]])) {
      identity <- diag(nlevels(result$frame[[name]]))
      dimnames(identity) <- list(levels(result$frame[[name]]),levels(result$frame[[name]]))
      attr(result$frame[[name]],"contrasts") <- identity
    }
    expanded <- rsdeseq2:::.rsdeseq2_model_matrix(result$formula,result$frame)
    expect_identical(bits(expanded[,unlist(case$r_expanded_columns),drop=FALSE]),case$expanded_expected,info=case$name)
  }
  for(case in reference$factor_expansions) {
    run <- function()rsdeseq2:::.rsdeseq2_formula_request(list(operation="expanded_factors",values=list(.formula_decode(case$input)),null_level=case$null_level,collation="C"))[[1L]]
    if(is.null(case$expected))expect_error(run(),info=case$name)
    else expect_identical(.formula_canonical(.formula_descriptor(run())),.formula_canonical(case$expected),info=case$name)
  }
})

test_that("ordered fitted state preserves polynomial selectors and rejects treatment selectors", {
  counts <- .rsdeseq2_test_counts()
  cd <- data.frame(g=factor(rep(c(".L",".Q"),each=4)),row.names=colnames(counts))
  formula <- ~ordered(g)
  dds <- DESeq2::DESeqDataSetFromMatrix(counts,cd,formula)
  object <- DESeq(dds,fitType="mean",sizeFactors=rep(1,ncol(counts)),minReplicatesForReplace=Inf)
  state <- S4Vectors::metadata(object)$rsdeseq2RustFit
  expect_identical(state$orderedFactors,"ordered(g)")
  expectedMatrix <- stats::model.matrix(formula,cd)
  retainedMatrix <- expectedMatrix
  attr(retainedMatrix,"assign") <- NULL
  attr(retainedMatrix,"contrasts") <- NULL
  expect_identical(state$design,retainedMatrix)
  explicit <- DESeq(list(counts=counts,design=expectedMatrix),fitType="mean",sizeFactors=rep(1,ncol(counts)),minReplicatesForReplace=Inf)
  expected <- results(explicit,contrast=c(0,1),cooksCutoff=FALSE,independentFiltering=FALSE)
  actual <- results(object,name="ordered.g..L",cooksCutoff=FALSE,independentFiltering=FALSE)
  for(column in c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))expect_identical(unname(actual[[column]]),unname(expected[[column]]),info=column)
  for(reference in list(NULL,".L")) {
    request <- list(contrast=c("ordered.g.",".Q",".L"),reference=reference)
    expect_error(do.call(results,c(list(object),request)),"ordered factor")
    expect_error(resultsList(object,list(request)),"ordered factor")
  }
  lrtSource <- state
  lrtSource$lrtStat <- rep(1,nrow(counts));lrtSource$lrtDf <- 1
  for(primitive in list(nbinomWaldTestRust(object=state),nbinomLRTRust(object=lrtSource))) {
    expect_identical(primitive$orderedFactors,state$orderedFactors)
    expect_error(results(primitive,contrast=c("ordered.g.",".Q",".L")),"ordered factor")
  }
  changed <- object
  metadata <- S4Vectors::metadata(changed);metadata$rsdeseq2RustFit$orderedFactors<-"g";S4Vectors::metadata(changed)<-metadata
  expect_error(results(changed),"ordered factor metadata changed")
  # Older valid fits recover class metadata from their already-validated frame.
  metadata <- S4Vectors::metadata(object);metadata$rsdeseq2RustFit$orderedFactors<-NULL;S4Vectors::metadata(object)<-metadata
  expect_error(results(object,contrast=c("ordered.g.",".Q",".L")),"ordered factor")
})
