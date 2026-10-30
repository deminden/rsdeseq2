test_that("bounded primitive batches preserve selectors, zero rules and complete single-result payloads", {
  input <- .rsdeseq2_primitive_fixture()
  input$counts <- matrix(c(0,0,0,0,5,6,0,1,1,0,5,6),2,6,byrow=TRUE)
  input$sampleLevels <- rep(c("A","B","C"),each=2)
  input$modelMatrix <- cbind(1,rep(c(0,1,0),each=2),rep(c(0,0,1),each=2))
  requests <- list(named=list(name="condition_B_vs_A",independentFiltering=FALSE),
    factor=list(contrast=c("condition","B","A"),independentFiltering=FALSE),
    numeric=list(contrast=c(0,1,0),alpha=.2),
    list=list(contrast=list("condition_C_vs_A","condition_B_vs_A"),reference="A"))
  for(method in c("Wald","LRT")) {
    fields <- input
    if(method=="LRT") fields <- c(fields,list(lrtStat=c(6,8),lrtPvalue=c(.02,.01)))
    fit <- do.call(if(method=="Wald")waldFitRust else lrtFitRust,fields)
    expected <- lapply(requests,function(args)do.call(results,c(list(object=fit),args)))
    actual <- resultsList(fit,requests,batchSize=2)
    expect_identical(actual,expected,info=method)
    expect_equal(actual$factor$log2FoldChange[1],0,info=method)
    expect_equal(actual$numeric$log2FoldChange[1],2,info=method)
    if(method=="LRT") expect_identical(actual$factor$pvalue,c(.02,.01))
  }
  # Explicit NULL in a request clears its common selector.
  fit <- do.call(waldFitRust,input)
  actual <- resultsList(fit,list(clear=list(name=NULL),keep=list()),name="condition_B_vs_A")
  expect_identical(actual$clear,results(fit))
  expect_identical(actual$keep,results(fit,name="condition_B_vs_A"))
  expect_identical(resultsList(fit,list()),list())
})

test_that("dataset validation and primitive extraction happen once across bounded chunks", {
  ct <- .rsdeseq2_test_counts();ct[1,8] <- 2000L
  cd <- data.frame(condition=factor(rep(c("A","B"),each=4)),row.names=colnames(ct))
  fit <- DESeq(DESeq2::DESeqDataSetFromMatrix(ct,cd,~condition),
    fitType="mean",minReplicatesForReplace=Inf)
  coefficient <- tail(resultsNames(fit),1)
  requests <- list(default=list(),named=list(name=coefficient),numeric=list(contrast=c(0,1)),
    unmasked=list(name=coefficient,cooksCutoff=FALSE),threshold=list(name=coefficient,lfcThreshold=.75,altHypothesis="greater"))
  expected <- lapply(requests,function(args)do.call(results,c(list(object=fit,independentFiltering=FALSE),args)))
  original_view <- getFromNamespace(".rsdeseq2_fitted_view","rsdeseq2")
  original_primitive <- nbinomWaldTestRust
  original_batch <- getFromNamespace(".rsdeseq2_evaluate_result_batch","rsdeseq2")
  views <- primitives <- batches <- 0L
  local_mocked_bindings(
    .rsdeseq2_fitted_view=function(...) { views <<- views+1L;original_view(...) },
    nbinomWaldTestRust=function(...) { primitives <<- primitives+1L;original_primitive(...) },
    .rsdeseq2_evaluate_result_batch=function(...) { batches <<- batches+1L;original_batch(...) },
    .package="rsdeseq2")
  actual <- resultsList(fit,requests,independentFiltering=FALSE,batchSize=2L)
  expect_identical(actual,expected)
  expect_identical(c(views,primitives,batches),c(1L,1L,2L))
  expect_true(is.na(actual$named$pvalue[1]))
  expect_true(is.finite(actual$unmasked$pvalue[1]))
  subset <- fit[c(5,3,1),]
  expect_identical(resultsList(subset,list(named=list(name=coefficient)),independentFiltering=FALSE)[[1]],
    results(subset,name=coefficient,independentFiltering=FALSE))
  changed <- fit;SummarizedExperiment::assay(changed,"counts")[1,1] <- ct[1,1]+1L
  before <- batches
  expect_error(resultsList(changed,requests),"counts changed")
  expect_identical(batches,before)
})

test_that("batch requests independently select filtering and validate request options", {
  set.seed(20261029);n <- 2002L
  filter <- c(rep(0,102),seq_len(n-102L))
  p <- runif(n)^ifelse(seq_len(n)>800,5,1)
  fit <- waldFitRust(filter,cbind(Intercept=1,effect=stats::qnorm(p/2,lower.tail=FALSE)*.5),
    betaSE=matrix(.5,n,2))
  requests <- list(tight=list(name="effect",alpha=.001),loose=list(name="effect",alpha=.2),
    reverse=list(contrast=c(0,-1),altHypothesis="greater",lfcThreshold=.75),
    disabled=list(name="effect",independentFiltering=FALSE),
    t=list(name="effect",useT=TRUE,degreesOfFreedom=10))
  expected <- lapply(requests,function(args)do.call(results,c(list(object=fit),args)))
  actual <- resultsList(fit,requests,batchSize=3)
  expect_identical(actual,expected)
  expect_false(identical(attr(actual$tight,"independentFiltering")$numRej,
    attr(actual$loose,"independentFiltering")$numRej))
  expect_false(identical(actual$reverse$pvalue,actual$loose$pvalue))
  expect_error(resultsList(fit,list(list(unknownOption=TRUE))),"unknownOption")
  expect_error(resultsList(fit,list(list(alpha=1))),"alpha")
  expect_error(resultsList(fit,list(list(cooksCutoff=TRUE))),"retained Cook")
  expect_error(resultsList(fit,list(list()),batchSize=65),"batchSize")
  expect_error(resultsList(fit,list(list(name="missing"))),"missing")
  expect_error(resultsList(fit,list(c("effect"))),"named.*list")
})

test_that("replacement state and LRT inference survive batched and cached requests", {
  counts <- .rsdeseq2_test_counts(); counts[1L,8L] <- 2000L
  counts[2L,] <- 0L; counts[2L,8L] <- 1000L
  coldata <- data.frame(condition=factor(rep(c("A","B"),each=4)),row.names=colnames(counts))
  for(method in c("Wald","LRT")) {
    dds <- DESeq2::DESeqDataSetFromMatrix(counts,coldata,~condition)
    fit <- DESeq(dds,test=method,reduced=if(method=="LRT")~1 else NULL,
      fitType="mean",minReplicatesForReplace=3L)
    state <- S4Vectors::metadata(fit)$rsdeseq2RustFit
    expect_equal(state$nRefit,2L,info=method)
    coefficient <- tail(resultsNames(fit),1)
    requests <- list(default=list(),named=list(name=coefficient),
      numeric=list(contrast=c(0,1)),factor=list(contrast=c("condition","B","A")),
      unmasked=list(name=coefficient,cooksCutoff=FALSE),numeric_cutoff=list(name=coefficient,cooksCutoff=.1))
    expected <- lapply(requests,function(args)do.call(results,c(list(object=fit),args)))
    actual <- resultsList(fit,requests,batchSize=2)
    expect_identical(actual,expected,info=method)
    expect_equal(actual$numeric$log2FoldChange[2],0,info=method)
    expect_equal(actual$numeric$pvalue[2],1,info=method)
    expect_equal(S4Vectors::metadata(actual$default)$cooksReplacement$nRefit,2L,info=method)
    # Cache-only queries keep supporting legacy states lacking coefficient arrays.
    legacy <- state; legacy$beta <- legacy$betaSE <- legacy$betaCovariance <- NULL
    expect_identical(resultsList(legacy,list(list(),list(cooksCutoff=FALSE)))[[2]],
      results(legacy,cooksCutoff=FALSE),info=method)
  }
})
