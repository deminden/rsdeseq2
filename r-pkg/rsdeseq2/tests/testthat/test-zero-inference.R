testthat::test_that("shared zero policy matches independent count/design rules for Wald and LRT display", {
  counts <- matrix(c(0,0,0,0,5,6, 0,0,0,0,0,0, 0,1,0,0,5,6, 0,0,1,0,0,0),4,6,byrow=TRUE)
  design <- diag(3)[c(1,1,2,2,3,3),]
  beta <- matrix(rep(c(1,2,4),4),4,3,byrow=TRUE)
  se <- matrix(rep(c(.5,.75,1),4),4,3,byrow=TRUE)
  covariance <- array(0,c(4,3,3))
  for(gene in 1:4) covariance[gene,,] <- diag(c(.25,.5625,1))
  cases <- list(
    factor=list(type="character",contrast=c(-1,1,0)),
    numeric=list(type="numeric",contrast=c(-1,1,0)),
    one_sided=list(type="numeric",contrast=c(0,1,0),noDesign=TRUE),
    cancelled=list(type="numeric",contrast=c(-1,1,0),cancel=TRUE)
  )
  for(name in names(cases)) {
    case <- cases[[name]]
    model <- if(isTRUE(case$noDesign))NULL else design
    if(isTRUE(case$cancel)) model <- cbind(rep(-1,6),rep(1,6),rep(0,6))
    fit <- list(beta=beta,betaSE=se,betaCovariance=covariance,counts=counts,
      modelMatrix=model,sampleLevels=rep(c("A","B","C"),each=2))
    rule <- list(type=case$type,numerator="B",denominator="A")
    # Independent R reference calculations belong only in tests.
    flags <- if(case$type=="character") rowSums(counts[,1:4,drop=FALSE])==0 else
      if(all(case$contrast>=0)||all(case$contrast<=0))rep(FALSE,4) else
      rowSums(counts[,as.numeric(model %*% as.numeric(case$contrast!=0))!=0,drop=FALSE])==0
    # Exclude the globally zero row even though its supplied estimates are finite.
    eligible <- flags & rowSums(counts)!=0
    effect <- as.numeric(beta %*% case$contrast)
    standardError <- sqrt(as.numeric(t(case$contrast) %*% covariance[1,,] %*% case$contrast))
    raw <- rsdeseq2:::.rsdeseq2_lrt_from_primitive_fit(fit,case$contrast,rule)
    expectedEffect <- effect; expectedEffect[eligible] <- 0
    testthat::expect_equal(raw$log2FoldChange,expectedEffect,tolerance=1e-14,info=name)
    testthat::expect_equal(raw$lfcSE,rep(standardError,4),tolerance=1e-14,info=name)
    df <- c(NA_real_,NA_real_,2,100)
    options <- list(lfcThreshold=.75,altHypothesis="greater",useT=TRUE,degreesOfFreedom=df)
    actual <- rsdeseq2:::.rsdeseq2_wald_from_primitive_fit(fit,case$contrast,options,rule)
    statistic <- pmax((effect-.75)/standardError,0)
    pvalue <- stats::pt((effect-.75)/standardError,df,lower.tail=FALSE)
    statistic[eligible] <- 0; pvalue[eligible] <- 1
    testthat::expect_equal(actual$log2FoldChange,expectedEffect,tolerance=1e-14,info=name)
    testthat::expect_equal(actual$lfcSE,rep(standardError,4),tolerance=1e-14,info=name)
    testthat::expect_equal(actual$stat,statistic,tolerance=1e-14,info=name)
    testthat::expect_equal(actual$pvalue,pvalue,tolerance=1e-13,info=name)
  }
})
