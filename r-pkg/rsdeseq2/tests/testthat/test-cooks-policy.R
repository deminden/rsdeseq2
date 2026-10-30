test_that("fitted row policy matches independent masks without R reductions or intermediate BH", {
  counts <- matrix(c(1,5,6,7,0, 1,5,6,0,0, 1,9,5,6,7,
    1,5,6,7,0, 9,5,6,7,0, 1,2,3,4,5),6,5,byrow=TRUE)
  cooks <- matrix(c(10,.1,.2,.3,.4, 10,.1,.2,.3,.4, 10,10,.2,.3,.4,
    NaN,NaN,NaN,NaN,NaN, 5,.1,.2,.3,.4, 10,.1,.2,.3,.4),6,5,byrow=TRUE)
  table <- data.frame(log2FoldChange=seq_len(6)*1.,lfcSE=rep(.5,6),
    stat=seq_len(6)*2.,pvalue=c(.01,.02,.03,.04,.05,.06),cooksOutlier=rep(FALSE,6))
  design <- cbind(1,c(0,0,1,1,1))
  groups <- factor(c("A","A","B","B","B"))
  base <- list(test="Wald",schemaVersion=1L,counts=counts,cooks=cooks,design=design,
    results=list(maxCooks=c(10,10,10,10,5,NA_real_)),
    refitRows=c(FALSE,TRUE,FALSE,FALSE,FALSE,FALSE),
    newAllZeroRows=c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE))
  variants <- list(matching=list(design=design,colData=data.frame(group=groups)),
    releveled=list(design=design,colData=data.frame(group=relevel(groups,"B"))),
    mismatched=list(design=design,colData=data.frame(group=factor(c("A","B","A","B","B")))),
    numeric=list(design=design,colData=data.frame(group=as.integer(groups))),
    additive=list(design=cbind(design,c(0,1,0,1,0)),colData=data.frame(group=groups)),
    missing=list(design=design,colData=data.frame(group=factor(c("A",NA,"B","B","B")))),
    third_level=list(design=design,colData=data.frame(group=factor(c("A","A","B","B","C")))))
  reference <- function(fit,cutoff) {
    out <- table; zero <- which(fit$newAllZeroRows)
    for(field in c("log2FoldChange","lfcSE","stat")) out[[field]][zero] <- 0
    out$pvalue[zero] <- 1
    out$cooksOutlier[] <- NA
    if(identical(cutoff,FALSE)) return(out)
    if(identical(cutoff,TRUE)) cutoff <- qf(.99,ncol(fit$design),nrow(fit$design)-ncol(fit$design))
    eligible <- ncol(fit$design)==2L && any(vapply(fit$colData,function(x)
      is.factor(x)&&nlevels(x)==2L&&identical(outer(unname(as.character(x)),unname(as.character(x)),"=="),
        outer(unname(fit$design[,2]),unname(fit$design[,2]),"==")),logical(1)))
    maximum <- fit$results$maxCooks; maximum[fit$refitRows|fit$newAllZeroRows] <- 0
    for(gene in seq_len(nrow(counts))) if(!is.na(maximum[gene])) {
      flag <- maximum[gene]>cutoff
      if(flag&&eligible) {
        finite <- which(is.finite(cooks[gene,]))
        index <- if(length(finite)) finite[which.max(cooks[gene,finite])] else integer()
        if(length(index)&&sum(counts[gene,]>counts[gene,index])>=3L) flag <- FALSE
      }
      out$cooksOutlier[gene] <- flag
      if(flag) out$pvalue[gene] <- NA_real_
    }
    out
  }
  policy <- getFromNamespace(".rsdeseq2_fit_cooks_policy","rsdeseq2")
  local_mocked_bindings(.rsdeseq2_bh_adjust=function(...)stop("unexpected preliminary BH"),.package="rsdeseq2")
  cases <- lapply(names(variants), function(name) list(variant=name,cutoff=5))
  cases <- c(cases,list(list(variant="matching",cutoff=0),list(variant="matching",cutoff=FALSE),
    list(variant="matching",cutoff=TRUE),list(variant="additive",cutoff=TRUE)))
  for(case in cases) {
    name <- case$variant;cutoff <- case$cutoff
    fit <- modifyList(base,variants[[name]])
    expected <- reference(fit,cutoff)
    # A production outer product must fail inside the guarded call; the
    # independent reference above deliberately retains the established rule.
    suppressMessages(trace("outer",where=baseenv(),print=FALSE,tracer=quote(stop("R outer used"))))
    actual <- tryCatch(policy(table,fit,cutoff),finally=suppressMessages(untrace("outer",where=baseenv())))
    expect_equal(actual,expected,tolerance=0,info=paste(name,cutoff))
  }
  fit <- modifyList(base,variants$matching)
  corrupt <- fit; corrupt$cooks[1,1] <- Inf
  expect_error(policy(table,corrupt,5),"cooks must contain finite")
  expect_equal(policy(table,fit,5)$cooksOutlier,c(FALSE,FALSE,FALSE,FALSE,FALSE,NA))
  cached <- fit; cached$rawPvalue <- table$pvalue/2
  expect_equal(policy(table,cached,FALSE,cached=TRUE)$pvalue,c(.005,.01,.015,1,.025,.03))
  legacy <- fit; legacy$rawPvalue <- NULL
  expect_error(policy(table,legacy,FALSE,cached=TRUE),"lacks raw p-values")
  legacy$test <- "LRT"; legacy$schemaVersion <- NULL
  expect_error(policy(table,legacy,5),"verified raw statistics")
})

test_that("fitted policy validates probabilities and maxima after replacement overrides", {
  policy <- getFromNamespace(".rsdeseq2_fit_cooks_policy","rsdeseq2")
  out <- data.frame(log2FoldChange=1.,lfcSE=.5,stat=2.,pvalue=.1)
  fit <- list(test="Wald",schemaVersion=1L,results=list(maxCooks=Inf),refitRows=FALSE,newAllZeroRows=FALSE)
  expect_error(policy(out,fit,5),"maxCooks")
  fit$refitRows <- TRUE
  expect_equal(policy(out,fit,5)$pvalue,.1)
  out$pvalue <- 1.2
  expect_error(policy(out,fit,5),"pvalue")
  fit$newAllZeroRows <- TRUE
  expect_equal(policy(out,fit,5)$pvalue,1)
  # Disabled filtering preserves a missing payload and never mutates input.
  out$pvalue <- NA_real_; before <- out
  fit$newAllZeroRows <- FALSE
  expect_identical(policy(out,fit,FALSE),before)
  expect_identical(out,before)
})

test_that("shared fitted policy matches frozen original DESeq2 formula and matrix decisions", {
  r <- jsonlite::fromJSON(test_path("data", "cooks_policy_r462.json"), simplifyVector=FALSE)
  values <- function(x) as.double(unlist(x))
  matrixRows <- function(x, ncol) matrix(values(x), ncol=ncol, byrow=TRUE)
  decode <- function(x) vapply(x,function(v) {
    if(is.null(v)) return(NA_real_)
    bytes <- as.raw(strtoi(substring(v,seq(1,15,2),seq(2,16,2)),16L))
    readBin(bytes,"double",n=1,size=8,endian="big")
  },numeric(1))
  bits <- function(x) lapply(x,function(v) if(is.na(v))NULL else writeBin(v,raw(),size=8,endian="big"))
  expect_length(r$cases,12L)
  for(case in r$cases) {
    fields <- lapply(case$factors,function(f) factor(as.integer(unlist(f$codes)),levels=seq_len(f$levels)))
    names(fields) <- vapply(case$factors,`[[`,character(1),"name")
    cd <- as.data.frame(fields)
    # Original numeric fields are absent from this factor-only context. They
    # must never acquire eligibility from an unrelated two-level factor.
    fit <- list(test="Wald",schemaVersion=2L,
      designOrigin=list(kind=case$kind,variables=as.character(unlist(case$variables))),
      colData=cd,design=matrixRows(case$design,case$coefficients),
      counts=matrixRows(r$counts,r$samples),cooks=matrixRows(r$cooks,r$samples),
      results=list(maxCooks=values(r$maximum)))
    table <- data.frame(log2FoldChange=rep(1,r$genes),lfcSE=rep(.5,r$genes),
      stat=rep(2,r$genes),pvalue=decode(case$raw_pvalue),cooksOutlier=rep(FALSE,r$genes))
    original <- table
    actual <- rsdeseq2:::.rsdeseq2_fit_cooks_policy(table,fit,as.double(r$cutoff))
    expect_identical(bits(actual$pvalue),bits(decode(case$pvalue)),info=case$name)
    expect_identical(bits(rsdeseq2:::.rsdeseq2_bh_adjust(actual$pvalue)),bits(decode(case$padj)),info=case$name)
    expect_identical(is.na(actual$pvalue)&!is.na(table$pvalue),as.logical(unlist(case$masked)),info=case$name)
    expect_identical(actual$cooksOutlier,as.logical(unlist(case$masked)),info=case$name)
    expect_identical(table,original,info=case$name)
  }
})
