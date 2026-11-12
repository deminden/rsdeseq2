#!/usr/bin/env Rscript
# R is a test/reference process here, never a production numerical backend.
args <- commandArgs(TRUE)
if (length(args) != 3L) stop("usage: STANDALONE_BINARY R_PACKAGE_LIBRARY OUTPUT_DIR")
binary <- normalizePath(args[[1L]])
library(rsdeseq2, lib.loc = normalizePath(args[[2L]]))
symbols <- system2("nm", c("-u", shQuote(getLoadedDLLs()[["rsdeseq2"]][["path"]])), stdout = TRUE, stderr = TRUE)
if (!is.null(attr(symbols, "status"))) stop("Cannot inspect native bridge numerical dependencies")
if (any(grepl("(^|[ _])(Rf_)?(pnorm5|pt|pchisq|qchisq|qf|dnorm4|lgammafn|digamma)$", trimws(symbols))))
  stop("Native bridge imports an R numerical backend")
output <- args[[3L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
hex <- function(value) {
  if (is.na(value)) return(NULL)
  paste(format(writeBin(as.double(value), raw(), size = 8, endian = "big")), collapse = "")
}
hexVector <- function(values) lapply(unname(values), hex)
withoutRCalculations <- function(expression, projection = FALSE, zeros = FALSE, model = FALSE) {
  # Fail at runtime if the production wrapper reintroduces statistical R calls.
  functions <- list(stats = c("p.adjust", "quantile", "lowess", "pnorm", "pt", "pchisq", "qf", "dnorm"), base = c("mean","exp"))
  if (projection) functions$base <- c(functions$base, "%*%", "sqrt")
  if (zeros) functions$base <- c(functions$base, "rowSums", "colSums", "sum")
  if (model) functions$stats <- c(functions$stats,"model.matrix","model.matrix.default","contr.treatment","contr.sum","contr.helmert","contr.poly","contr.SAS")
  on.exit(for (package in names(functions)) for (name in functions[[package]])
    suppressMessages(untrace(name, where = asNamespace(package))), add = TRUE)
  for (package in names(functions)) for (name in functions[[package]])
    suppressMessages(trace(name, where = asNamespace(package), print = FALSE,
      tracer = quote(stop("Production wrapper executed an R numerical function", call. = FALSE))))
  force(expression)
}
guard <- tryCatch(withoutRCalculations(stats::p.adjust(c(.1, .2))), error = identity)
stopifnot(inherits(guard, "error"), grepl("Production wrapper executed", conditionMessage(guard)))
guard <- tryCatch(withoutRCalculations(base::exp(0)),error=identity)
stopifnot(inherits(guard,"error"),grepl("Production wrapper executed",conditionMessage(guard)))
guard <- tryCatch(withoutRCalculations(stats::model.matrix(~1,data.frame(x=1:2)),model=TRUE),error=identity)
stopifnot(inherits(guard,"error"),grepl("Production wrapper executed",conditionMessage(guard)))
set.seed(20261029)
n <- 2002L
filter <- c(rep(0, 102), seq_len(n - 102L))
p <- runif(n)^ifelse(seq_len(n) > 800, 5, 1)
p[c(3, 120, 1750)] <- NA_real_
cases <- list(
  default = list(filter = filter, p = p),
  custom = list(filter = filter, p = p, theta = c(0, .15, .3, .5, .75, .95)),
  repeated = list(filter = filter, p = p, theta = c(1, 1)),
  unsorted = list(filter = filter, p = p, theta = c(.9, .15, .5, .3)),
  all_zero = list(filter = rep(0, n), p = p),
  all_missing = list(filter = filter, p = rep(NA_real_, n)),
  disabled = list(filter = filter, p = p, enabled = FALSE),
  small_alpha = list(filter = filter, p = p, alpha = .01),
  tails_ties = list(filter = seq_len(20L), p = c(NA_real_, 0, 1e-320, 1e-200, rep(.01, 16L)))
)
compareIdentity <- function(name, input, actual) {
  path <- file.path(output, paste0(name, "-input.json"))
  jsonlite::write_json(input, path, auto_unbox = TRUE, null = "null")
  standalone <- system2(binary, shQuote(path), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(standalone, "status"))) stop(paste(standalone, collapse = "\n"))
  expected <- jsonlite::fromJSON(paste(standalone, collapse = "\n"), simplifyVector = FALSE)
  # JSON normalizes list ordering and scalar integer transport, never f64 values:
  # all floating-point outputs above are exact IEEE-754 hex, not rounded decimals.
  actual <- jsonlite::fromJSON(jsonlite::toJSON(actual, auto_unbox = TRUE, null = "null"), simplifyVector = FALSE)
  actual <- actual[sort(names(actual))]; expected <- expected[sort(names(expected))]
  writeLines(paste(standalone, collapse = "\n"), file.path(output, paste0(name, "-standalone.json")))
  jsonlite::write_json(actual, file.path(output, paste0(name, "-wrapped.json")), auto_unbox = TRUE, null = "null", pretty = TRUE)
  if (!identical(actual, expected)) stop(sprintf("Cross-interface bit mismatch: %s\n%s", name,
    paste(all.equal(actual, expected, tolerance = 0), collapse = "\n")))
  cat(name, ": every output and metadata value identical\n", sep = "")
}

for (name in names(cases)) {
  case <- cases[[name]]
  enabled <- if (is.null(case$enabled)) TRUE else case$enabled
  alpha <- if (is.null(case$alpha)) .1 else case$alpha
  input <- list(filter = as.list(sprintf("%.17g", case$filter)),
    pvalues = lapply(case$p, function(v) if (is.na(v)) NULL else sprintf("%.17g", v)),
    theta = if (is.null(case$theta)) NULL else as.list(sprintf("%.17g", case$theta)),
    alpha = alpha |> sprintf(fmt = "%.17g"), enabled = enabled)
  compareIdentity(paste0("filter_", name), input, {
  result <- withoutRCalculations({
    table <- rsdeseq2::resultsTableRust(case$filter, pvalue = case$p)
    rsdeseq2::applyIndependentFilteringRust(table, alpha = alpha, theta = case$theta, enabled = enabled)
  })
  meta <- attr(result, "independentFiltering")
  actual <- list(padj = hexVector(result$padj),
    filtered = lapply(result$filtered, function(v) if (is.na(v)) NULL else v),
    enabled = meta$enabled, theta = hexVector(meta$theta), numRej = as.list(meta$numRej),
    selected = if (is.na(meta$selected)) NULL else meta$selected,
    filterTheta = hex(meta$filterTheta), filterThreshold = hex(meta$filterThreshold),
    lowessX = if (is.null(meta$lo.fit)) NULL else hexVector(meta$lo.fit$x),
    lowessY = if (is.null(meta$lo.fit)) NULL else hexVector(meta$lo.fit$y), alpha = hex(meta$alpha))
    actual
  })
}

baseCounts <- matrix(c(10L,40L,90L,30L,40L,50L,20L,200L,100L,0L,0L,0L),4L,byrow=TRUE)
zeroCounts <- matrix(c(0L,4L,16L,2L,9L,0L,36L,4L,1L,2L,3L,4L,0L,0L,0L,0L,7L,11L,13L,17L),5L,byrow=TRUE)
normalizationCases <- list(
  ratio = list(counts = baseCounts),
  poscounts = list(counts = zeroCounts, method = "poscounts"),
  geometric = list(counts = baseCounts, geo = c(5,25,50,0)),
  geometric_zero = list(counts = baseCounts, geo = c(0,25,50,0)),
  geometric_infinite = list(counts = baseCounts, geo = c(5,Inf,50,0)),
  geometric_large = list(counts = baseCounts, geo = c(1000,2000,3000,0)),
  controls = list(counts = zeroCounts, control = c(1L,3L,5L)),
  named_controls = list(counts = zeroCounts, control = c(1L,3L,5L), named = TRUE),
  supplied_factors = list(counts = baseCounts, sizeFactors = c(.7,1.1,1.8)),
  normalization_factors = list(counts = baseCounts, sizeFactors = c(.7,1.1,1.8),
    nf = matrix(seq_len(12L)/7+.5,4L)),
  weights = list(counts = baseCounts, nf = matrix(seq_len(12L)/7+.5,4L),
    weights = matrix(rep(c(0,.3,.8,1),3L),4L))
)
optionalNumbers <- function(values) if (is.null(values)) NULL else as.list(sprintf("%.17g",values))
for (name in names(normalizationCases)) {
  case <- normalizationCases[[name]]; counts <- case$counts
  rownames(counts) <- paste0("gene",seq_len(nrow(counts)))
  colnames(counts) <- paste0("sample",seq_len(ncol(counts)))
  method <- if (is.null(case$method)) "ratio" else case$method
  control <- if (isTRUE(case$named)) rownames(counts)[case$control] else case$control
  input <- list(genes = nrow(counts), samples = ncol(counts),
    counts = as.list(as.integer(t(counts))), method = method,
    geo_means = optionalNumbers(case$geo),
    control_genes = if (is.null(case$control)) NULL else as.list(case$control - 1L),
    size_factors = optionalNumbers(case$sizeFactors),
    normalization_factors = optionalNumbers(if (is.null(case$nf)) NULL else as.vector(t(case$nf))),
    weights = optionalNumbers(if (is.null(case$weights)) NULL else as.vector(t(case$weights))))
  compareIdentity(paste0("normalization_",name),input,{
    primitive <- withoutRCalculations({
      estimated <- rsdeseq2::estimateSizeFactorsRust(counts,method=method,geoMeans=case$geo,controlGenes=control)
      factors <- if (is.null(case$sizeFactors)) estimated else case$sizeFactors
      normalized <- rsdeseq2::normalizedCountsRust(counts,sizeFactors=factors,normalizationFactors=case$nf)
      metadata <- rsdeseq2::baseMetadataRust(counts,sizeFactors=factors,normalizationFactors=case$nf,weights=case$weights)
      list(estimated=estimated,normalized=normalized,metadata=metadata)
    })
    list(sizeFactors=hexVector(primitive$estimated),normalized=hexVector(as.vector(t(primitive$normalized))),
      baseMean=hexVector(primitive$metadata$baseMean),baseVar=hexVector(primitive$metadata$baseVar),
      allZero=as.list(primitive$metadata$allZero))
  })
}
capture.output(sessionInfo(), file = file.path(output, "R-session.txt"))

# Inference primitives exercise every alternative, invalid/missing rows, scalar
# and per-gene t df, LRT tails and both present/absent Cook's cutoff metadata.
values <- c(-4,-.25,0,.25,2,NaN,Inf,1e308)
errors <- c(.25,.5,1,.5,1,1,1,1e-308)
strings <- function(values) as.list(sprintf("%.17g",values))
alternatives <- c("greaterAbs","greaterAbsUPSHOT","greaterAbs2014","lessAbs","greater","less")
for (t in c(FALSE,TRUE)) for (alternative in alternatives) {
  if (t && alternative=="greaterAbsUPSHOT") next
  df <- if (t) c(.5,1,5,10,25,NA_real_,-1,Inf) else NULL
  # Represent missing df as IEEE NaN for both transports, independently of JSON null.
  if (t) df[is.na(df)] <- NaN
  input <- list(operation="wald",beta=strings(values),se=strings(errors),threshold="0.75",
                alternative=alternative,df=if(t) strings(df) else NULL)
  actual <- withoutRCalculations(rsdeseq2:::.rsdeseq2_wald_stat_pvalue_with_options(values,errors,
    list(lfcThreshold=.75,altHypothesis=alternative,useT=t,degreesOfFreedom=df)))
  compareIdentity(paste0("wald_",alternative,"_",t),input,
    list(stat=hexVector(actual$stat),pvalue=hexVector(actual$pvalue)))
}
for (df in c(.5,2,10,100,1e6,1e12)) {
  input <- list(operation="wald",beta=strings(values),se=strings(errors),threshold="0",
                alternative="greaterAbs",df=strings(df))
  actual <- withoutRCalculations(.Call("rsdeseq2_wald_statistics",values,errors,0,0L,TRUE,df,PACKAGE="rsdeseq2"))
  compareIdentity(paste0("wald_scalar_df_",df),input,list(stat=hexVector(actual$stat),pvalue=hexVector(actual$pvalue)))
}
for (df in list(2,c(.5,1,2,3,10,100,1e6,NaN))) {
  stat <- c(0,.1,2,6,40,100,1e6,NaN)
  input <- list(operation="lrt",statistics=strings(stat),df=strings(df))
  actual <- withoutRCalculations(.Call("rsdeseq2_lrt_pvalues",stat,df,PACKAGE="rsdeseq2"))
  compareIdentity(paste0("lrt_",length(df)),input,list(pvalue=hexVector(actual)))
}
for (pair in list(c(8L,2L),c(4L,4L),c(61L,3L),c(52L,2L),c(1000L,500L),c(1000000L,500000L))) {
  input <- list(operation="cooks",samples=pair[1],coefficients=pair[2])
  actual <- withoutRCalculations(.Call("rsdeseq2_default_cooks_cutoff",pair[1],pair[2],PACKAGE="rsdeseq2"))
  compareIdentity(paste0("cooks_",pair[1],"_",pair[2]),input,list(cutoff=if(is.null(actual)) NULL else hex(actual)))
}

emptyInput <- list(operation="wald",beta=list(),se=list(),threshold="0",alternative="greaterAbs",df=list())
emptyActual <- withoutRCalculations(.Call("rsdeseq2_wald_statistics",numeric(),numeric(),0,0L,TRUE,numeric(),PACKAGE="rsdeseq2"))
compareIdentity("wald_empty_t",emptyInput,list(stat=hexVector(emptyActual$stat),pvalue=hexVector(emptyActual$pvalue)))

# The standalone input is row-major; the wrapper borrows R's column-major
# beta/covariance storage. Logical operation order must agree exactly.
projectionBeta <- matrix(c(1,2,4,-2,1,3,5,-1,2,7,2,-3),4,3,byrow=TRUE)
projectionSE <- matrix(rep(c(.5,.75,1),4),4,3,byrow=TRUE)
projectionCov <- array(0,c(4,3,3))
for(gene in 1:4) projectionCov[gene,,] <- matrix(c(1,.125,-.25,.125,2,.375,-.25,.375,3),3,3)*gene
projectionCases <- list(
  covariance=list(contrast=c(0,-1,1)),
  fractional=list(contrast=c(.25,-.5,2)),
  marginal=list(contrast=c(0,-2,0)),
  covariance_only=list(contrast=c(0,-2,0),noSE=TRUE),
  se_only=list(contrast=c(0,-2,0),noCov=TRUE),
  invalid_marginal=list(contrast=c(0,-2,0),invalidSE=TRUE)
)
for(name in names(projectionCases)) {
  case <- projectionCases[[name]]
  se <- if(isTRUE(case$noSE)) NULL else projectionSE
  covariance <- if(isTRUE(case$noCov)) NULL else projectionCov
  if(isTRUE(case$invalidSE)) se[,2] <- c(.75,0,NaN,2)
  packedCov <- if(is.null(covariance)) NULL else unlist(lapply(1:4,function(gene)as.vector(t(covariance[gene,,]))))
  for(mode in c("raw","normal","t")) {
    useT <- identical(mode,"t")
    df <- if(useT) c(.5,10,NaN,100) else NULL
    options <- list(lfcThreshold=.75,altHypothesis="greater",useT=useT,degreesOfFreedom=df)
    fit <- list(beta=projectionBeta,betaSE=se,betaCovariance=covariance)
    input <- list(operation="projection",genes=4L,coefficients=3L,
      beta=strings(as.vector(t(projectionBeta))),se=if(is.null(se))NULL else strings(as.vector(t(se))),
      covariance=if(is.null(packedCov))NULL else strings(packedCov),contrast=strings(case$contrast),
      threshold=if(mode=="raw")NULL else ".75",alternative="greater",df=if(useT)strings(df)else NULL)
    actual <- withoutRCalculations(if(mode=="raw")
      rsdeseq2:::.rsdeseq2_lrt_from_primitive_fit(fit,case$contrast) else
      rsdeseq2:::.rsdeseq2_wald_from_primitive_fit(fit,case$contrast,options), projection=TRUE)
    compareIdentity(paste0("projection_",name,"_",mode),input,lapply(actual,hexVector))
  }
}

# Zero rules are evaluated within the same Rust projection call. R only aligns
# metadata. Include a globally zero gene with supplied finite estimates to
# distinguish global exclusion from the contrast-specific override.
zeroCounts <- matrix(c(0,0,0,0,5,6, 0,0,0,0,0,0, 0,1,0,0,5,6, 0,0,1,0,0,0),4,6,byrow=TRUE)
zeroDesign <- diag(3)[c(1,1,2,2,3,3),]
zeroCases <- list(
  factor=list(rule="character",contrast=c(-1,1,0)),
  factor_reference=list(rule="character",contrast=c(0,1,0),noDesign=TRUE),
  numeric=list(rule="numeric",contrast=c(-1,1,0)),
  fractional=list(rule="numeric",contrast=c(-2,.5,0)),
  one_sided=list(rule="numeric",contrast=c(0,1,0),noDesign=TRUE),
  cancelled_design=list(rule="numeric",contrast=c(-1,1,0),cancel=TRUE),
  unicode=list(rule="character",contrast=c(-1,1,0),unicode=TRUE)
)
for(name in names(zeroCases)) for(mode in c("raw","normal","t")) {
  case <- zeroCases[[name]]
  design <- if(isTRUE(case$noDesign))NULL else zeroDesign
  if(isTRUE(case$cancel)) design <- cbind(rep(-1,6),rep(1,6),rep(0,6))
  levels <- if(isTRUE(case$unicode))rep(c("A β","B β","C β"),each=2)else rep(c("A","B","C"),each=2)
  numerator <- levels[3]; denominator <- levels[1]
  fit <- list(beta=projectionBeta,betaSE=projectionSE,betaCovariance=projectionCov,
    counts=zeroCounts,modelMatrix=design,sampleLevels=levels)
  allZero <- list(type=case$rule,numerator=numerator,denominator=denominator)
  df <- if(mode=="t")c(NaN,NaN,2,100)else NULL
  options <- list(lfcThreshold=.75,altHypothesis="greater",useT=mode=="t",degreesOfFreedom=df)
  input <- list(operation="projection",genes=4L,coefficients=3L,
    beta=strings(as.vector(t(projectionBeta))),se=strings(as.vector(t(projectionSE))),
    covariance=strings(unlist(lapply(1:4,function(gene)as.vector(t(projectionCov[gene,,]))))),
    contrast=strings(case$contrast),threshold=if(mode=="raw")NULL else ".75",alternative="greater",
    df=if(is.null(df))NULL else strings(df),
    zero=list(samples=6L,counts=as.list(as.integer(t(zeroCounts))),
      design=if(is.null(design))NULL else strings(as.vector(t(design))),levels=as.list(levels),
      numerator=numerator,denominator=denominator,rule=case$rule))
  actual <- withoutRCalculations(if(mode=="raw")
    rsdeseq2:::.rsdeseq2_lrt_from_primitive_fit(fit,case$contrast,allZero) else
    rsdeseq2:::.rsdeseq2_wald_from_primitive_fit(fit,case$contrast,options,allZero),
    projection=TRUE,zeros=TRUE)
  compareIdentity(paste0("zero_",name,"_",mode),input,lapply(actual,hexVector))
}

# Cross-interface prepared batches, including mixed raw/Normal/t requests and
# both layouts. A 65-request public call additionally crosses a native chunk.
batchRequests <- lapply(seq_len(65L),function(i) {
  mode <- (i-1L) %% 3L
  list(contrast=c(-1,1,0)*if(i %% 2L)1 else -2,
    threshold=if(mode==0L)NULL else .75,alternative="greater",
    df=if(mode==2L)c(.5,10,NaN,100)else NULL,
    zero=c("character","numeric","none")[[mode+1L]],numerator="B",denominator="A")
})
batchInput <- list(genes=4L,coefficients=3L,beta=strings(as.vector(t(projectionBeta))),
  se=strings(as.vector(t(projectionSE))),
  covariance=strings(unlist(lapply(1:4,function(g)as.vector(t(projectionCov[g,,]))))),
  samples=6L,counts=as.list(as.integer(t(zeroCounts))),design=strings(as.vector(t(zeroDesign))),
  levels=as.list(rep(c("A","B","C"),each=2)),
  requests=lapply(batchRequests,function(r) {
    r$contrast <- strings(r$contrast);r["threshold"] <- list(if(is.null(r$threshold))NULL else sprintf("%.17g",r$threshold))
    r["df"] <- list(if(is.null(r$df))NULL else strings(r$df));r
  }))
transport <- lapply(batchRequests,function(r)list(contrast=r$contrast,
  wald=if(is.null(r$threshold))NULL else list(threshold=r$threshold,alternative=4L,useT=!is.null(r$df),df=r$df),
  zero=if(r$zero=="none")NULL else if(r$zero=="numeric")list(rule=1L,counts=zeroCounts,design=zeroDesign)else
    list(rule=2L,counts=zeroCounts,sampleLevels=rep(c("A","B","C"),each=2),numerator="B",denominator="A")))
actual <- withoutRCalculations(c(
  .Call("rsdeseq2_contrast_batch",projectionBeta,projectionSE,projectionCov,transport[1:64],PACKAGE="rsdeseq2"),
  .Call("rsdeseq2_contrast_batch",projectionBeta,projectionSE,projectionCov,transport[65],PACKAGE="rsdeseq2")),projection=TRUE,zeros=TRUE)
compareIdentity("prepared_batch_65",batchInput,list(results=lapply(actual,function(r)lapply(r[sort(names(r))],hexVector))))

publicFit <- waldFitRust(rep(20,4),projectionBeta,betaSE=projectionSE,betaCovariance=projectionCov,
  counts=zeroCounts,modelMatrix=zeroDesign,sampleLevels=rep(c("A","B","C"),each=2),
  resultsNames=c("Intercept","condition_B_vs_A","condition_C_vs_A"))
publicRequests <- lapply(seq_len(65L),function(i)list(contrast=c(-1,1,0)*if(i%%2L)1 else -2,
  lfcThreshold=.75,altHypothesis="greater",independentFiltering=FALSE,
  useT=i%%3L==0L,degreesOfFreedom=if(i%%3L==0L)c(.5,10,NaN,100)else NULL))
publicInput <- batchInput
publicInput$requests <- lapply(publicRequests,function(r)list(contrast=strings(r$contrast),threshold=".75",
  alternative="greater",df=if(r$useT)strings(r$degreesOfFreedom)else NULL,
  zero="numeric",numerator="B",denominator="A"))
publicActual <- withoutRCalculations(resultsList(publicFit,publicRequests,batchSize=32L),projection=TRUE,zeros=TRUE)
compareIdentity("public_batch_65",publicInput,list(results=lapply(publicActual,function(r)
  lapply(as.list(r[c("lfcSE","log2FoldChange","pvalue","stat")]),hexVector))))

# Borrowed Cook/refit row policy: compare mask-only and fitted presentation with
# row-major unsigned standalone inputs and column-major double wrapper inputs.
cpCounts <- matrix(c(1,5,6,7,0, 1,5,6,0,0, 1,9,5,6,7, 1,5,6,7,0, 9,5,6,7,0, 1,2,3,4,5),6,5,byrow=TRUE)
cpCooks <- matrix(c(10,.1,.2,.3,.4,10,.1,.2,.3,.4,10,10,.2,.3,.4,
  NaN,NaN,NaN,NaN,NaN,5,.1,.2,.3,.4,10,.1,.2,.3,.4),6,5,byrow=TRUE)
cpP <- c(.01,.02,.03,.04,.05,.06);cpMax <- c(10,10,10,10,5,NA_real_)
cpDesign <- cbind(1,c(0,0,1,1,1))
cpTable <- data.frame(log2FoldChange=as.double(1:6),lfcSE=rep(.5,6),stat=as.double(2*(1:6)),pvalue=cpP,cooksOutlier=rep(FALSE,6))
cpBase <- list(test="Wald",schemaVersion=1L,results=list(maxCooks=cpMax),counts=cpCounts,cooks=cpCooks,design=cpDesign,
  refitRows=c(FALSE,TRUE,FALSE,FALSE,FALSE,FALSE),newAllZeroRows=c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE))
cpInput <- function(fit,cutoff,mode,low_count=NULL) list(cooks_policy=list(mode=mode,samples=5L,
  pvalue=lapply(cpP,function(v)if(is.na(v))NULL else sprintf("%.17g",v)),maximum=lapply(cpMax,function(v)if(is.na(v))NULL else sprintf("%.17g",v)),
  effect=strings(cpTable$log2FoldChange),se=strings(cpTable$lfcSE),stat=strings(cpTable$stat),cutoff=if(is.null(cutoff))NULL else sprintf("%.17g",cutoff),
  counts=as.list(as.integer(t(cpCounts))),cooks=strings(as.vector(t(cpCooks))),
  design=strings(as.vector(t(fit$design))),coefficients=ncol(fit$design),
  factors=unname(lapply(Filter(is.factor,fit$colData),function(x)list(codes=as.list(replace(as.integer(x),is.na(x),-2147483648)),levels=nlevels(x)))),
  low_count=low_count,refitted=as.list(fit$refitRows),new_zero=as.list(fit$newAllZeroRows)))
cpGroups <- factor(c("A","A","B","B","B"))
for(low_count in c(FALSE,TRUE)) for(cutoff in c(0,5,10)) {
  fit <- cpBase;fit$refitRows[] <- FALSE;fit$newAllZeroRows[] <- FALSE
  fit$colData <- data.frame(group=cpGroups)
  actual <- withoutRCalculations(applyCooksCutoffRust(cpP,cpMax,cutoff,cpCounts,cpCooks,low_count))
  compareIdentity(paste0("cooks_mask_",low_count,"_",cutoff),cpInput(fit,cutoff,"mask",low_count),
    list(pvalue=hexVector(actual$pvalue),cooksOutlier=as.list(actual$cooksOutlier)))
}
cpVariants <- list(matching=list(design=cpDesign,colData=data.frame(group=cpGroups)),
  releveled=list(design=cpDesign,colData=data.frame(group=relevel(cpGroups,"B"))),
  mismatched=list(design=cpDesign,colData=data.frame(group=factor(c("A","B","A","B","B")))),
  additive=list(design=cbind(cpDesign,c(0,1,0,1,0)),colData=data.frame(group=cpGroups)),
  missing=list(design=cpDesign,colData=data.frame(group=factor(c("A",NA,"B","B","B")))),
  third_level=list(design=cpDesign,colData=data.frame(group=factor(c("A","A","B","B","C")))))
for(name in names(cpVariants)) for(enabled in if(name %in% c("matching","additive"))c(FALSE,TRUE)else TRUE) {
  fit <- modifyList(cpBase,cpVariants[[name]]); cutoff <- if(enabled)5 else FALSE
  actual <- withoutRCalculations(rsdeseq2:::.rsdeseq2_fit_cooks_policy(cpTable,fit,cutoff),projection=TRUE,zeros=TRUE)
  actual <- c(lapply(as.list(actual[c("log2FoldChange","lfcSE","stat","pvalue")]),hexVector),
    list(cooksOutlier=lapply(actual$cooksOutlier,function(v)if(is.na(v))NULL else v)))
  compareIdentity(paste0("cooks_present_",name,"_",enabled),cpInput(fit,if(enabled)5 else NULL,"present"),actual)
}

# Strict formula-origin policies share the same original-DESeq2 oracle as the
# standalone and installed R tests. No R numerical routines may run in policy.
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[[1L]])
script <- gsub("~+~"," ",script,fixed=TRUE)
reference <- jsonlite::fromJSON(file.path(dirname(normalizePath(script)),"..","..",
  "crates/rsdeseq2/tests/data/cooks_policy_r462.json"),simplifyVector=FALSE)
decodeProbability <- function(x) vapply(x,function(v) {
  if(is.null(v)) return(NA_real_)
  readBin(as.raw(strtoi(substring(v,seq(1,15,2),seq(2,16,2)),16L)),"double",n=1,size=8,endian="big")
},numeric(1))
for(case in reference$cases) {
  matrixRows <- function(x,ncol) matrix(as.double(unlist(x)),ncol=ncol,byrow=TRUE)
  fields <- lapply(case$factors,function(f)factor(as.integer(unlist(f$codes)),levels=seq_len(f$levels)))
  names(fields) <- vapply(case$factors,`[[`,character(1),"name")
  p <- decodeProbability(case$raw_pvalue)
  table <- data.frame(log2FoldChange=rep(1,reference$genes),lfcSE=rep(.5,reference$genes),stat=rep(2,reference$genes),pvalue=p,cooksOutlier=rep(FALSE,reference$genes))
  fit <- list(test="Wald",schemaVersion=2L,
    designOrigin=list(kind=case$kind,variables=as.character(unlist(case$variables))),
    design=matrixRows(case$design,case$coefficients),colData=as.data.frame(fields),
    counts=matrixRows(reference$counts,reference$samples),cooks=matrixRows(reference$cooks,reference$samples),
    results=list(maxCooks=as.double(unlist(reference$maximum))))
  actual <- withoutRCalculations(rsdeseq2:::.rsdeseq2_fit_cooks_policy(table,fit,5),projection=TRUE,zeros=TRUE)
  stopifnot(identical(hexVector(actual$pvalue),hexVector(decodeProbability(case$pvalue))))
  input <- list(cooks_policy=list(mode="present",samples=reference$samples,
    origin=list(kind=case$kind,variables=case$variables),
    pvalue=lapply(p,function(v)if(is.na(v))NULL else sprintf("%.17g",v)),
    maximum=reference$maximum,effect=strings(table$log2FoldChange),se=strings(table$lfcSE),stat=strings(table$stat),cutoff="5",
    counts=reference$counts,cooks=reference$cooks,design=case$design,coefficients=case$coefficients,
    factors=case$factors,low_count=NULL,refitted=as.list(rep(FALSE,reference$genes)),new_zero=as.list(rep(FALSE,reference$genes))))
  compareIdentity(paste0("cooks_origin_",case$name),input,c(
    lapply(as.list(actual[c("log2FoldChange","lfcSE","stat","pvalue")]),hexVector),
    list(cooksOutlier=lapply(actual$cooksOutlier,function(v)if(is.na(v))NULL else v))))
}

# Offset products and normal-shrink preparation share the standalone kernels.
preparationReference <- jsonlite::fromJSON(file.path(dirname(normalizePath(script)),"..","..",
  "crates/rsdeseq2/tests/data/row_preparation_r462.json"),simplifyVector=FALSE)
preparationNumbers <- function(values)vapply(values,function(value)
  readBin(as.raw(strtoi(substring(value,seq(1,15,2),seq(2,16,2)),16L)),"double",n=1,size=8,endian="big"),numeric(1))
for(case in preparationReference$offsets) {
  counts <- matrix(1,case$genes,case$samples)
  sf <- if(is.null(case$size_factors))NULL else preparationNumbers(case$size_factors)
  nf <- if(is.null(case$factors))NULL else matrix(preparationNumbers(case$factors),case$genes,case$samples,byrow=TRUE)
  actual <- tryCatch(withoutRCalculations(.Call("rsdeseq2_offset_factors",counts,sf,nf,
    preparationNumbers(case$offsets),PACKAGE="rsdeseq2"),projection=TRUE,zeros=TRUE),error=identity)
  stopifnot(inherits(actual,"error")==is.null(case$expected))
  actual <- if(inherits(actual,"error"))list(error=TRUE)else list(factors=hexVector(as.vector(t(actual))))
  compareIdentity(paste0("offset_",case$name),list(preparation=c(list(kind="Offset"),case)),actual)
}
for(case in preparationReference$rows) {
  counts <- matrix(as.double(unlist(case$counts)),case$genes,case$samples,byrow=TRUE)
  actual <- tryCatch(withoutRCalculations(.Call("rsdeseq2_normal_shrinkage_rows",counts,
    preparationNumbers(case$dispersions),preparationNumbers(case$fitted),
    if(is.null(case$excluded))NULL else as.logical(unlist(case$excluded)),PACKAGE="rsdeseq2"),projection=TRUE,zeros=TRUE),error=identity)
  stopifnot(inherits(actual,"error")==is.null(case$expected))
  actual <- if(inherits(actual,"error"))list(error=TRUE)else list(indices=as.list(actual))
  compareIdentity(paste0("shrink_rows_",case$name),list(preparation=c(list(kind="Rows"),case)),actual)
}

# Model-frame expression evaluation is separate; matrix coding and products
# must execute wholly in the same standalone kernel for these evaluated inputs.
testSources <- file.path(dirname(normalizePath(script)),"..","..","r-pkg/rsdeseq2/tests/testthat")
source(file.path(testSources,"helper-reference-values.R"))
source(file.path(testSources,"helper-model-matrix.R"))
matrixReference <- jsonlite::fromJSON(.formula_reference_path(file.path(testSources,"data"),"model_matrix_r462.json"),simplifyVector=FALSE)
for(case in matrixReference$cases) {
  frame <- .model_reference_frame(case)
  actual <- suppressWarnings(withoutRCalculations(rsdeseq2:::.rsdeseq2_model_matrix(~1,frame),projection=TRUE,zeros=TRUE,model=TRUE))
  compareIdentity(paste0("model_matrix_",case$name),case,
    list(names=as.list(colnames(actual)),assignments=as.list(attr(actual,"assign")),values=hexVector(as.vector(t(actual)))))
}
for(case in matrixReference$average_contrasts) {
  positive<-as.integer(unlist(case$positive))+1L;negative<-as.integer(unlist(case$negative))+1L
  actual<-tryCatch(withoutRCalculations(.Call("rsdeseq2_average_contrast",as.integer(case$coefficients),positive,negative,PACKAGE="rsdeseq2"),projection=TRUE,zeros=TRUE,model=TRUE),error=identity)
  stopifnot(inherits(actual,"error")==is.null(case$expected))
  actual<-if(inherits(actual,"error"))list(error=TRUE)else list(values=hexVector(actual))
  compareIdentity(paste0("average_contrast_",case$name),list(preparation=c(list(kind="AverageContrast"),case)),actual)
}
