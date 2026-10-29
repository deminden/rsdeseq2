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
withoutRCalculations <- function(expression) {
  # Fail at runtime if the production wrapper reintroduces statistical R calls.
  functions <- list(stats = c("p.adjust", "quantile", "lowess"), base = "mean")
  on.exit(for (package in names(functions)) for (name in functions[[package]])
    suppressMessages(untrace(name, where = asNamespace(package))), add = TRUE)
  for (package in names(functions)) for (name in functions[[package]])
    suppressMessages(trace(name, where = asNamespace(package), print = FALSE,
      tracer = quote(stop("Production wrapper executed an R numerical function", call. = FALSE))))
  force(expression)
}
guard <- tryCatch(withoutRCalculations(stats::p.adjust(c(.1, .2))), error = identity)
stopifnot(inherits(guard, "error"), grepl("Production wrapper executed", conditionMessage(guard)))
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
