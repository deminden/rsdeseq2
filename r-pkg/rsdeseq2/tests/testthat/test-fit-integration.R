testthat::test_that("Rust DESeq fit exposes coefficients for name and contrast results", {
  counts <- .rsdeseq2_test_counts()
  condition <- factor(rep(c("A", "B"), each = 4), levels = c("A", "B"))
  design <- stats::model.matrix(~ condition)

  fit <- DESeq(list(counts = counts, design = design), fitType = "mean")
  default <- results(fit)
  by_name <- results(fit, name = "conditionB", independentFiltering = FALSE)
  by_contrast <- results(fit, contrast = c(0, 1), independentFiltering = FALSE)

  testthat::expect_s3_class(fit, "rsdeseq2RustFit")
  testthat::expect_equal(dim(fit$beta), c(nrow(counts), ncol(design)))
  testthat::expect_equal(dim(fit$betaCovariance), c(nrow(counts), ncol(design), ncol(design)))
  testthat::expect_equal(names(fit$sizeFactors), colnames(counts))
  testthat::expect_equal(rownames(default), rownames(counts))
  testthat::expect_named(default, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "dispersion", "converged", "maxCooks", "cooksOutlier", "filtered"))
  testthat::expect_equal(unname(by_name$log2FoldChange), unname(by_contrast$log2FoldChange))
  testthat::expect_equal(unname(by_name$stat), unname(by_contrast$stat))
  testthat::expect_named(by_name, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "filtered"))
  testthat::expect_equal(rownames(by_name), rownames(counts))

  supplied_size_factors <- stats::setNames(rep(1, ncol(counts)), colnames(counts))
  supplied_fit <- DESeq(
    list(counts = counts, design = design),
    fitType = "mean",
    sizeFactors = supplied_size_factors
  )
  testthat::expect_equal(unname(supplied_fit$sizeFactors), unname(supplied_size_factors))

  normalization_factors <- matrix(
    1,
    nrow(counts),
    ncol(counts),
    dimnames = list(rev(rownames(counts)), rev(colnames(counts)))
  )
  observation_weights <- normalization_factors
  normalized_fit <- DESeq(
    list(counts = counts, design = design),
    fitType = "mean",
    normalizationFactors = normalization_factors,
    weights = observation_weights
  )
  testthat::expect_equal(unname(normalized_fit$normalizationFactors), matrix(1, nrow(counts), ncol(counts)))
  testthat::expect_equal(unname(normalized_fit$weights), matrix(1, nrow(counts), ncol(counts)))
  testthat::expect_equal(dimnames(normalized_fit$normalizationFactors), dimnames(counts))
  testthat::expect_equal(dimnames(normalized_fit$weights), dimnames(counts))
})

testthat::test_that("Rust LRT fit returns model-comparison results and selected coefficients", {
  counts <- .rsdeseq2_test_counts()
  condition <- factor(rep(c("A", "B"), each = 4), levels = c("A", "B"))
  data <- data.frame(condition = condition)

  fit <- DESeq(
    list(counts = counts, design = stats::model.matrix(~ condition)),
    fitType = "mean",
    test = "LRT",
    reduced = stats::model.matrix(~ 1, data = data)
  )
  default <- results(fit, independentFiltering = FALSE)
  selected <- results(fit, name = "conditionB", independentFiltering = FALSE)

  testthat::expect_equal(fit$lrtDf, 1L)
  testthat::expect_true(all(is.finite(default$pvalue)))
  testthat::expect_true(all(is.finite(selected$log2FoldChange)))
  testthat::expect_equal(unname(default$stat), unname(fit$lrtStat))
})

testthat::test_that("Rust DESeq accepts a Bioconductor DESeqDataSet", {
  testthat::skip_if_not_installed("DESeq2")
  counts <- .rsdeseq2_test_counts()
  col_data <- data.frame(
    condition = factor(rep(c("A", "B"), each = 4), levels = c("A", "B")),
    exposure = seq(0.8, 1.2, length.out = ncol(counts)),
    row.names = colnames(counts)
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(
    counts,
    col_data,
    design = ~ condition + offset(log(exposure))
  )
  supplied_norm_factors <- matrix(1.2, nrow(counts), ncol(counts), dimnames = dimnames(counts))
  DESeq2::normalizationFactors(dds) <- supplied_norm_factors
  SummarizedExperiment::assay(dds, "weights") <- matrix(1, nrow(counts), ncol(counts), dimnames = dimnames(counts))

  fit <- DESeq(dds, fitType = "mean")
  fit_state <- S4Vectors::metadata(fit)$rsdeseq2RustFit

  testthat::expect_true(inherits(fit, "DESeqDataSet"))
  testthat::expect_s3_class(fit_state, "rsdeseq2RustFit")
  testthat::expect_equal(resultsNames(fit), c("(Intercept)", "conditionB"))
  testthat::expect_equal(unname(DESeq2::sizeFactors(fit)), unname(fit_state$sizeFactors))
  expected_norm_factors <- supplied_norm_factors * matrix(
    rep(col_data$exposure, each = nrow(counts)),
    nrow = nrow(counts),
    ncol = ncol(counts)
  )
  testthat::expect_equal(unname(DESeq2::normalizationFactors(fit)), unname(expected_norm_factors))
  testthat::expect_equal(unname(fit_state$weights), matrix(1, nrow(counts), ncol(counts)))
  testthat::expect_equal(rownames(results(fit, name = "conditionB")), rownames(counts))
  shrunk <- lfcShrink(fit, coef = "conditionB")
  testthat::expect_true(all(is.finite(shrunk$log2FoldChange)))
  testthat::expect_true(all(is.finite(shrunk$lfcSE)))
})

testthat::test_that("Rlog dataset wrapper applies supplied frozen intercepts", {
  testthat::skip_if_not_installed("DESeq2")
  counts <- rbind(.rsdeseq2_test_counts(), all_zero = rep(0, ncol(.rsdeseq2_test_counts())))
  col_data <- data.frame(
    condition = factor(rep(c("A", "B"), each = 4), levels = c("A", "B")),
    row.names = colnames(counts)
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  supplied_intercept <- stats::setNames(seq(3.5, 4.5, length.out = nrow(counts)), rownames(counts))
  supplied_intercept <- rev(supplied_intercept)

  transformed <- rlog(
    dds,
    blind = TRUE,
    intercept = supplied_intercept,
    betaPriorVar = 0.5,
    fitType = "mean"
  )

  testthat::expect_s4_class(transformed, "DESeqTransform")
  testthat::expect_equal(
    S4Vectors::mcols(transformed)$rlogIntercept,
    unname(supplied_intercept[rownames(counts)])
  )
  testthat::expect_equal(
    unname(SummarizedExperiment::assay(transformed)["all_zero", ]),
    rep(unname(supplied_intercept[["all_zero"]]), ncol(counts))
  )
  testthat::expect_equal(attr(transformed, "betaPriorVar"), 0.5)
})

testthat::test_that("dataset VST and rlog return Bioconductor transform objects", {
  testthat::skip_if_not_installed("DESeq2")
  counts <- .rsdeseq2_test_counts()
  col_data <- data.frame(
    condition = factor(rep(c("A", "B"), each = 4), levels = c("A", "B")),
    row.names = colnames(counts)
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)

  fast <- vst(dds, blind = TRUE, nsub = 10L, fitType = "mean")
  full <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "mean")
  regularized <- rlog(dds, blind = FALSE, fitType = "mean")
  fast_default <- vst(dds, nsub = 10L)
  full_default <- varianceStabilizingTransformation(dds)
  regularized_default <- rlog(dds)

  for (transformed in list(fast, full, regularized, fast_default, full_default, regularized_default)) {
    testthat::expect_s4_class(transformed, "DESeqTransform")
    testthat::expect_equal(dim(SummarizedExperiment::assay(transformed)), dim(counts))
    testthat::expect_true(all(is.finite(SummarizedExperiment::assay(transformed))))
    testthat::expect_equal(rownames(transformed), rownames(counts))
    transformed_col_data <- as.data.frame(SummarizedExperiment::colData(transformed))
    testthat::expect_equal(transformed_col_data$condition, col_data$condition)
    testthat::expect_true(all(is.finite(transformed_col_data$sizeFactor)))
  }
  testthat::expect_equal(S4Vectors::metadata(fast)$rsdeseq2Transform$type, "vst")
  testthat::expect_equal(S4Vectors::metadata(full)$rsdeseq2Transform$type, "vstFull")
  testthat::expect_equal(S4Vectors::metadata(regularized)$rsdeseq2Transform$type, "rlog")
  reference_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  reference_full <- DESeq2::varianceStabilizingTransformation(
    reference_dds, blind = FALSE, fitType = "mean"
  )
  reference_rlog <- DESeq2::rlog(reference_dds, blind = FALSE, fitType = "mean")
  local_full <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "local")
  reference_local_full <- DESeq2::varianceStabilizingTransformation(
    reference_dds, blind = FALSE, fitType = "local"
  )
  local_fast <- vst(dds, blind = FALSE, nsub = 10L, fitType = "local")
  reference_local_fast <- DESeq2::vst(
    reference_dds, blind = FALSE, nsub = 10L, fitType = "local"
  )
  testthat::expect_equal(
    SummarizedExperiment::assay(full),
    SummarizedExperiment::assay(reference_full),
    tolerance = 4e-5
  )
  testthat::expect_equal(
    attr(regularized, "betaPriorVar"),
    attr(reference_rlog, "betaPriorVar"),
    tolerance = 1e-10
  )
  testthat::expect_lt(
    max(abs(
      SummarizedExperiment::assay(regularized) -
        SummarizedExperiment::assay(reference_rlog)
    )),
    0.015
  )
  testthat::expect_lt(
    max(abs(
      SummarizedExperiment::assay(local_full) -
        SummarizedExperiment::assay(reference_local_full)
    )),
    0.005
  )
  testthat::expect_lt(
    max(abs(
      SummarizedExperiment::assay(local_fast) -
        SummarizedExperiment::assay(reference_local_fast)
    )),
    0.005
  )
  for (transformed in list(fast_default, full_default, regularized_default)) {
    testthat::expect_equal(S4Vectors::metadata(transformed)$rsdeseq2Transform$fitTypeRequested, "parametric")
    testthat::expect_equal(S4Vectors::metadata(transformed)$rsdeseq2Transform$fitType, "local")
  }
})

testthat::test_that("dataset transforms retain DESeq2 parity on a large seeded fixture", {
  testthat::skip_if_not_installed("DESeq2")
  fixture <- .rsdeseq2_large_transform_fixture()
  rust_dds <- DESeq2::DESeqDataSetFromMatrix(
    fixture$counts, fixture$colData, design = ~ condition
  )
  reference_dds <- DESeq2::DESeqDataSetFromMatrix(
    fixture$counts, fixture$colData, design = ~ condition
  )

  tolerances <- list(
    mean = c(full = 0.002, fast = 0.002),
    parametric = c(full = 0.0025, fast = 0.0025),
    local = c(full = 0.006, fast = 0.004)
  )
  for (fit_type in names(tolerances)) {
    rust_full <- rsdeseq2::varianceStabilizingTransformation(
      rust_dds, blind = FALSE, fitType = fit_type
    )
    reference_full <- DESeq2::varianceStabilizingTransformation(
      reference_dds, blind = FALSE, fitType = fit_type
    )
    testthat::expect_lt(
      max(abs(
        SummarizedExperiment::assay(rust_full) -
          SummarizedExperiment::assay(reference_full)
      )),
      tolerances[[fit_type]][["full"]]
    )
    testthat::expect_equal(
      S4Vectors::metadata(rust_full)$rsdeseq2Transform$fitType,
      fit_type
    )

    rust_fast <- rsdeseq2::vst(
      rust_dds, blind = FALSE, nsub = 1000L, fitType = fit_type
    )
    reference_fast <- DESeq2::vst(
      reference_dds, blind = FALSE, nsub = 1000L, fitType = fit_type
    )
    testthat::expect_lt(
      max(abs(
        SummarizedExperiment::assay(rust_fast) -
          SummarizedExperiment::assay(reference_fast)
      )),
      tolerances[[fit_type]][["fast"]]
    )
    testthat::expect_equal(
      S4Vectors::metadata(rust_fast)$rsdeseq2Transform$fitType,
      fit_type
    )
  }

  rust_default_fast <- rsdeseq2::vst(rust_dds, nsub = 1000L)
  reference_default_fast <- DESeq2::vst(reference_dds, nsub = 1000L)
  testthat::expect_true(
    S4Vectors::metadata(rust_default_fast)$rsdeseq2Transform$blind
  )
  testthat::expect_lt(
    max(abs(
      SummarizedExperiment::assay(rust_default_fast) -
        SummarizedExperiment::assay(reference_default_fast)
    )),
    0.012
  )

  rust_default_full <- rsdeseq2::varianceStabilizingTransformation(rust_dds)
  reference_default_full <- DESeq2::varianceStabilizingTransformation(reference_dds)
  testthat::expect_true(
    S4Vectors::metadata(rust_default_full)$rsdeseq2Transform$blind
  )
  testthat::expect_lt(
    max(abs(
      SummarizedExperiment::assay(rust_default_full) -
        SummarizedExperiment::assay(reference_default_full)
    )),
    0.01
  )

  rust_default_rlog <- rsdeseq2::rlog(rust_dds)
  reference_default_rlog <- DESeq2::rlog(reference_dds)
  testthat::expect_true(
    S4Vectors::metadata(rust_default_rlog)$rsdeseq2Transform$blind
  )
  testthat::expect_lt(
    max(abs(
      SummarizedExperiment::assay(rust_default_rlog) -
        SummarizedExperiment::assay(reference_default_rlog)
    )),
    0.03
  )

  for (fit_type in names(tolerances)) {
    rust_rlog <- rsdeseq2::rlog(rust_dds, blind = FALSE, fitType = fit_type)
    reference_rlog <- DESeq2::rlog(reference_dds, blind = FALSE, fitType = fit_type)
    testthat::expect_lt(
      max(abs(
        SummarizedExperiment::assay(rust_rlog) -
          SummarizedExperiment::assay(reference_rlog)
      )),
      0.007
    )
    testthat::expect_equal(
      attr(rust_rlog, "betaPriorVar"),
      attr(reference_rlog, "betaPriorVar"),
      tolerance = 1e-10
    )
  }
})

testthat::test_that("full VST reuses frozen DESeq2 dispersion functions", {
  testthat::skip_if_not_installed("DESeq2")
  counts <- .rsdeseq2_test_counts()
  col_data <- data.frame(
    condition = factor(rep(c("A", "B"), each = 4), levels = c("A", "B")),
    row.names = colnames(counts)
  )
  for (fit_type in c("parametric", "mean", "local")) {
    rust_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
    reference_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
    DESeq2::sizeFactors(rust_dds) <- rep(1, ncol(counts))
    DESeq2::sizeFactors(reference_dds) <- rep(1, ncol(counts))
    if (fit_type == "parametric") {
      trend <- function(mu) 0.15 + 1.5 / mu
      attr(trend, "fitType") <- "parametric"
      attr(trend, "coefficients") <- c(extraPois = 1.5, asymptDisp = 0.15)
    } else if (fit_type == "mean") {
      trend <- function(mu) rep(0.2, length(mu))
      attr(trend, "fitType") <- "mean"
      attr(trend, "mean") <- 0.2
    } else {
      trend <- function(mu) 0.15 + 0.3 / sqrt(pmax(mu, 1e-8))
      attr(trend, "fitType") <- "local"
    }
    DESeq2::dispersionFunction(rust_dds) <- trend
    DESeq2::dispersionFunction(reference_dds) <- trend

    transformed <- rsdeseq2::varianceStabilizingTransformation(
      rust_dds, blind = FALSE, fitType = "local"
    )
    reference <- DESeq2::varianceStabilizingTransformation(
      reference_dds, blind = FALSE, fitType = "local"
    )
    testthat::expect_equal(
      SummarizedExperiment::assay(transformed),
      SummarizedExperiment::assay(reference),
      tolerance = 2e-6
    )
    transform_metadata <- S4Vectors::metadata(transformed)$rsdeseq2Transform
    testthat::expect_equal(transform_metadata$fitType, fit_type)
    testthat::expect_equal(transform_metadata$fitTypeRequested, "local")
    testthat::expect_equal(transform_metadata$trendSource, "frozen_dispersion_function")
    testthat::expect_equal(transform_metadata$frozenFitType, fit_type)
  }

  rust_nf_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  reference_nf_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  normalization_factors <- outer(
    seq(0.8, 1.2, length.out = nrow(counts)),
    seq(0.9, 1.1, length.out = ncol(counts))
  )
  dimnames(normalization_factors) <- dimnames(counts)
  DESeq2::normalizationFactors(rust_nf_dds) <- normalization_factors
  DESeq2::normalizationFactors(reference_nf_dds) <- normalization_factors
  local_trend <- function(mu) 0.15 + 0.3 / sqrt(pmax(mu, 1e-8))
  attr(local_trend, "fitType") <- "local"
  DESeq2::dispersionFunction(rust_nf_dds) <- local_trend
  DESeq2::dispersionFunction(reference_nf_dds) <- local_trend
  transformed_nf <- rsdeseq2::varianceStabilizingTransformation(
    rust_nf_dds, blind = FALSE, fitType = "local"
  )
  reference_nf <- DESeq2::varianceStabilizingTransformation(
    reference_nf_dds, blind = FALSE, fitType = "local"
  )
  testthat::expect_equal(
    SummarizedExperiment::assay(transformed_nf),
    SummarizedExperiment::assay(reference_nf),
    tolerance = 2e-6
  )
})

testthat::test_that("design-aware rlog reuses stored DESeq2 dispFit values", {
  testthat::skip_if_not_installed("DESeq2")
  counts <- .rsdeseq2_test_counts()
  col_data <- data.frame(
    condition = factor(rep(c("A", "B"), each = 4), levels = c("A", "B")),
    row.names = colnames(counts)
  )
  rust_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  reference_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  rust_dds <- DESeq2::estimateSizeFactors(rust_dds)
  reference_dds <- DESeq2::estimateSizeFactors(reference_dds)
  rust_dds <- DESeq2::estimateDispersions(rust_dds, fitType = "mean", quiet = TRUE)
  reference_dds <- DESeq2::estimateDispersions(reference_dds, fitType = "mean", quiet = TRUE)

  transformed <- rsdeseq2::rlog(rust_dds, blind = FALSE, fitType = "local")
  reference <- DESeq2::rlog(reference_dds, blind = FALSE)
  testthat::expect_lt(
    max(abs(
      as.vector(SummarizedExperiment::assay(transformed)) -
        as.vector(SummarizedExperiment::assay(reference))
    )),
    0.02
  )
  testthat::expect_equal(
    attr(transformed, "betaPriorVar"),
    attr(reference, "betaPriorVar"),
    tolerance = 1e-8
  )
  transform_metadata <- S4Vectors::metadata(transformed)$rsdeseq2Transform
  testthat::expect_equal(transform_metadata$trendSource, "frozen_dispFit")
  testthat::expect_equal(transform_metadata$fitType, "mean")
})

testthat::test_that("dataset transformations preserve normalization offsets and weights", {
  testthat::skip_if_not_installed("DESeq2")
  counts <- .rsdeseq2_test_counts()
  exposure <- seq(0.8, 1.2, length.out = ncol(counts))
  col_data <- data.frame(
    condition = factor(rep(c("A", "B"), each = 4), levels = c("A", "B")),
    exposure = exposure,
    row.names = colnames(counts)
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(
    counts, col_data, design = ~ condition + offset(log(exposure))
  )
  supplied_factors <- matrix(
    1.2, nrow(counts), ncol(counts),
    dimnames = list(rev(rownames(counts)), rev(colnames(counts)))
  )
  SummarizedExperiment::assay(dds, "weights", withDimnames = FALSE) <- matrix(
    1, nrow(counts), ncol(counts),
    dimnames = list(rev(rownames(counts)), rev(colnames(counts)))
  )
  DESeq2::normalizationFactors(dds) <- supplied_factors

  transformed <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "mean")
  expected_factors <- supplied_factors * matrix(
    rep(exposure, each = nrow(counts)), nrow = nrow(counts), ncol = ncol(counts)
  )

  testthat::expect_equal(
    unname(S4Vectors::metadata(transformed)$rsdeseq2Transform$normalizationFactors),
    unname(expected_factors)
  )
  testthat::expect_equal(
    S4Vectors::metadata(transformed)$rsdeseq2Transform$formulaOffsets,
    log(exposure)
  )
  testthat::expect_true(S4Vectors::metadata(transformed)$rsdeseq2Transform$weightsUsed)
  testthat::expect_true(all(is.finite(SummarizedExperiment::assay(transformed))))
})

testthat::test_that("dataset fast VST matches weighted subset and eligibility behavior", {
  testthat::skip_if_not_installed("DESeq2")
  counts <- .rsdeseq2_test_counts()
  col_data <- data.frame(
    condition = factor(rep(c("A", "B"), each = 4), levels = c("A", "B")),
    row.names = colnames(counts)
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  weights <- matrix(0.2, nrow(counts), ncol(counts), dimnames = dimnames(counts))
  for (gene in seq_len(nrow(weights))) {
    weights[gene, ((gene - 1L) %% ncol(weights)) + 1L] <- 1
  }
  SummarizedExperiment::assay(dds, "weights") <- weights

  observed <- vst(dds, blind = TRUE, nsub = 10L, fitType = "mean")
  reference <- DESeq2::vst(dds, blind = TRUE, nsub = 10L, fitType = "mean")
  testthat::expect_equal(
    SummarizedExperiment::assay(observed),
    SummarizedExperiment::assay(reference),
    tolerance = 5e-5
  )

  small_dds <- DESeq2::DESeqDataSetFromMatrix(
    counts[seq_len(5L), , drop = FALSE],
    col_data,
    design = ~ condition
  )
  testthat::expect_error(vst(small_dds, nsub = 10L, fitType = "mean"))
  testthat::expect_error(DESeq2::vst(small_dds, nsub = 10L, fitType = "mean"))
})

testthat::test_that("normal lfcShrink returns Rust beta-prior estimates for a dataset coefficient", {
  testthat::skip_if_not_installed("DESeq2")
  counts <- rbind(.rsdeseq2_test_counts(), all_zero = rep(0L, 8L))
  col_data <- data.frame(
    condition = factor(rep(c("A", "B"), each = 4), levels = c("A", "B")),
    row.names = colnames(counts)
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  fitted <- rsdeseq2::DESeq(dds, fitType = "mean")
  shrunk_by_name <- rsdeseq2::lfcShrink(fitted, coef = "conditionB")
  shrunk_by_index <- rsdeseq2::lfcShrink(fitted, coef = 2L)
  thresholded <- rsdeseq2::lfcShrink(
    fitted, coef = "conditionB", lfcThreshold = 0.5,
    independentFiltering = FALSE
  )
  mle <- rsdeseq2::results(fitted, name = "conditionB", independentFiltering = FALSE)
  reference_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  reference_dds <- DESeq2::DESeq(reference_dds, fitType = "mean", betaPrior = FALSE, quiet = TRUE)
  reference <- DESeq2::lfcShrink(
    reference_dds,
    coef = utils::tail(DESeq2::resultsNames(reference_dds), 1L),
    type = "normal"
  )
  reference_thresholded <- DESeq2::lfcShrink(
    reference_dds,
    coef = utils::tail(DESeq2::resultsNames(reference_dds), 1L),
    type = "normal",
    lfcThreshold = 0.5
  )

  testthat::expect_s4_class(shrunk_by_name, "DESeqResults")
  testthat::expect_equal(rownames(shrunk_by_name), rownames(counts))
  testthat::expect_equal(
    unname(shrunk_by_name$log2FoldChange),
    unname(shrunk_by_index$log2FoldChange)
  )
  nonzero <- rownames(counts) != "all_zero"
  testthat::expect_true(all(is.finite(shrunk_by_name$log2FoldChange[nonzero])))
  testthat::expect_true(all(is.finite(shrunk_by_name$lfcSE[nonzero])))
  testthat::expect_true(is.na(shrunk_by_name$log2FoldChange[!nonzero]))
  testthat::expect_true(any(abs(shrunk_by_name$log2FoldChange) < abs(mle$log2FoldChange)))
  testthat::expect_equal(shrunk_by_name$log2FoldChange, reference$log2FoldChange, tolerance = 3e-5)
  testthat::expect_equal(shrunk_by_name$lfcSE, reference$lfcSE, tolerance = 5e-6)
  testthat::expect_equal(thresholded$log2FoldChange, reference_thresholded$log2FoldChange, tolerance = 3e-5)
  testthat::expect_equal(thresholded$lfcSE, reference_thresholded$lfcSE, tolerance = 5e-6)
  testthat::expect_equal(thresholded$stat, reference_thresholded$stat, tolerance = 1.2e-4)
  testthat::expect_equal(thresholded$pvalue, reference_thresholded$pvalue, tolerance = 4.2e-5)
  testthat::expect_equal(
    thresholded$padj,
    stats::p.adjust(thresholded$pvalue, method = "BH")
  )
  testthat::expect_equal(DESeq2::priorInfo(shrunk_by_name)$type, "normal")
  testthat::expect_true(is.finite(DESeq2::priorInfo(shrunk_by_name)$betaPriorVar[["conditionB"]]))
})

testthat::test_that("ashr lfcShrink matches DESeq2 posterior summaries and s-values", {
  testthat::skip_if_not_installed("DESeq2")
  testthat::skip_if_not_installed("ashr")
  fixture <- .rsdeseq2_large_transform_fixture()
  rust_dds <- DESeq2::DESeqDataSetFromMatrix(
    fixture$counts, fixture$colData, design = ~ condition
  )
  reference_dds <- DESeq2::DESeqDataSetFromMatrix(
    fixture$counts, fixture$colData, design = ~ condition
  )
  rust_fit <- rsdeseq2::DESeq(rust_dds, fitType = "mean")
  reference_fit <- DESeq2::DESeq(reference_dds, fitType = "mean", quiet = TRUE)
  reference_coef <- utils::tail(DESeq2::resultsNames(reference_fit), 1L)
  rust_shrunk <- rsdeseq2::lfcShrink(
    rust_fit, coef = "conditionB", type = "ashr", independentFiltering = FALSE
  )
  reference_shrunk <- DESeq2::lfcShrink(
    reference_fit, coef = reference_coef, type = "ashr"
  )
  rust_svalues <- rsdeseq2::lfcShrink(
    rust_fit, coef = "conditionB", type = "ashr", svalue = TRUE,
    independentFiltering = FALSE
  )
  reference_svalues <- DESeq2::lfcShrink(
    reference_fit, coef = reference_coef, type = "ashr", svalue = TRUE
  )
  rust_thresholded <- rsdeseq2::lfcShrink(
    rust_fit, coef = "conditionB", type = "ashr", lfcThreshold = 0.5,
    independentFiltering = FALSE
  )
  reference_thresholded <- DESeq2::lfcShrink(
    reference_fit, coef = reference_coef, type = "ashr", lfcThreshold = 0.5
  )

  testthat::expect_s4_class(rust_shrunk, "DESeqResults")
  testthat::expect_equal(
    rust_shrunk$log2FoldChange, reference_shrunk$log2FoldChange,
    tolerance = 1e-3
  )
  testthat::expect_equal(rust_shrunk$lfcSE, reference_shrunk$lfcSE, tolerance = 1e-3)
  testthat::expect_equal(rust_svalues$svalue, reference_svalues$svalue, tolerance = 0.02)
  testthat::expect_equal(
    rust_thresholded$svalue, reference_thresholded$svalue, tolerance = 0.02
  )
  testthat::expect_equal(
    rust_thresholded$log2FoldChange,
    reference_thresholded$log2FoldChange,
    tolerance = 1e-3
  )
  testthat::expect_equal(DESeq2::priorInfo(rust_shrunk)$type, "ashr")
  testthat::expect_true(!is.null(DESeq2::priorInfo(rust_shrunk)$fitted_g))
})

testthat::test_that("apeglm lfcShrink matches DESeq2 coefficient posteriors and FSOS s-values", {
  testthat::skip_if_not_installed("DESeq2")
  testthat::skip_if_not_installed("apeglm")
  fixture <- .rsdeseq2_large_transform_fixture()
  rust_dds <- DESeq2::DESeqDataSetFromMatrix(
    fixture$counts, fixture$colData, design = ~ condition
  )
  reference_dds <- DESeq2::DESeqDataSetFromMatrix(
    fixture$counts, fixture$colData, design = ~ condition
  )
  rust_fit <- rsdeseq2::DESeq(rust_dds, fitType = "mean")
  reference_fit <- DESeq2::DESeq(reference_dds, fitType = "mean", quiet = TRUE)
  reference_coef <- utils::tail(DESeq2::resultsNames(reference_fit), 1L)
  rust_shrunk <- rsdeseq2::lfcShrink(
    rust_fit, coef = "conditionB", type = "apeglm", quiet = TRUE
  )
  reference_shrunk <- DESeq2::lfcShrink(
    reference_fit, coef = reference_coef, type = "apeglm", quiet = TRUE
  )
  rust_svalues <- rsdeseq2::lfcShrink(
    rust_fit, coef = "conditionB", type = "apeglm", svalue = TRUE,
    quiet = TRUE
  )
  reference_svalues <- DESeq2::lfcShrink(
    reference_fit, coef = reference_coef, type = "apeglm", svalue = TRUE,
    quiet = TRUE
  )
  rust_thresholded <- rsdeseq2::lfcShrink(
    rust_fit, coef = "conditionB", type = "apeglm", lfcThreshold = 0.5,
    quiet = TRUE
  )
  reference_thresholded <- DESeq2::lfcShrink(
    reference_fit, coef = reference_coef, type = "apeglm",
    lfcThreshold = 0.5, quiet = TRUE
  )

  testthat::expect_s4_class(rust_shrunk, "DESeqResults")
  testthat::expect_equal(
    unname(rust_shrunk$log2FoldChange), unname(reference_shrunk$log2FoldChange),
    tolerance = 0.01
  )
  testthat::expect_equal(
    unname(rust_shrunk$lfcSE), unname(reference_shrunk$lfcSE), tolerance = 0.01
  )
  testthat::expect_equal(rust_svalues$svalue, reference_svalues$svalue, tolerance = 0.02)
  testthat::expect_equal(
    unname(rust_thresholded$log2FoldChange),
    unname(reference_thresholded$log2FoldChange),
    tolerance = 0.01
  )
  testthat::expect_equal(
    rust_thresholded$svalue, reference_thresholded$svalue, tolerance = 0.05
  )
  testthat::expect_equal(
    names(rust_shrunk), c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
  )
  testthat::expect_equal(DESeq2::priorInfo(rust_shrunk)$type, "apeglm")
  testthat::expect_true(!is.null(DESeq2::priorInfo(rust_shrunk)$prior.control))
})

testthat::test_that("normal lfcShrink applies an expanded Rust fit to factor contrasts", {
  testthat::skip_if_not_installed("DESeq2")
  counts <- .rsdeseq2_test_counts()
  col_data <- data.frame(
    condition = factor(rep(c("A", "B"), each = 4L), levels = c("A", "B")),
    row.names = colnames(counts)
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  fitted <- rsdeseq2::DESeq(dds, fitType = "mean")
  shrunk <- rsdeseq2::lfcShrink(
    fitted, contrast = c("condition", "B", "A"), independentFiltering = FALSE
  )
  thresholded <- rsdeseq2::lfcShrink(
    fitted, contrast = c("condition", "B", "A"), lfcThreshold = 0.5,
    independentFiltering = FALSE
  )
  reference_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  reference_dds <- DESeq2::DESeq(reference_dds, fitType = "mean", quiet = TRUE)
  reference <- DESeq2::lfcShrink(
    reference_dds, contrast = c("condition", "B", "A"), type = "normal"
  )
  reference_thresholded <- DESeq2::lfcShrink(
    reference_dds, contrast = c("condition", "B", "A"), type = "normal",
    lfcThreshold = 0.5
  )
  unshrunk <- rsdeseq2::results(
    fitted, contrast = c("condition", "B", "A"),
    independentFiltering = FALSE
  )

  testthat::expect_s4_class(shrunk, "DESeqResults")
  testthat::expect_equal(shrunk$log2FoldChange, reference$log2FoldChange, tolerance = 1e-4)
  testthat::expect_equal(shrunk$lfcSE, reference$lfcSE, tolerance = 2e-5)
  testthat::expect_equal(shrunk$stat, unshrunk$stat, tolerance = 1e-10)
  testthat::expect_equal(shrunk$pvalue, unshrunk$pvalue, tolerance = 1e-10)
  testthat::expect_equal(thresholded$log2FoldChange, reference_thresholded$log2FoldChange, tolerance = 1e-4)
  testthat::expect_equal(thresholded$lfcSE, reference_thresholded$lfcSE, tolerance = 2e-5)
  testthat::expect_equal(thresholded$stat, reference_thresholded$stat, tolerance = 2e-4)
  testthat::expect_equal(thresholded$pvalue, reference_thresholded$pvalue, tolerance = 1e-4)
  testthat::expect_equal(
    names(DESeq2::priorInfo(shrunk)$betaPriorVar),
    c("Intercept", "conditionA", "conditionB")
  )
})

testthat::test_that("normal contrast shrinkage estimates factor priors across levels and covariates", {
  testthat::skip_if_not_installed("DESeq2")
  set.seed(41)
  n_genes <- 80L
  col_data <- data.frame(
    condition = factor(rep(c("A", "B", "C"), each = 4L), levels = c("A", "B", "C")),
    batch = factor(rep(c("x", "y"), 6L), levels = c("x", "y"))
  )
  rownames(col_data) <- paste0("sample", seq_len(nrow(col_data)))
  base_mean <- exp(seq(log(20), log(800), length.out = n_genes))
  fold_change <- rep(c(1, 1.4, 0.7, 2.2), length.out = n_genes)
  counts <- t(vapply(seq_len(n_genes), function(gene) {
    stats::rnbinom(12L, mu = base_mean[[gene]] * rep(c(1, fold_change[[gene]], 1 / fold_change[[gene]]), each = 4L), size = 5)
  }, numeric(12L)))
  rownames(counts) <- paste0("gene", seq_len(n_genes))
  colnames(counts) <- rownames(col_data)
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ batch + condition)
  fitted <- rsdeseq2::DESeq(dds, fitType = "mean")
  shrunk <- rsdeseq2::lfcShrink(
    fitted, contrast = c("condition", "C", "A"), independentFiltering = FALSE
  )
  thresholded <- rsdeseq2::lfcShrink(
    fitted, contrast = c("condition", "C", "A"), lfcThreshold = 0.5,
    independentFiltering = FALSE
  )
  reference_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ batch + condition)
  reference_dds <- DESeq2::DESeq(reference_dds, fitType = "mean", quiet = TRUE)
  reference <- DESeq2::lfcShrink(
    reference_dds, contrast = c("condition", "C", "A"), type = "normal"
  )
  reference_thresholded <- DESeq2::lfcShrink(
    reference_dds, contrast = c("condition", "C", "A"), type = "normal",
    lfcThreshold = 0.5
  )

  testthat::expect_equal(shrunk$log2FoldChange, reference$log2FoldChange, tolerance = 2e-4)
  testthat::expect_equal(shrunk$lfcSE, reference$lfcSE, tolerance = 2e-6)
  testthat::expect_equal(thresholded$log2FoldChange, reference_thresholded$log2FoldChange, tolerance = 2e-4)
  testthat::expect_equal(thresholded$lfcSE, reference_thresholded$lfcSE, tolerance = 2e-6)
  testthat::expect_equal(thresholded$stat, reference_thresholded$stat, tolerance = 5e-4)
  testthat::expect_equal(thresholded$pvalue, reference_thresholded$pvalue, tolerance = 2e-4)
  testthat::expect_equal(
    unname(DESeq2::priorInfo(shrunk)$betaPriorVar),
    unname(DESeq2::priorInfo(reference)$betaPriorVar),
    tolerance = 2e-5
  )
})

testthat::test_that("normal lfcShrink accepts numeric expanded-design contrasts with supplied results", {
  testthat::skip_if_not_installed("DESeq2")
  counts <- .rsdeseq2_test_counts()
  col_data <- data.frame(
    condition = factor(rep(c("A", "B"), each = 4L), levels = c("A", "B")),
    row.names = colnames(counts)
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  fitted <- rsdeseq2::DESeq(dds, fitType = "mean")
  unshrunk <- rsdeseq2::results(
    fitted, contrast = c("condition", "B", "A"), independentFiltering = FALSE
  )
  numeric_contrast <- c(0, -1, 1)
  shrunk <- rsdeseq2::lfcShrink(
    fitted, contrast = numeric_contrast, res = unshrunk,
    independentFiltering = FALSE
  )
  thresholded <- rsdeseq2::lfcShrink(
    fitted, contrast = numeric_contrast, res = unshrunk,
    lfcThreshold = 0.5, independentFiltering = FALSE
  )
  reference_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  reference_dds <- DESeq2::DESeq(reference_dds, fitType = "mean", quiet = TRUE)
  reference_unshrunk <- DESeq2::results(
    reference_dds, contrast = c("condition", "B", "A"),
    independentFiltering = FALSE
  )
  reference <- DESeq2::lfcShrink(
    reference_dds, contrast = numeric_contrast, res = reference_unshrunk,
    type = "normal"
  )
  reference_thresholded <- DESeq2::lfcShrink(
    reference_dds, contrast = numeric_contrast, res = reference_unshrunk,
    type = "normal", lfcThreshold = 0.5
  )

  testthat::expect_equal(shrunk$log2FoldChange, reference$log2FoldChange, tolerance = 1e-4)
  testthat::expect_equal(shrunk$lfcSE, reference$lfcSE, tolerance = 2e-5)
  testthat::expect_equal(shrunk$stat, unshrunk$stat, tolerance = 1e-10)
  testthat::expect_equal(thresholded$log2FoldChange, reference_thresholded$log2FoldChange, tolerance = 1e-4)
  testthat::expect_equal(thresholded$lfcSE, reference_thresholded$lfcSE, tolerance = 2e-5)
  testthat::expect_equal(thresholded$stat, reference_thresholded$stat, tolerance = 2e-4)
  testthat::expect_equal(thresholded$pvalue, reference_thresholded$pvalue, tolerance = 1e-4)
})

testthat::test_that("normal lfcShrink resolves list contrasts on expanded coefficients", {
  testthat::skip_if_not_installed("DESeq2")
  counts <- .rsdeseq2_test_counts()
  col_data <- data.frame(
    condition = factor(rep(c("A", "B"), each = 4L), levels = c("A", "B")),
    row.names = colnames(counts)
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  fitted <- rsdeseq2::DESeq(dds, fitType = "mean")
  contrast <- list("conditionB", "conditionA")
  unshrunk <- rsdeseq2::results(
    fitted, contrast = c("condition", "B", "A"), independentFiltering = FALSE
  )
  shrunk <- rsdeseq2::lfcShrink(
    fitted, contrast = contrast, res = unshrunk, independentFiltering = FALSE
  )
  thresholded <- rsdeseq2::lfcShrink(
    fitted, contrast = contrast, res = unshrunk, lfcThreshold = 0.5,
    independentFiltering = FALSE
  )
  reference_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  reference_dds <- DESeq2::DESeq(reference_dds, fitType = "mean", quiet = TRUE)
  reference_unshrunk <- DESeq2::results(
    reference_dds, contrast = c("condition", "B", "A"),
    independentFiltering = FALSE
  )
  reference <- DESeq2::lfcShrink(
    reference_dds, contrast = contrast, res = reference_unshrunk, type = "normal"
  )
  reference_thresholded <- DESeq2::lfcShrink(
    reference_dds, contrast = contrast, res = reference_unshrunk,
    type = "normal", lfcThreshold = 0.5
  )

  testthat::expect_equal(shrunk$log2FoldChange, reference$log2FoldChange, tolerance = 1e-4)
  testthat::expect_equal(shrunk$lfcSE, reference$lfcSE, tolerance = 2e-5)
  testthat::expect_equal(thresholded$log2FoldChange, reference_thresholded$log2FoldChange, tolerance = 1e-4)
  testthat::expect_equal(thresholded$lfcSE, reference_thresholded$lfcSE, tolerance = 2e-5)
  testthat::expect_equal(thresholded$stat, reference_thresholded$stat, tolerance = 2e-4)
  testthat::expect_equal(thresholded$pvalue, reference_thresholded$pvalue, tolerance = 1e-4)
})

testthat::test_that("normal lfcShrink handles transformed-factor coefficients and matches DESeq2 contrast rejection", {
  testthat::skip_if_not_installed("DESeq2")
  counts <- .rsdeseq2_test_counts()
  col_data <- data.frame(
    condition = factor(rep(c("A", "B"), each = 4), levels = c("A", "B")),
    row.names = colnames(counts)
  )
  rust_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ factor(condition))
  reference_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ factor(condition))
  fitted <- rsdeseq2::DESeq(rust_dds, fitType = "mean")
  reference <- DESeq2::DESeq(reference_dds, fitType = "mean", quiet = TRUE)
  coefficient <- "factor.condition.B"

  shrunk <- rsdeseq2::lfcShrink(fitted, coef = coefficient)
  reference_shrunk <- DESeq2::lfcShrink(reference, coef = coefficient, type = "normal")
  testthat::expect_equal(shrunk$log2FoldChange, reference_shrunk$log2FoldChange, tolerance = 1e-4)
  testthat::expect_equal(shrunk$lfcSE, reference_shrunk$lfcSE, tolerance = 2e-5)
  testthat::expect_equal(shrunk$stat, reference_shrunk$stat, tolerance = 2e-4)

  testthat::expect_error(
    rsdeseq2::lfcShrink(fitted, contrast = c("condition", "B", "A")),
    "not found|could not be resolved|must name a factor"
  )
  testthat::expect_error(
    DESeq2::lfcShrink(reference, contrast = c("condition", "B", "A"), type = "normal"),
    "expecting condition_B_vs_A"
  )
})
