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
    if (fit_type == "mean") {
      native_fit <- rsdeseq2::DESeq(rust_dds, fitType = "mean", minReplicatesForReplace = Inf)
      reused <- rsdeseq2::rlog(native_fit, blind = FALSE, fitType = "local")
      testthat::expect_lt(max(abs(SummarizedExperiment::assay(reused) - SummarizedExperiment::assay(reference_rlog))), 0.007)
      testthat::expect_equal(S4Vectors::metadata(reused)$rsdeseq2Transform$fitType, "mean")
      changed <- native_fit
      SummarizedExperiment::assay(changed, "counts")[1L, 1L] <- fixture$counts[1L, 1L] + 1L
      testthat::expect_error(rsdeseq2::rlog(changed, blind = FALSE), "counts changed")
      changed <- native_fit
      S4Vectors::mcols(changed)$dispFit[1L] <- S4Vectors::mcols(changed)$dispFit[1L] * 2
      testthat::expect_error(rsdeseq2::rlog(changed, blind = FALSE), "dispFit changed")
    }
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
