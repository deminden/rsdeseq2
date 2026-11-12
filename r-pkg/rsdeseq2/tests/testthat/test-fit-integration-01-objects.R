testthat::test_that("Rust DESeq fit exposes coefficients for name and contrast results", {
  counts <- .rsdeseq2_test_counts()
  condition <- factor(rep(c("A", "B"), each = 4), levels = c("A", "B"))
  design <- stats::model.matrix(~ condition)

  fit <- DESeq(list(counts = counts, design = design), fitType = "mean")
  default <- results(fit)
  by_name <- results(fit, name = "conditionB", independentFiltering = FALSE)
  by_contrast <- results(fit, contrast = c(0, 1), independentFiltering = FALSE)

  testthat::expect_s3_class(fit, "rsdeseq2RustFit")
  testthat::expect_null(fit$lrtDf)
  testthat::expect_equal(dim(fit$beta), c(nrow(counts), ncol(design)))
  testthat::expect_equal(dim(fit$betaCovariance), c(nrow(counts), ncol(design), ncol(design)))
  testthat::expect_equal(names(fit$sizeFactors), colnames(counts))
  testthat::expect_equal(rownames(default), rownames(counts))
  testthat::expect_named(default, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "dispersion", "converged", "maxCooks", "cooksOutlier", "filtered"))
  result_metadata <- attr(default, "metadata")
  testthat::expect_equal(result_metadata$resultName, "conditionB")
  testthat::expect_equal(result_metadata$test, "Wald")
  testthat::expect_true(result_metadata$independentFiltering$enabled)
  testthat::expect_equal(attr(default, "comparison"), "coefficient conditionB")
  by_name_metadata <- attr(by_name, "metadata")
  testthat::expect_equal(unname(by_name_metadata$contrast), c(0, 1))
  testthat::expect_false(by_name_metadata$cooksReplacement$enabled)
  testthat::expect_equal(by_name_metadata$cooksReplacement$minReplicatesForReplace, 7L)
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
  data <- data.frame(condition = condition, dose = rep(1:4, 2))

  for (df in 1:2) {
    fit <- DESeq(
      list(counts = counts, design = stats::model.matrix(if (df == 1L) ~ condition else ~ condition + dose, data)),
      fitType = "mean",
      test = "LRT",
      reduced = stats::model.matrix(~ 1, data = data)
    )
    default <- results(fit, independentFiltering = FALSE)
    selected <- results(fit, name = "conditionB", independentFiltering = FALSE)

    testthat::expect_identical(fit$lrtDf, df)
    testthat::expect_true(all(is.finite(default$pvalue)))
    testthat::expect_equal(unname(default$pvalue), stats::pchisq(default$stat, df, lower.tail = FALSE), tolerance = 2e-12)
    testthat::expect_true(all(is.finite(selected$log2FoldChange)))
    testthat::expect_equal(unname(default$stat), unname(fit$lrtStat))
  }
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
  testthat::expect_identical(SummarizedExperiment::assay(fit, "mu"), fit_state$mu)
  testthat::expect_identical(SummarizedExperiment::assay(fit, "cooks"), fit_state$cooks)
  testthat::expect_equal(dim(fit_state$mu), dim(counts))
  testthat::expect_equal(dim(fit_state$cooks), dim(counts))
  row_metadata <- S4Vectors::mcols(fit)
  base_metadata <- baseMetadataRust(
    counts,
    normalizationFactors = expected_norm_factors,
    weights = SummarizedExperiment::assay(fit, "weights")
  )
  testthat::expect_equal(row_metadata$baseMean, fit_state$results$baseMean)
  testthat::expect_equal(row_metadata$baseVar, base_metadata$baseVar)
  testthat::expect_equal(row_metadata$allZero, base_metadata$allZero)
  testthat::expect_equal(row_metadata$dispFit, fit_state$dispFit)
  testthat::expect_equal(row_metadata$dispersion, fit_state$dispersions)
  testthat::expect_equal(row_metadata$betaConv, fit_state$results$converged)
  testthat::expect_equal(row_metadata$maxCooks, fit_state$results$maxCooks)
  testthat::expect_equal(row_metadata$MLE_conditionB, fit_state$beta[, "conditionB"])
  testthat::expect_equal(row_metadata$SE_conditionB, fit_state$betaSE[, "conditionB"])
  testthat::expect_equal(rownames(results(fit, name = "conditionB")), rownames(counts))
  shrunk <- lfcShrink(fit, coef = "conditionB")
  testthat::expect_true(all(is.finite(shrunk$log2FoldChange)))
  testthat::expect_true(all(is.finite(shrunk$lfcSE)))
})

testthat::test_that("DESeqDataSet Cook's replacement refits match DESeq2", {
  testthat::skip_if_not_installed("DESeq2")
  for (scenario in c("mixed", "all_zero_only")) {
    counts <- .rsdeseq2_test_counts()
    if (scenario == "mixed") counts[1L, 8L] <- 2000L
    expectedRefit <- if(scenario=="mixed")2L else 1L
    expectedRefitRows <- if(scenario=="mixed")rownames(counts)[1L] else character()
    counts[2L, ] <- 0L
    counts[2L, 8L] <- 1000L
    col_data <- data.frame(
      condition = factor(rep(c("A", "B"), each = 4L), levels = c("A", "B")),
      row.names = colnames(counts)
    )
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, ~ condition)

    for (test in c("Wald", "LRT")) {
      reference_args <- list(
        object = dds,
        fitType = "mean",
        test = test,
        minReplicatesForReplace = 3L,
        quiet = TRUE
      )
      if (test == "LRT") reference_args$reduced <- ~ 1
      reference <- do.call(DESeq2::DESeq, reference_args)
      fitted <- DESeq(
        dds,
        fitType = "mean",
        test = test,
        reduced = if (test == "LRT") ~ 1 else NULL,
        minReplicatesForReplace = 3L
      )
      fit_state <- S4Vectors::metadata(fitted)$rsdeseq2RustFit
      reference_results <- DESeq2::results(reference, independentFiltering = FALSE)
      fitted_results <- results(fitted, independentFiltering = FALSE)

      testthat::expect_equal(fit_state$nRefit, expectedRefit)
      testthat::expect_true(methods::is(fitted_results, "DESeqResults"))
      testthat::expect_equal(S4Vectors::metadata(fitted_results)$resultName, "conditionB")
      testthat::expect_equal(S4Vectors::metadata(fitted_results)$cooksReplacement$nRefit, expectedRefit)
      testthat::expect_equal(S4Vectors::mcols(fitted_results)$type, rep("results", ncol(fitted_results)))
      testthat::expect_equal(S4Vectors::mcols(fitted_results)$description[[2L]], "log2 fold change (MLE): conditionB")
      testthat::expect_equal(S4Vectors::mcols(fitted)$replace, S4Vectors::mcols(reference)$replace)
      testthat::expect_equal(S4Vectors::mcols(fitted)$baseVar, S4Vectors::mcols(reference)$baseVar, tolerance = 1e-4)
      testthat::expect_equal(S4Vectors::mcols(fitted)$allZero, S4Vectors::mcols(reference)$allZero)
      testthat::expect_equal(
        unname(SummarizedExperiment::colData(fitted)$replaceable),
        unname(SummarizedExperiment::colData(reference)$replaceable)
      )
      testthat::expect_equal(
        SummarizedExperiment::assay(fitted, "replaceCounts"),
        SummarizedExperiment::assay(reference, "replaceCounts"),
        tolerance = 0
      )
      testthat::expect_equal(
        SummarizedExperiment::assay(fitted, "replaceCooks"),
        SummarizedExperiment::assay(reference, "replaceCooks"),
        tolerance = 1e-5
      )
      testthat::expect_equal(
        SummarizedExperiment::assay(fitted, "counts"),
        SummarizedExperiment::assay(dds, "counts")
      )
      testthat::expect_equal(fitted_results$baseMean, reference_results$baseMean, tolerance = 1e-10)
      testthat::expect_equal(fitted_results$log2FoldChange, reference_results$log2FoldChange, tolerance = 1e-4)
      testthat::expect_equal(fitted_results$lfcSE, reference_results$lfcSE, tolerance = 1e-4)
      testthat::expect_equal(unname(fitted_results$stat), unname(reference_results$stat), tolerance = 2e-4)
      testthat::expect_equal(unname(fitted_results$pvalue), unname(reference_results$pvalue), tolerance = 1e-4)
      result_metadata <- attr(fitted_results, "metadata")
      testthat::expect_true(result_metadata$cooksReplacement$enabled)
      testthat::expect_equal(result_metadata$cooksReplacement$nRefit, expectedRefit)
      testthat::expect_equal(result_metadata$cooksReplacement$refitRows, expectedRefitRows)
      testthat::expect_equal(result_metadata$cooksReplacement$newAllZeroRows, rownames(counts)[2L])
      for (assayName in c("replaceCounts", "replaceCooks")) {
        edited <- fitted
        SummarizedExperiment::assay(edited, assayName)[1, 1] <-
          SummarizedExperiment::assay(edited, assayName)[1, 1] + 1
        testthat::expect_error(results(edited), paste(assayName, "changed"))
      }

    }

    filtered_result <- results(fitted)
    filtered_metadata <- S4Vectors::metadata(filtered_result)
    testthat::expect_equal(filtered_metadata$alpha, 0.1)
    testthat::expect_identical(
      filtered_metadata$independentFiltering$enabled,
      attr(filtered_result, "independentFiltering")$enabled
    )
    if (filtered_metadata$independentFiltering$enabled) {
      testthat::expect_true(is.finite(filtered_metadata$filterThreshold))
      testthat::expect_true(is.data.frame(filtered_metadata$filterNumRej))
    }

    disabled <- DESeq(dds, fitType = "mean", minReplicatesForReplace = Inf)
    testthat::expect_false("replaceCounts" %in% SummarizedExperiment::assayNames(disabled))
    testthat::expect_false(attr(results(disabled), "metadata")$cooksReplacement$enabled)
    testthat::expect_error(DESeq(dds, fitType = "mean", minReplicatesForReplace = 2L), "at least 3")
  }
})

testthat::test_that("DESeqDataSet fits attach fitted means, Cook's distances, and row metadata", {
  testthat::skip_if_not_installed("DESeq2")
  counts <- .rsdeseq2_test_counts()
  col_data <- data.frame(
    condition = factor(rep(c("A", "B"), each = 4L), levels = c("A", "B")),
    row.names = colnames(counts)
  )
  rust_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
  reference_dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)

  rust_fit <- rsdeseq2::DESeq(rust_dds, fitType = "mean")
  reference_fit <- DESeq2::DESeq(reference_dds, fitType = "mean", quiet = TRUE)
  rust_state <- S4Vectors::metadata(rust_fit)$rsdeseq2RustFit
  row_metadata <- S4Vectors::mcols(rust_fit)
  reference_metadata <- S4Vectors::mcols(reference_fit)

  for (assay_name in c("mu", "cooks")) {
    testthat::expect_equal(
      SummarizedExperiment::assay(rust_fit, assay_name),
      SummarizedExperiment::assay(reference_fit, assay_name),
      tolerance = 1e-5
    )
  }
  for (column in c("baseMean", "baseVar", "allZero", "dispFit", "dispersion", "maxCooks")) {
    testthat::expect_equal(
      unname(row_metadata[[column]]),
      unname(reference_metadata[[column]]),
      tolerance = 1e-3
    )
  }
  testthat::expect_equal(row_metadata$MLE_conditionB, rust_state$beta[, "conditionB"])
  testthat::expect_equal(row_metadata$SE_conditionB, rust_state$betaSE[, "conditionB"])
  testthat::expect_equal(row_metadata$betaConv, rust_state$results$converged)
})
