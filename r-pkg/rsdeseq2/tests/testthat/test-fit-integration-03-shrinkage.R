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
  .rsdeseq2_expect_shrinkage_parity(
    rust_shrunk$log2FoldChange, reference_shrunk$log2FoldChange,
    tolerance = 1e-3
  )
  .rsdeseq2_expect_shrinkage_parity(rust_shrunk$lfcSE, reference_shrunk$lfcSE, tolerance = 1e-3)
  .rsdeseq2_expect_shrinkage_parity(rust_svalues$svalue, reference_svalues$svalue, tolerance = 0.02)
  .rsdeseq2_expect_shrinkage_parity(
    rust_thresholded$svalue, reference_thresholded$svalue, tolerance = 0.02
  )
  .rsdeseq2_expect_shrinkage_parity(
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
  .rsdeseq2_expect_shrinkage_parity(
    unname(rust_shrunk$log2FoldChange), unname(reference_shrunk$log2FoldChange),
    tolerance = 0.01
  )
  .rsdeseq2_expect_shrinkage_parity(
    unname(rust_shrunk$lfcSE), unname(reference_shrunk$lfcSE), tolerance = 0.01
  )
  .rsdeseq2_expect_shrinkage_parity(rust_svalues$svalue, reference_svalues$svalue, tolerance = 0.02)
  .rsdeseq2_expect_shrinkage_parity(
    unname(rust_thresholded$log2FoldChange),
    unname(reference_thresholded$log2FoldChange),
    tolerance = 0.01
  )
  .rsdeseq2_expect_shrinkage_parity(
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
