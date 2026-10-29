test_that("fitted selectors share Cook policy and expose raw Wald and LRT values", {
  counts <- .rsdeseq2_test_counts()
  counts[1, 8] <- 2000L
  cd <- data.frame(condition = factor(rep(c("A", "B"), each = 4)), row.names = colnames(counts))
  for (test in c("Wald", "LRT")) {
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, cd, ~ condition)
    fit <- DESeq(dds, fitType = "mean", test = test,
                 reduced = if (test == "LRT") ~ 1 else NULL,
                 minReplicatesForReplace = Inf)
    selected <- tail(resultsNames(fit), 1)
    get <- function(...) as.data.frame(results(fit, independentFiltering = FALSE, ...))
    default <- get()
    named <- get(name = selected)
    expect_equal(default$pvalue, named$pvalue, tolerance = 1e-9, info = test)
    expect_true(is.na(default$pvalue[1]), info = test)
    expect_equal(default$pvalue, get(contrast = c(0, 1))$pvalue, tolerance = 1e-9)
    expect_true(is.na(get(cooksCutoff = 0)$pvalue[1]))
    expect_error(results(fit, cooksCutoff = -1), "non-negative")
    unmasked <- get(cooksCutoff = FALSE)
    expect_true(is.finite(unmasked$pvalue[1]), info = test)
    expect_equal(unmasked$pvalue, get(name = selected, cooksCutoff = FALSE)$pvalue,
                 tolerance = 1e-9, info = test)
    expect_equal(unmasked$padj, p.adjust(unmasked$pvalue, "BH"), info = test)
    expect_error(results(fit, totallyUnknownArgument = TRUE), "totallyUnknownArgument")
    args <- list(object = dds, fitType = "mean", test = test,
                 minReplicatesForReplace = Inf, quiet = TRUE)
    if (test == "LRT") args$reduced <- ~ 1
    reference <- do.call(DESeq2::DESeq, args)
    referenceMasked <- DESeq2::results(reference, independentFiltering = FALSE)$pvalue
    expect_identical(is.na(default$pvalue), unname(is.na(referenceMasked)))
    rp <- DESeq2::results(reference, cooksCutoff = FALSE, independentFiltering = FALSE)$pvalue
    expect_identical(is.na(unmasked$pvalue), unname(is.na(rp)))
    expect_lt(max(abs(unmasked$pvalue - rp), na.rm = TRUE), 1e-6)
  }
})
