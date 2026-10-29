test_that("formula normalization is aligned and composed exactly once", {
  ct <- .rsdeseq2_test_counts()
  cd <- data.frame(condition = factor(rep(c("A", "B"), each = 4)),
                   exposure = seq(.8, 1.2, length.out = 8), row.names = colnames(ct))
  dds <- DESeq2::DESeqDataSetFromMatrix(ct, cd, ~ condition + offset(log(exposure)))
  sf <- setNames(seq(.7, 1.4, length.out = 8), colnames(ct))
  fit <- DESeq(dds, fitType = "mean", sizeFactors = sf, minReplicatesForReplace = Inf)
  again <- DESeq(fit, fitType = "mean", minReplicatesForReplace = Inf)
  nf <- DESeq2::normalizationFactors(fit)
  expect_equal(nf, matrix(rep(unname(sf) * cd$exposure, each = nrow(ct)),
                         nrow(ct), dimnames = dimnames(ct)), tolerance = 0)
  expect_equal(DESeq2::normalizationFactors(again), nf, tolerance = 0)
  expect_equal(as.data.frame(results(again)), as.data.frame(results(fit)), tolerance = 0)
  legacy <- fit
  S4Vectors::metadata(legacy)$rsdeseq2Normalization <- NULL
  expect_error(results(legacy), "offset provenance changed")
  restored <- unserialize(serialize(fit, NULL))
  expect_equal(as.data.frame(results(restored)), as.data.frame(results(fit)), tolerance = 0)
  base <- matrix(rep(unname(sf), each = nrow(ct)), nrow(ct), dimnames = dimnames(ct))
  nfFit <- DESeq(dds, fitType = "mean", normalizationFactors = base,
                 minReplicatesForReplace = Inf)
  permuted <- DESeq(dds, fitType = "mean", normalizationFactors = base[, 8:1],
                    minReplicatesForReplace = Inf)
  expect_equal(DESeq2::normalizationFactors(nfFit), DESeq2::normalizationFactors(permuted), tolerance = 0)
  expect_equal(as.data.frame(results(nfFit)), as.data.frame(results(permuted)), tolerance = 0)
  # Blindness changes trend/design fitting, not supplied exposure normalization.
  getFactors <- function(object) {
    rsdeseq2:::.rsdeseq2_prepare_normalization(
      object, ct, DESeq2::sizeFactors(object), DESeq2::normalizationFactors(object), log(cd$exposure))$normalizationFactors
  }
  expect_equal(getFactors(fit), nf, tolerance = 0)
  explicit <- DESeq2::DESeqDataSetFromMatrix(ct, cd, ~ condition)
  manual <- DESeq(explicit, fitType = "mean", normalizationFactors = nf,
                  minReplicatesForReplace = Inf)
  expect_equal(as.data.frame(results(fit)), as.data.frame(results(manual)), tolerance = 0)
  plainDesign <- stats::model.matrix(~ condition, cd)
  input <- list(counts = ct, design = plainDesign, sizeFactors = unname(sf),
    normalizationFactors = base * 2, weights = matrix(0.8, nrow(ct), ncol(ct)))
  cleared <- DESeq(input, fitType = "mean", normalizationFactors = NULL,
    sizeFactors = NULL, weights = NULL, minReplicatesForReplace = Inf)
  plain <- DESeq(ct, design = plainDesign, fitType = "mean", minReplicatesForReplace = Inf)
  expect_null(cleared$normalizationFactors)
  expect_null(cleared$weights)
  expect_equal(cleared$sizeFactors, plain$sizeFactors, tolerance = 0)
  expect_equal(results(cleared), results(plain), tolerance = 0)
  weighted <- DESeq(dds, fitType = "mean", sizeFactors = sf,
    weights = matrix(rep(c(1, 0.8), 4L), nrow(ct), ncol(ct), byrow = TRUE), minReplicatesForReplace = Inf)
  expect_s4_class(results(weighted), "DESeqResults")
  unweighted <- DESeq(weighted, fitType = "mean", weights = NULL, minReplicatesForReplace = Inf)
  expect_false("weights" %in% SummarizedExperiment::assayNames(unweighted))
  expect_equal(as.data.frame(results(unweighted)), as.data.frame(results(fit)), tolerance = 0)
})
