test_that("fitted views map gene order and reject mutations of fitted inputs", {
  ct <- .rsdeseq2_test_counts()
  cd <- data.frame(condition = factor(rep(c("A", "B"), each = 4)), row.names = colnames(ct))
  fit <- DESeq(DESeq2::DESeqDataSetFromMatrix(ct, cd, ~ condition),
               fitType = "mean", minReplicatesForReplace = Inf)
  full <- as.data.frame(results(fit, independentFiltering = FALSE))
  subset <- as.data.frame(results(fit[c(3, 1), ], independentFiltering = FALSE))
  expect_identical(rownames(subset), rownames(ct)[c(3, 1)])
  expect_equal(subset$pvalue, full$pvalue[c(3, 1)], tolerance = 0)
  expect_equal(subset$padj, p.adjust(subset$pvalue, "BH"), tolerance = 0)
  reversed <- as.data.frame(results(fit[nrow(fit):1L, ], independentFiltering = FALSE))
  expect_equal(reversed, full[nrow(full):1L, ], tolerance = 0)
  alteredSubset <- fit[c(3, 1), ]
  SummarizedExperiment::assay(alteredSubset, "counts")[1, 1] <- ct[3, 1] + 1L
  expect_error(results(alteredSubset), "counts changed")
  changed <- fit
  SummarizedExperiment::assay(changed, "counts")[1, 1] <- ct[1, 1] + 1L
  consumers <- list(results, resultsNames,
    function(x) lfcShrink(x, coef = 2, type = "normal"),
    function(x) lfcShrink(x, coef = 2, type = "ashr"),
    function(x) lfcShrink(x, coef = 2, type = "apeglm"))
  for (consumer in consumers) expect_error(consumer(changed), "counts changed")
  expect_error(results(fit[, 8:1]), "samples changed")
  changed <- fit
  DESeq2::sizeFactors(changed) <- DESeq2::sizeFactors(changed) * 2
  expect_error(results(changed), "size factors changed")
  changed <- fit
  SummarizedExperiment::colData(changed)$condition <- rev(cd$condition)
  expect_error(results(changed), "design changed")
  for (field in c("maxCooks", "betaConv")) {
    changed <- fit
    values <- S4Vectors::mcols(changed)[[field]]
    values[[1L]] <- if (field == "betaConv") !isTRUE(values[[1L]]) else 100
    S4Vectors::mcols(changed)[[field]] <- values
    expect_error(results(changed), paste(field, "changed"))
  }
  reference <- DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(ct, cd, ~ condition),
    fitType = "mean", minReplicatesForReplace = Inf, quiet = TRUE)
  rows <- seq.int(35L, 1L, by = -2L)
  selectedFit <- fit[rows, ]; selectedReference <- reference[rows, ]
  for (type in c("normal", "ashr", "apeglm")) {
    actual <- suppressMessages(lfcShrink(selectedFit, coef = 2, type = type))
    expected <- suppressMessages(DESeq2::lfcShrink(selectedReference, coef = 2, type = type, quiet = TRUE))
    expect_identical(rownames(actual), rownames(ct)[rows], info = type)
    expect_equal(as.numeric(actual$log2FoldChange), as.numeric(expected$log2FoldChange), tolerance = 1e-5, info = type)
    expect_equal(as.numeric(actual$lfcSE), as.numeric(expected$lfcSE), tolerance = 1e-5, info = type)
    expect_equal(actual$padj, p.adjust(actual$pvalue, "BH"), info = type)
  }
})

test_that("saved fit origin distinguishes formula semantics from identical numeric designs", {
  ct <- .rsdeseq2_test_counts()
  cd <- data.frame(condition=factor(rep(c("A","B"),each=4)),nuisance=seq_len(8),row.names=colnames(ct))
  fit <- DESeq(DESeq2::DESeqDataSetFromMatrix(ct,cd,~condition),fitType="mean",minReplicatesForReplace=Inf)
  state <- S4Vectors::metadata(fit)$rsdeseq2RustFit
  expect_identical(state$schemaVersion,2L)
  expect_identical(state$designOrigin,list(kind="formula",variables="condition"))
  path <- tempfile(fileext=".rds");on.exit(unlink(path))
  saveRDS(fit,path)
  expected <- as.data.frame(results(fit,independentFiltering=FALSE))
  expect_identical(as.data.frame(results(readRDS(path),independentFiltering=FALSE)),expected)
  changed <- fit
  DESeq2::design(changed) <- ~condition+nuisance-nuisance
  expect_identical(unname(model.matrix(DESeq2::design(changed),cd)),unname(model.matrix(~condition,cd)))
  expect_error(results(changed),"formula origin changed")
  for(bad in list(NULL,list(kind="legacy",variables=character()),
      list(kind="formula",variables=c("condition","condition")),
      list(kind="matrix",variables="condition"))) {
    changed <- fit; metadata <- S4Vectors::metadata(changed)
    metadata$rsdeseq2RustFit$designOrigin <- bad;S4Vectors::metadata(changed) <- metadata
    expect_error(results(changed),"origin")
  }
  legacy <- fit; metadata <- S4Vectors::metadata(legacy)
  metadata$rsdeseq2RustFit$designOrigin <- NULL
  metadata$rsdeseq2RustFit$schemaVersion <- 1L
  S4Vectors::metadata(legacy) <- metadata
  expect_identical(as.data.frame(results(legacy,independentFiltering=FALSE)),expected)
  expect_null(S4Vectors::metadata(legacy)$rsdeseq2RustFit$designOrigin)
  metadata$rsdeseq2RustFit$schemaVersion <- 3L
  S4Vectors::metadata(legacy) <- metadata
  expect_error(results(legacy),"Unsupported fitted schema")
  matrixFit <- DESeq(ct,design=model.matrix(~condition,cd),fitType="mean",minReplicatesForReplace=Inf)
  expect_identical(matrixFit$designOrigin,list(kind="matrix",variables=character()))
})
