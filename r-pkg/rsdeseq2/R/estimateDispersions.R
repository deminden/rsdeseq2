estimateDispersionsRust <- function(counts,
                                    design,
                                    fitType = c("parametric", "local", "mean"),
                                    nproc = 1L) {
  fit <- DESeq(
    object = list(counts = counts, design = design),
    fitType = fitType,
    test = "Wald",
    nproc = nproc
  )
  dispersion <- fit$results$dispersion
  names(dispersion) <- rownames(fit$results)
  dispersion
}
