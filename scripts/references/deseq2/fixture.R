counts <- matrix(
  c(
    10L, 12L, 20L, 24L,
    0L, 0L, 5L, 7L,
    100L, 80L, 90L, 120L,
    3L, 6L, 9L, 12L
  ),
  nrow = 4,
  byrow = TRUE
)
rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
colnames(counts) <- paste0("sample", seq_len(ncol(counts)))

col_data <- S4Vectors::DataFrame(
  condition = factor(c("A", "A", "B", "B")),
  row.names = colnames(counts)
)

full_design <- stats::model.matrix(~ condition, data = as.data.frame(col_data))
reduced_design <- stats::model.matrix(~ 1, data = as.data.frame(col_data))
rownames(full_design) <- rownames(col_data)
rownames(reduced_design) <- rownames(col_data)
fixed_dispersions <- c(0.10, 0.15, 0.05, 0.20)
names(fixed_dispersions) <- rownames(counts)
bioconductor_version <- if (requireNamespace("BiocManager", quietly = TRUE)) {
  as.character(BiocManager::version())
} else {
  NA_character_
}

write_matrix_tsv(counts, "gene", file.path(data_dir, "counts.tsv"))
write_tsv(
  data.frame(
    sample = rownames(col_data),
    condition = as.character(col_data$condition),
    check.names = FALSE
  ),
  file.path(data_dir, "col_data.tsv")
)
write_matrix_tsv(full_design, "sample", file.path(data_dir, "design_full.tsv"))
write_matrix_tsv(reduced_design, "sample", file.path(data_dir, "design_reduced.tsv"))
write_tsv(
  data.frame(gene = names(fixed_dispersions), dispersion = fixed_dispersions),
  file.path(data_dir, "fixed_dispersions.tsv")
)

metadata <- data.frame(
  key = c(
    "DESeq2_version",
    "Bioconductor_version",
    "R_version",
    "platform",
    "reference_fixture",
    "fixed_glm_reference_mode",
    "beta_prior_reference_mode",
    "normalization_factors_reference_mode",
    "native_nf_dispersion_reference_mode",
    "fixed_glm_force_optim_reference_mode",
    "weighted_reference_mode",
    "native_glm_mu_cr_reference_mode",
    "native_weighted_glm_mu_cr_reference_mode",
    "native_weighted_glm_mu_reference_mode",
    "native_weighted_glm_mu_mean_cr_map_reference_mode",
    "native_glm_mu_mean_reference_mode",
    "native_glm_mu_mean_cr_map_reference_mode",
    "native_glm_mu_local_reference_mode",
    "native_glm_mu_local_cr_map_reference_mode",
    "native_weighted_glm_mu_local_reference_mode",
    "native_weighted_glm_mu_local_cr_map_reference_mode",
    "dispersion_trend_reference_mode",
    "dispersion_prior_reference_mode",
    "map_dispersion_reference_mode"
  ),
  value = c(
    as.character(utils::packageVersion("DESeq2")),
    bioconductor_version,
    paste(R.version$major, R.version$minor, sep = "."),
    R.version$platform,
    "tiny_two_group_four_gene",
    "DESeq2:::fitNbinomGLMs with supplied dispersions, default 1e-6 beta ridge, useQR=FALSE, useOptim=FALSE",
    "DESeq2:::estimateBetaPriorVar from supplied-dispersion MLE betas with weighted and quantile methods",
    "counts(dds, normalized=TRUE) with normalizationFactors(dds), which preempt sizeFactors(dds)",
    "DESeq2 roughDispEstimate, momentsDispEstimate, and estimateDispersionsGeneEst stored mu with normalizationFactors(dds)",
    "DESeq2:::fitNbinomGLMs with supplied dispersions, default 1e-6 beta ridge, useQR=FALSE, useOptim=TRUE, forceOptim=TRUE",
    "getBaseMeansAndVariances with raw weights, getAndCheckWeights row-normalized weights, and fitNbinomGLMs with supplied dispersions",
    "estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=TRUE) without observation weights for Cox-Reid gene-wise dispersion references",
    "estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=TRUE) with observation weights for weighted Cox-Reid gene-wise dispersion references",
    "estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=FALSE), estimateDispersionsFit(fitType='mean'), estimateDispersionsMAP(useCR=FALSE), and full/reduced fitNbinomGLMs(useQR=FALSE,useOptim=FALSE) with observation weights",
    "estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=TRUE), estimateDispersionsFit(fitType='mean'), and estimateDispersionsMAP(useCR=TRUE) with observation weights",
    "estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=FALSE), estimateDispersionsFit(fitType='mean'), and estimateDispersionsMAP(useCR=FALSE) without observation weights",
    "estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=TRUE), estimateDispersionsFit(fitType='mean'), and estimateDispersionsMAP(useCR=TRUE) without observation weights",
    "estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=FALSE), estimateDispersionsFit(fitType='local'), estimateDispersionsMAP(useCR=FALSE), and full/reduced fitNbinomGLMs(useQR=FALSE,useOptim=FALSE) without observation weights",
    "estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=TRUE), estimateDispersionsFit(fitType='local'), and estimateDispersionsMAP(useCR=TRUE) without observation weights",
    "estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=FALSE), estimateDispersionsFit(fitType='local'), estimateDispersionsMAP(useCR=FALSE), and full/reduced fitNbinomGLMs(useQR=FALSE,useOptim=FALSE) with observation weights",
    "estimateDispersionsGeneEst(linearMu=FALSE,niter=2,useCR=TRUE), estimateDispersionsFit(fitType='local'), and estimateDispersionsMAP(useCR=TRUE) with observation weights",
    "parametricDispersionFit, estimateDispersionsFit(fitType='mean'), and localDispersionFit fixtures including prediction and threshold edge cases",
    "estimateDispersionsPriorVar high-df and low-df fixtures; low-df uses DESeq2's seeded Monte Carlo/R loess branch",
    "DESeq2:::fitDisp MAP dispersion fixture with explicit prior and Cox-Reid settings"
  )
)
write_tsv(metadata, file.path(data_dir, "metadata.tsv"))
