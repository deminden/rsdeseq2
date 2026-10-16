# Small explicit inputs for wrapper dispatch tests; live DESeq2 oracles live in
# test-fit-integration.R and must remain independent of these fixtures.
.rsdeseq2_primitive_fixture <- function(
    values = c(1, 2, 3, 2, 4, 7),
    coefficients = c("Intercept", "condition_B_vs_A", "condition_C_vs_A"),
    variances = c(0, 0.25, 1)) {
  beta <- matrix(values, nrow = 2, byrow = TRUE,
                 dimnames = list(c("g1", "g2"), coefficients))
  covariance <- array(0, dim = c(2, ncol(beta), ncol(beta)))
  for (i in seq_along(variances)) covariance[, i, i] <- variances[[i]]
  list(baseMean = c(g1 = 10, g2 = 20), beta = beta, betaCovariance = covariance)
}

.rsdeseq2_count_fixture <- function(values) {
  counts <- matrix(values, nrow = 2, byrow = TRUE)
  dimnames(counts) <- list(c("g1", "g2"), paste0("s", seq_len(ncol(counts))))
  counts
}

.rsdeseq2_mcols_fixture <- function() {
  data.frame(
    baseMean = c(10, 20), Intercept = c(1, 1), condition_B_vs_A = c(2, -2),
    SE_Intercept = c(0.1, 0.1), SE_condition_B_vs_A = c(0.5, 0.5),
    row.names = c("g1", "g2"), check.names = FALSE
  )
}

.rsdeseq2_s4_fixture <- function(class, slots, fields) {
  if (!methods::isClass(class)) methods::setClass(class, slots = slots)
  do.call(methods::new, c(list(Class = class), fields))
}

# Explicit expected effects/SEs and supplied LRT statistics are the oracle.
# Wald tails and BH are computed independently with base R's statistics API.
.rsdeseq2_expect_result <- function(out, lfc, se, stat = lfc / se,
                                   pvalue = 2 * stats::pnorm(abs(stat), lower.tail = FALSE),
                                   info = NULL) {
  expected <- data.frame(
    baseMean = c(10, 20), log2FoldChange = lfc, lfcSE = se,
    stat = stat, pvalue = pvalue, padj = stats::p.adjust(pvalue, "BH")
  )
  testthat::expect_identical(names(out), c(names(expected), "filtered"), info = info)
  testthat::expect_identical(rownames(out), c("g1", "g2"), info = info)
  testthat::expect_true(all(is.finite(as.matrix(out[names(expected)]))), info = info)
  testthat::expect_identical(out$filtered, c(NA, NA), info = info)
  for (field in names(expected)) {
    testthat::expect_true(max(abs(out[[field]] - expected[[field]])) <= 1e-12, info = paste(info, field))
  }
}
