.rsdeseq2_test_counts <- function() {
  set.seed(20261014)
  n_genes <- 40L
  base_means <- exp(seq(log(10), log(500), length.out = n_genes))
  dispersions <- seq(0.03, 0.35, length.out = n_genes)
  fold_changes <- rep(c(1.5, 1, 2, 0.75), length.out = n_genes)
  counts <- t(vapply(seq_len(n_genes), function(gene) {
    means <- c(rep(base_means[[gene]], 4L), rep(base_means[[gene]] * fold_changes[[gene]], 4L))
    stats::rnbinom(8L, mu = means, size = 1 / dispersions[[gene]])
  }, numeric(8L)))
  rownames(counts) <- paste0("gene", seq_len(n_genes))
  colnames(counts) <- paste0("sample", seq_len(ncol(counts)))
  counts
}

.rsdeseq2_large_transform_fixture <- function() {
  set.seed(314159)
  n_genes <- 1200L
  n_samples <- 12L
  condition <- factor(rep(c("A", "B"), each = n_samples / 2L), levels = c("A", "B"))
  base_mean <- stats::rlnorm(n_genes, log(100), 1.35)
  dispersion <- 0.02 + 0.35 / base_mean
  log_fold_change <- stats::rnorm(n_genes, 0, 0.7)
  size_factor <- exp(seq(log(0.75), log(1.3), length.out = n_samples))
  means <- outer(base_mean, size_factor)
  means[, condition == "B"] <- means[, condition == "B"] * exp(log_fold_change)
  counts <- matrix(
    stats::rnbinom(
      n_genes * n_samples,
      mu = as.vector(means),
      size = rep(1 / dispersion, n_samples)
    ),
    nrow = n_genes
  )
  rownames(counts) <- paste0("gene", seq_len(n_genes))
  colnames(counts) <- paste0("sample", seq_len(n_samples))
  col_data <- data.frame(condition = condition, row.names = colnames(counts))
  list(counts = counts, colData = col_data)
}
