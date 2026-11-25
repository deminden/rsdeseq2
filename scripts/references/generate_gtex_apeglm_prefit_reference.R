#!/usr/bin/env Rscript
# Original apeglm only: compiled prefit (nbinomC) and refined posterior
# (nbinomCR) for the two genes of the gtex_compact golden input whose
# refined intercepts are the most sensitive to last-bit differences in the
# prefit. Counts and design come from the gtex_v8_real_wald fixture; size
# factors and dispersions are the values of the native parametric fit of
# that input (golden workflow wald_parametric, fields fit.size_factors and
# fit.dispersion), written here as IEEE-754 bit patterns. The prior matches
# the golden shrink_apeglm call: Cauchy scale 0.001 on all but the
# intercept, Normal(0, 15^2) on the intercept, coefficient 4 reported.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("usage: FRESH_OUTPUT_JSON")
stopifnot(as.character(getRversion()) == "4.6.2", as.character(packageVersion("apeglm")) == "1.34.0", !file.exists(args[[1L]]))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~", " ", script, fixed = TRUE)), "..", ".."))
decode <- function(h) vapply(h, function(x) readBin(as.raw(strtoi(substring(x, seq(1, 15, 2), seq(2, 16, 2)), 16L)), "double", n = 1, size = 8, endian = "big"), numeric(1), USE.NAMES = FALSE)
bits <- function(v) lapply(as.numeric(v), function(x) paste0("bits:", paste(format(writeBin(x, raw(), size = 8, endian = "big")), collapse = "")))
expected <- function(v) lapply(as.numeric(v), function(x) if (is.na(x) && !is.nan(x)) "NA" else sprintf("%.17g", x))
size_factors <- decode(strsplit(paste0(
  "3ff4d0f40e1f8571,3fe9888ac43f7e6e,3ff333c1e8574263,3ff15e1fbf36ab10,3ff0fd8f8c5c28c2,3ff062e3281c7539,",
  "3ff74c87ec630b31,3ff20fb0ef47bc27,3feb5ccb31ac174e,3fe82ac9bbd8e865,3ff5abffcee8ca40,3feec6a44ba4248e,",
  "3ff061571a479106,3fefaebd8e466036,3fe9afb4c49b960e,3fe48443babddb06,3ff1cdaecc4ef5ad,3fe9c5bd373ae63e,",
  "40005d2b1872dc2e,3fe78feeb92a0ece,3ff4c483ece044e4,3fee23963c6ca604,3fe461a4c0691c54,3fec8042ce3aa04d,",
  "3fe9054950599b77,3fe42148bbeb32ce,3fe4c719510dd7f2,4001e69c150b868a,3ff36c68a8f8c798,3fedbc589c276095,",
  "3fea0fa27f28271d,3ffd6bd49c7c0718,3ff6aaee60cb48c1,3feb4ae363d3a9cc,3fe5b24a7a1637eb,3febc95f6aaf7517,",
  "3fee85e5c86b34c9,400025519939386f,3febd7949ed78734,3fe200713de268bc,3fee577850ff34b8,3ff0c3fb49880067,",
  "3ff5dcdc3a1cfd48,3ff0c5cebfcb0ce3,4001cc946e57c043,3fe317ab6c3dadc7,3ff1827a07e53926,3ff2d61b08e02e07,",
  "3ff1e7fd5b4374d5,3fe83c7422ba0cd1,3ff6e0d41d90a975,3fdf6b65ed8c753a,3ff66edd571cd829,3fed38cfc234829c,",
  "3ff4e76699e9ace0,3feed4ba3371075a,3fe741527388875d,3fe91ffb420adbff,3ff642ba47768a68,3fec55125ef0e7f6,",
  "3ff2df410ddaad56"), ",")[[1L]])
genes <- c(392L, 188L)  # zero-based rows of the compact counts
dispersions <- decode(c("401324acc984ddd7", "40017114910e0b49"))
counts <- as.matrix(read.delim(file.path(root, "crates/rsdeseq2/tests/data/gtex_v8_real_wald/compact/counts.tsv"), row.names = 1, check.names = FALSE))
design <- as.matrix(read.delim(file.path(root, "crates/rsdeseq2/tests/data/gtex_v8_real_wald/design_matrix.tsv"), row.names = 1, check.names = FALSE))
stopifnot(identical(colnames(counts), rownames(design)), length(size_factors) == ncol(counts))
Y <- counts[genes + 1L, , drop = FALSE]
offset <- matrix(log(size_factors), nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
pc <- list(no.shrink = 1, prior.mean = 0, prior.scale = 0.001, prior.df = 1, prior.no.shrink.mean = 0, prior.no.shrink.scale = 15)
fit <- function(method) apeglm::apeglm(Y = Y, x = design, param = dispersions, coef = 4, offset = offset, threshold = 0.3,
                                      prior.control = pc, method = method)
prefit <- fit("nbinomC"); refined <- fit("nbinomCR")
dimnames(design) <- NULL
case <- list(name = "gtex_compact_sensitive_intercepts", genes = genes,
  counts = list(rows = nrow(Y), columns = ncol(Y), values = bits(Y)),
  design = list(rows = nrow(design), columns = ncol(design), values = bits(design)),
  dispersions = bits(dispersions), size_factors = bits(size_factors), offsets = NULL, weights = NULL, mle = NULL,
  coefficient = 3L, scale = bits(0.001)[[1L]], nuisance = list(0L), nuisance_scale = bits(15)[[1L]],
  threshold = bits(0.3)[[1L]], flip_sign = TRUE, map_only = FALSE,
  prefit_mode = expected(prefit$map), mode = expected(refined$map), sd = expected(refined$sd),
  fsr = expected(refined$fsr), sign_svalues = expected(refined$svalue))
jsonlite::write_json(list(source = list(R = as.character(getRversion()), apeglm = as.character(packageVersion("apeglm"))),
                          cases = list(case)), args[[1L]], auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
