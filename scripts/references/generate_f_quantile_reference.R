#!/usr/bin/env Rscript
# Independent qf(0.99, df1, df2) references for the Cook's distance cutoff.
# Runs base R only; never loads rsdeseq2.
stopifnot(getRversion() == "4.6.2")
args <- commandArgs(TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~", " ", script, fixed = TRUE)), "..", ".."))
path <- if (length(args)) args[[1L]] else
  file.path(root, "crates/rcompat-numeric/tests/data/f_quantile_r462.tsv")
if (file.exists(path)) stop("Use a fresh path; references are never overwritten")
hex <- function(v) paste(format(writeBin(as.double(v), raw(), size = 8, endian = "big")), collapse = "")
df1 <- c(1:12, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200)
df2 <- c(1:12, 15, 20, 30, 50, 100, 200, 500, 1000, 2000, 5000, 1e4, 2e4, 5e4, 1e5)
grid <- expand.grid(df1 = df1, df2 = df2)
grid$quantile <- vapply(seq_len(nrow(grid)), function(i)
  hex(qf(0.99, grid$df1[[i]], grid$df2[[i]])), "")
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
writeLines(c(paste("#", R.version.string, "on", R.version$platform),
             "df1\tdf2\tquantile_bits"),
           path)
write.table(grid, path, sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = FALSE, append = TRUE)
