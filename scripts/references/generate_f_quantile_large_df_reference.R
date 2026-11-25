#!/usr/bin/env Rscript
# Independent qf(0.99, df1, df2) references for very large degrees of freedom
# (docs/provenance/specs/f_quantile_large_df.md): one of df1, df2 lies between
# 400001 and 1e7, the other between 1 and 200, plus the pairs from
# scripts/provenance/probe_qf_large_df.R where both are large.
# Runs base R only; never loads rsdeseq2.
stopifnot(getRversion() == "4.6.2")
args <- commandArgs(TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~", " ", script, fixed = TRUE)), "..", ".."))
path <- if (length(args)) args[[1L]] else
  file.path(root, "crates/rcompat-numeric/tests/data/f_quantile_large_df_r462.tsv")
if (file.exists(path)) stop("Use a fresh path; references are never overwritten")
hex <- function(v) paste(format(writeBin(as.double(v), raw(), size = 8, endian = "big")), collapse = "")
small <- c(1, 2, 3, 5, 10, 20, 50, 100, 200)
large <- c(400001, 400002, 5e5, 1e6, 3e6, 1e7)
grid <- rbind(expand.grid(df1 = small, df2 = large),
              setNames(expand.grid(df2 = small, df1 = large)[, c("df1", "df2")],
                       c("df1", "df2")),
              data.frame(df1 = c(5e5, 6e5, 5e5, 4e5), df2 = c(6e5, 5e5, 5e5, 400001)))
grid$quantile <- vapply(seq_len(nrow(grid)), function(i)
  hex(qf(0.99, grid$df1[[i]], grid$df2[[i]])), "")
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
writeLines(c(paste("#", R.version.string, "on", R.version$platform),
             "df1\tdf2\tquantile_bits"),
           path)
write.table(format(grid, scientific = FALSE, trim = TRUE), path, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
