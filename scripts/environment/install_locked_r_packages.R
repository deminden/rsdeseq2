# Invoked by install_r_reference.py after source hashes and dependency ordering pass.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
if (as.character(getRversion()) != args[[3L]]) stop("Reference R version must be ", args[[3L]])
rows <- read.delim(args[[1L]], stringsAsFactors = FALSE, check.names = FALSE)
lib <- args[[2L]]
if (args[[4L]] == "install") {
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(normalizePath(lib), .libPaths()))
}
for (i in seq_len(nrow(rows))) {
  name <- rows$package[[i]]
  if (args[[4L]] == "install") {
    utils::install.packages(rows$archive[[i]], repos = NULL, type = "source", lib = lib)
    description <- suppressWarnings(utils::packageDescription(name, lib.loc = lib))
  } else description <- suppressWarnings(utils::packageDescription(name))
  if (!is.list(description) || !identical(description$Version, rows$version[[i]])) {
    stop("Reference package version differs: ", name, " expected ", rows$version[[i]])
  }
}
cat(nrow(rows), "locked R package versions verified\n")
print(sessionInfo())
