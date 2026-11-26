#!/usr/bin/env Rscript
# Black-box oracles for every tailoring rcompat-collate supports: R's order
# of every code point (probe_collation_full.R) and of the expansion strings
# (expansion_input.tsv, written by collation_expansion_strings.py) under each
# locale, and the expansion strings under root.
#
# Only R's own string comparison is consulted, through the two probe scripts
# it runs. Usage:
#   python3 scripts/provenance/collation_expansion_strings.py
#   Rscript --vanilla scripts/provenance/probe_collation_supported.R [dir [cores]]
# Writes <locale>_single.tsv and expansion_<locale>.tsv into /tmp/rcollate by
# default (about 13 MB and 0.4 MB each; they are not committed, their SHA-256
# go to oracle_summary.tsv).
stopifnot(getRversion() == "4.6.2", capabilities("ICU"))
args <- commandArgs(TRUE)
dir <- if (length(args)) args[[1]] else "/tmp/rcollate"
cores <- if (length(args) > 1L) as.integer(args[[2]]) else 4L
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
full <- file.path(here, "probe_collation_full.R")
oracle <- file.path(here, "probe_collation_oracle.R")
input <- file.path(dir, "expansion_input.tsv")
stopifnot(file.exists(input))
locales <- c("sv", "fi", "nb", "nn", "es", "gl", "pl", "cs", "sk", "hu", "hr", "bs", "sl", "ro", "lt", "lv",
  "et", "tr", "az", "ru", "uk", "be", "bg", "sr", "mk", "de@collation=phonebook", "da", "is", "el", "fr_CA")
rscript <- file.path(R.home("bin"), "Rscript")
run <- function(l) {
  name <- gsub("[@=]", "_", l)
  status <- system2(rscript, c("--vanilla", shQuote(oracle), "sort", shQuote(l), shQuote(input),
    shQuote(file.path(dir, sprintf("expansion_%s.tsv", name)))))
  if (l != "root") {
    status <- status + system2(rscript, c("--vanilla", shQuote(full), shQuote(dir), shQuote(l)), stdout = FALSE)
  }
  status
}
status <- unlist(parallel::mclapply(c("root", locales), run, mc.cores = cores))
stopifnot(all(status == 0))
cat("locales", length(locales), "\n")
