#!/usr/bin/env Rscript
# Black-box probe: how R orders strings under the ICU collation of a locale,
# compared with the root collation.
#
# Only R's own string comparison is consulted (`icuSetCollate`,
# `icuGetCollate`, `<`, `order`). For each locale the same strings are sorted:
# single letters of the Latin, Greek and Cyrillic blocks, all two-letter
# strings over an alphabet of plain and accented letters, all three-letter
# lowercase ASCII strings and random words (fixed seed). The rules of
# rcompat-collate's tailorings are read off the differences from the root
# order, and the sorted lists are the verification oracle.
#
# Usage: Rscript --vanilla scripts/provenance/probe_collation_tailorings.R [dir]
# Writes tailoring_<locale>.tsv (sorted strings with the relation to the
# next, as probe_collation_oracle.R writes them) and tailoring_locales.tsv
# into /tmp/rcollate by default.
stopifnot(getRversion() == "4.6.2", capabilities("ICU"))
args <- commandArgs(TRUE)
dir <- if (length(args)) args[[1]] else "/tmp/rcollate"
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
oracle <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
  "probe_collation_oracle.R")
locales <- c("root", "sv", "da", "nb", "nn", "fi", "is", "es", "ca", "gl", "pl", "cs", "sk", "hu", "hr", "sl", "bs", "ro",
  "lt", "lv", "et", "tr", "az", "ru", "uk", "be", "bg", "sr", "mk", "el", "de", "de@collation=phonebook",
  "de_DE@collation=phonebook", "fr_CA", "pt", "it", "nl", "en", "sq", "mt", "cy", "af", "fo", "kl", "se",
  "es@collation=traditional", "sr_Latn", "en_US_POSIX", "vi", "ha", "yo", "eo")
# Which collation ICU reports for each name.
rep <- do.call(rbind, lapply(locales, function(l) {
  icuSetCollate(locale = l)
  data.frame(requested = l, actual = icuGetCollate("actual"))
}))
icuSetCollate(locale = "root")
writeLines(c(sprintf("# %s on %s; ICU %s", R.version.string, R.version$platform, extSoftVersion()[["ICU"]]),
  "# requested\tactual", paste(rep$requested, rep$actual, sep = "\t")), file.path(dir, "tailoring_locales.tsv"))
set.seed(20261127)
latin <- c(0x41:0x5A, 0x61:0x7A)
single <- c(0x20:0x7E, 0xA0:0x24F, 0x370:0x3FF, 0x400:0x52F, 0x1E00:0x1EFF)
accented <- c(0xC0:0xFF, 0x100:0x17F, 0x1A0, 0x1A1, 0x1AF, 0x1B0, 0x218:0x21B)
alpha <- c(latin, accented, 0x6C + 0 * 1, 0x387)
cyr <- c(0x400:0x45F, 0x490, 0x491, 0x4E9, 0x4AF, 0x4BB)
greek <- c(0x386:0x3CE)
strings <- c(intToUtf8(single, multiple = TRUE),
  as.vector(outer(intToUtf8(alpha, multiple = TRUE), intToUtf8(alpha, multiple = TRUE), paste0)),
  as.vector(outer(outer(letters, letters, paste0), letters, paste0)),
  as.vector(outer(intToUtf8(cyr, multiple = TRUE), intToUtf8(cyr, multiple = TRUE), paste0)),
  as.vector(outer(intToUtf8(greek, multiple = TRUE), intToUtf8(greek, multiple = TRUE), paste0)),
  paste0(c("l", "L", "n", "N", "d", "D", "c", "C", "s", "S"), "·", "a"),
  vapply(1:20000, function(i) intToUtf8(sample(c(alpha, latin, latin, 0x300:0x30C, 0x323, 0x327, 0x328, 0x20, 0x2D, 0x30:0x39),
    sample(1:6, 1L), replace = TRUE)), ""),
  vapply(1:5000, function(i) intToUtf8(sample(c(cyr, greek, latin, 0x306, 0x308, 0x301), sample(1:5, 1L), replace = TRUE)), ""))
strings <- unique(strings[nzchar(strings) & !grepl("[\t\n]", strings)])
hex <- vapply(strings, function(s) paste(sprintf("%04X", utf8ToInt(s)), collapse = " "), "", USE.NAMES = FALSE)
q <- file.path(dir, "tailoring_input.tsv")
writeLines(hex, q)
for (l in locales) {
  out <- file.path(dir, sprintf("tailoring_%s.tsv", gsub("[@=]", "_", l)))
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(oracle), "sort", shQuote(l), shQuote(q), shQuote(out)))
  stopifnot(status == 0)
}
cat("strings", length(strings), "locales", length(locales), "\n")
