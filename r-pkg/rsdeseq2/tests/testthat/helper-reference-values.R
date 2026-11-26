# Exact floating input transport shared by independent reference panels.
.row_reference_values <- function(values) vapply(values,function(value)
  readBin(as.raw(strtoi(substring(value,seq(1,15,2),seq(2,16,2)),16L)),"double",n=1,size=8,endian="big"),numeric(1))
.row_reference_bits <- function(values) lapply(unname(values),function(value)writeBin(value,raw(),size=8,endian="big"))

# Exact formula references were recorded with R 4.6.2 on macOS; other hosts
# have no recorded reference.
# Tests that compare expected values skip on other hosts; validation scripts
# that only take inputs from these files, and compare the Rust interfaces with
# each other, use them everywhere.
.formula_reference_path <- function(directory, name) {
  if (Sys.info()[["sysname"]] != "Darwin" && identical(Sys.getenv("TESTTHAT"), "true")) {
    testthat::skip("exact formula references are recorded for macOS R only")
  }
  file.path(directory, name)
}

# R's poly() reduces with the BLAS R links. rsdeseq2 reproduces the order of the
# recorded macOS R build exactly; on other hosts poly() results agree with live R
# to rounding only, so formulas containing poly() are compared numerically there.
.expect_same_frame <- function(actual, expected, formula, info = NULL) {
  exact <- Sys.info()[["sysname"]] == "Darwin" ||
    !grepl("poly(", paste(deparse(formula), collapse = ""), fixed = TRUE)
  if (exact) {
    testthat::expect_true(identical(actual, expected, num.eq = FALSE, single.NA = FALSE),
      info = info)
  } else {
    testthat::expect_equal(actual, expected, tolerance = 1e-12, info = info)
  }
}
