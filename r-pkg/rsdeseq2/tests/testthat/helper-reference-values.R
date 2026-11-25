# Exact floating input transport shared by independent reference panels.
.row_reference_values <- function(values) vapply(values,function(value)
  readBin(as.raw(strtoi(substring(value,seq(1,15,2),seq(2,16,2)),16L)),"double",n=1,size=8,endian="big"),numeric(1))
.row_reference_bits <- function(values) lapply(unname(values),function(value)writeBin(value,raw(),size=8,endian="big"))

# Exact formula references were recorded with R 4.6.2 on macOS; other hosts
# have no recorded reference.
.formula_reference_path <- function(directory, name) {
  if (Sys.info()[["sysname"]] != "Darwin") {
    message <- "exact formula references are recorded for macOS R only"
    if (identical(Sys.getenv("TESTTHAT"), "true")) testthat::skip(message)
    stop(message, call. = FALSE)
  }
  file.path(directory, name)
}
