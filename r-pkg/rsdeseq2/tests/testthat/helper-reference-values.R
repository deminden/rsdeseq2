# Exact floating input transport shared by independent reference panels.
.row_reference_values <- function(values) vapply(values,function(value)
  readBin(as.raw(strtoi(substring(value,seq(1,15,2),seq(2,16,2)),16L)),"double",n=1,size=8,endian="big"),numeric(1))
.row_reference_bits <- function(values) lapply(unname(values),function(value)writeBin(value,raw(),size=8,endian="big"))

# Linux x86 uses extended scalar accumulation and a fixed reference BLAS kernel.
.formula_reference_path <- function(directory, name) {
  if (Sys.info()[["sysname"]] == "Linux" && Sys.info()[["machine"]] == "x86_64")
    name <- sub("\\.json$", "_linux_x86_64.json", name)
  file.path(directory, name)
}
