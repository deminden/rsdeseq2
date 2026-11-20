# Exact floating input transport shared by independent reference panels.
.row_reference_values <- function(values) vapply(values,function(value)
  readBin(as.raw(strtoi(substring(value,seq(1,15,2),seq(2,16,2)),16L)),"double",n=1,size=8,endian="big"),numeric(1))
.row_reference_bits <- function(values) lapply(unname(values),function(value)writeBin(value,raw(),size=8,endian="big"))

# Formula expectations follow the R host's reductions and BLAS schedule.
.formula_reference_path <- function(directory, name) {
  if (Sys.info()[["sysname"]] == "Linux") {
    machine <- Sys.info()[["machine"]]
    if (machine %in% c("x86_64", "aarch64"))
      name <- sub("\\.json$", paste0("_linux_", machine, ".json"), name)
  }
  file.path(directory, name)
}
