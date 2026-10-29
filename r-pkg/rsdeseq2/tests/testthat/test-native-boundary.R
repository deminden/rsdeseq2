test_that("registered native calls reject malformed buffers in a child process", {
  dll <- normalizePath(getLoadedDLLs()[["rsdeseq2"]][["path"]])
  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    sprintf("dyn.load(%s)", deparse(dll)),
    'ct <- matrix(as.double(1:12), 3, 4)',
    'dm <- cbind(1, c(0,0,1,1))',
    'call <- function(name, ...) .Call(name, ..., PACKAGE="rsdeseq2")',
    'bad <- list(',
    'function() call("rsdeseq2_normalized_counts", ct, 1, NULL),',
    'function() call("rsdeseq2_base_mean", ct, NULL, matrix(1, 2, 4)),',
    'function() call("rsdeseq2_base_metadata", ct, NULL, NULL, matrix(1, 3, 3)),',
    'function() call("rsdeseq2_estimate_size_factors", ct, c(1,1), NULL, 0L),',
    'function() call("rsdeseq2_estimate_size_factors", ct, NULL, c(NA_integer_), 0L),',
    'function() call("rsdeseq2_estimate_size_factors", ct, NULL, 1L, 2^32),',
    'function() call("rsdeseq2_apply_cooks_cutoff", rep(.5,3), 1, 1, NULL, NULL, FALSE),',
    'function() call("rsdeseq2_apply_cooks_cutoff", rep(.5,3), rep(1,3), 1, ct, matrix(1,3,3), TRUE),',
    'function() call("rsdeseq2_fit", ct, dm, NULL, NULL, NULL, NULL, 2L, 0L, 0L, 1, NA_integer_, 0L),',
    'function() call("rsdeseq2_transform", ct, dm, NULL, NULL, NULL, 2L, 1L, FALSE, 0L, 1L, 1, 1, -1L, NULL, NULL, NULL),',
    'function() call("rsdeseq2_lfc_shrink_normal_contrast", ct, dm, NULL, NULL, NULL, rep(.1,3), rep(1,3), rep(.1,3), matrix(0,3,2), c(0,1), c(0L,1L), c(0L,1L), c(0,1), 1L, 1L, 0),',
    'function() call("rsdeseq2_lfc_shrink_normal", ct, dm, 1, NULL, NULL, rep(.1,3), rep(1,3), rep(.1,3), 2L, 0L, 0)',
    ')',
    'for (i in seq_along(bad)) {',
    '  ans <- tryCatch({bad[[i]](); NULL}, error=identity)',
    '  if (!inherits(ans,"error")) stop(sprintf("malformed case %d succeeded", i))',
    '  stopifnot(grepl("required type and length|dimensions|indices|integer range", conditionMessage(ans)))',
    '}',
    'answer <- call("rsdeseq2_normalized_counts", ct, rep(2,4), NULL)',
    'stopifnot(identical(answer, ct / 2))',
    'set.seed(42); ct <- matrix(as.double(rnbinom(320, mu=50, size=5)),40,8)',
    'dm <- cbind(1, rep(c(0,1),each=4))',
    'fit <- call("rsdeseq2_fit",ct,dm,NULL,rep(1,8),NULL,NULL,2L,0L,0L,10,1L,0L)',
    'stopifnot(is.null(fit$replaceCounts),fit$nRefit==0L,identical(dim(fit$mu),dim(ct)))'
  ), script)
  output <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
                    stdout = TRUE, stderr = TRUE, timeout = 60)
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
})
