testthat::test_that("registered native bridge primitives are callable when package DLL is loaded", {
  testthat::skip_if_not("rsdeseq2" %in% names(getLoadedDLLs()))

  schema <- .Call("rsdeseq2_diagnostic_schema", PACKAGE = "rsdeseq2")
  testthat::expect_identical(schema, rsdeseq2DiagnosticSchemaRust())

  counts <- matrix(
    c(
      10, 20,
      20, 40,
      0, 5
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g0"), c("s1", "s2"))
  )
  logGeoMeans <- c(mean(log(c(10, 20))), mean(log(c(20, 40))), -Inf)
  sizeFactors <- .Call(
    "rsdeseq2_estimate_size_factors",
    counts,
    logGeoMeans,
    as.integer(seq_len(nrow(counts))),
    FALSE,
    PACKAGE = "rsdeseq2"
  )
  testthat::expect_equal(unname(sizeFactors), c(1 / sqrt(2), sqrt(2)), tolerance = 1e-12)

  normalized <- .Call(
    "rsdeseq2_normalized_counts",
    counts,
    c(1, 2),
    NULL,
    PACKAGE = "rsdeseq2"
  )
  testthat::expect_equal(normalized, normalizedCountsRust(counts, c(1, 2)))

  normalizationFactors <- matrix(
    c(
      1, 4,
      2, 8,
      1, 5
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = dimnames(counts)
  )
  normalizedByFactors <- .Call(
    "rsdeseq2_normalized_counts",
    counts,
    c(100, 100),
    normalizationFactors,
    PACKAGE = "rsdeseq2"
  )
  testthat::expect_equal(
    normalizedByFactors,
    normalizedCountsRust(counts, sizeFactors = c(100, 100), normalizationFactors = normalizationFactors)
  )

  baseMean <- .Call(
    "rsdeseq2_base_mean",
    counts,
    c(1, 2),
    NULL,
    PACKAGE = "rsdeseq2"
  )
  testthat::expect_equal(baseMean, baseMeanRust(counts, c(1, 2)))

  baseMeanByFactors <- .Call(
    "rsdeseq2_base_mean",
    counts,
    c(100, 100),
    normalizationFactors,
    PACKAGE = "rsdeseq2"
  )
  testthat::expect_equal(
    baseMeanByFactors,
    baseMeanRust(counts, sizeFactors = c(100, 100), normalizationFactors = normalizationFactors)
  )

  metadata <- .Call(
    "rsdeseq2_base_metadata",
    counts,
    c(1, 2),
    NULL,
    NULL,
    PACKAGE = "rsdeseq2"
  )
  testthat::expect_equal(metadata, baseMetadataRust(counts, c(1, 2)))

  weights <- matrix(
    c(
      1, 0.5,
      0, 1,
      1, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = dimnames(counts)
  )
  weightedMetadata <- .Call(
    "rsdeseq2_base_metadata",
    counts,
    c(100, 100),
    normalizationFactors,
    weights,
    PACKAGE = "rsdeseq2"
  )
  testthat::expect_equal(
    weightedMetadata,
    baseMetadataRust(
      counts,
      sizeFactors = c(100, 100),
      normalizationFactors = normalizationFactors,
      weights = weights
    )
  )

  cooksMasked <- .Call(
    "rsdeseq2_apply_cooks_cutoff",
    c(g1 = 0.01, g2 = 0.02, g3 = 0.5),
    c(0, 10, NA),
    5,
    NULL,
    NULL,
    FALSE,
    PACKAGE = "rsdeseq2"
  )
  testthat::expect_equal(cooksMasked$pvalue, c(g1 = 0.01, g2 = NA, g3 = 0.5))
  testthat::expect_equal(cooksMasked$cooksOutlier, c(FALSE, TRUE, NA))

  heuristicCounts <- matrix(c(1, 5, 6, 7, 9, 5, 6, 7), nrow = 2, byrow = TRUE)
  heuristicCooks <- matrix(c(10, 0.1, 0.2, 0.3, 10, 0.1, 0.2, 0.3), nrow = 2, byrow = TRUE)
  heuristicMasked <- .Call(
    "rsdeseq2_apply_cooks_cutoff",
    c(0.01, 0.02),
    c(10, 10),
    5,
    heuristicCounts,
    heuristicCooks,
    TRUE,
    PACKAGE = "rsdeseq2"
  )
  testthat::expect_equal(heuristicMasked$pvalue, c(0.01, NA))
  testthat::expect_equal(heuristicMasked$cooksOutlier, c(FALSE, TRUE))
})
