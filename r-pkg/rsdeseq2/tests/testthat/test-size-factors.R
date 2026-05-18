testthat::test_that("R wrapper estimates ratio size factors on log-ratio scale", {
  counts <- matrix(
    c(
      10L, 20L,
      20L, 40L,
      0L, 5L
    ),
    nrow = 3,
    byrow = TRUE
  )

  testthat::expect_equal(
    unname(estimateSizeFactorsRust(counts, method = "ratio")),
    c(1 / sqrt(2), sqrt(2)),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    estimateSizeFactorsRust(counts, method = "ratio", native = TRUE),
    estimateSizeFactorsRust(counts, method = "ratio"),
    tolerance = 1e-12
  )
})

testthat::test_that("R wrapper estimates poscounts size factors with zeros", {
  counts <- matrix(
    c(
      0L, 4L, 16L,
      9L, 0L, 36L,
      0L, 0L, 0L
    ),
    nrow = 3,
    byrow = TRUE
  )
  sf <- estimateSizeFactorsRust(counts, method = "poscounts")
  nativeSf <- estimateSizeFactorsRust(counts, method = "poscounts", native = TRUE)

  testthat::expect_length(sf, 3)
  testthat::expect_true(all(is.finite(sf)))
  testthat::expect_true(all(sf > 0))
  testthat::expect_equal(nativeSf, sf, tolerance = 1e-12)
})

testthat::test_that("R wrapper supports supplied geoMeans and control genes", {
  counts <- matrix(
    c(
      10L, 20L,
      20L, 40L,
      100L, 100L
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("s1", "s2"))
  )
  sf <- estimateSizeFactorsRust(
    counts,
    geoMeans = c(sqrt(200), sqrt(800), 100),
    controlGenes = c("a", "b")
  )
  nativeSf <- estimateSizeFactorsRust(
    counts,
    geoMeans = c(sqrt(200), sqrt(800), 100),
    controlGenes = c("a", "b"),
    native = TRUE
  )

  testthat::expect_equal(names(sf), c("s1", "s2"))
  testthat::expect_equal(exp(mean(log(sf))), 1, tolerance = 1e-12)
  testthat::expect_equal(unname(sf), c(1 / sqrt(2), sqrt(2)), tolerance = 1e-12)
  testthat::expect_equal(nativeSf, sf, tolerance = 1e-12)
})

testthat::test_that("R wrapper matches DESeq2 size-factor validation cases", {
  # Behavioral port of the validation block in
  # external/DESeq2/tests/testthat/test_size_factor.R.
  counts <- matrix(as.integer(1:16), ncol = 4)

  testthat::expect_error(estimateSizeFactorsRust(counts, geoMeans = 1:5))
  testthat::expect_error(estimateSizeFactorsRust(counts, geoMeans = rep(0, 4)))
  testthat::expect_error(estimateSizeFactorsRust(counts, controlGenes = "foo"))
  testthat::expect_silent(estimateSizeFactorsRust(counts, geoMeans = 1:4))
  testthat::expect_silent(estimateSizeFactorsRust(counts, controlGenes = 1:2))
})

testthat::test_that("R wrapper normalized counts and base means are available", {
  counts <- matrix(
    c(
      10L, 20L,
      30L, 60L
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("s1", "s2"))
  )
  normalized <- normalizedCountsRust(counts, c(1, 2))
  nativeNormalized <- normalizedCountsRust(counts, c(1, 2), native = TRUE)

  testthat::expect_equal(dimnames(normalized), dimnames(counts))
  testthat::expect_equal(unname(normalized[, "s2"]), c(10, 30))
  testthat::expect_equal(nativeNormalized, normalized)
  testthat::expect_equal(baseMeanRust(counts, c(1, 2)), c(g1 = 10, g2 = 30))
  testthat::expect_equal(baseMeanRust(counts, c(1, 2), native = TRUE), c(g1 = 10, g2 = 30))
})

testthat::test_that("R wrapper normalization factors preempt size factors", {
  counts <- matrix(
    c(
      10L, 20L,
      30L, 60L
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("s1", "s2"))
  )
  normalizationFactors <- matrix(
    c(
      1, 4,
      3, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = dimnames(counts)
  )

  normalized <- normalizedCountsRust(
    counts,
    sizeFactors = c(100, 100),
    normalizationFactors = normalizationFactors
  )
  nativeNormalized <- normalizedCountsRust(
    counts,
    sizeFactors = c(100, 100),
    normalizationFactors = normalizationFactors,
    native = TRUE
  )

  testthat::expect_equal(dimnames(normalized), dimnames(counts))
  testthat::expect_equal(as.vector(normalized), c(10, 10, 5, 30))
  testthat::expect_equal(nativeNormalized, normalized)
  testthat::expect_equal(
    baseMeanRust(counts, normalizationFactors = normalizationFactors),
    c(g1 = 7.5, g2 = 20)
  )
  testthat::expect_equal(
    baseMeanRust(counts, normalizationFactors = normalizationFactors, native = TRUE),
    c(g1 = 7.5, g2 = 20)
  )
})

testthat::test_that("R wrapper base metadata includes base variance and all-zero state", {
  counts <- matrix(
    c(
      10L, 20L, 0L,
      30L, 60L, 0L,
      0L, 0L, 0L
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g0"), c("s1", "s2", "s3"))
  )
  normalized <- sweep(counts, 2L, c(1, 2, 1), "/")
  meta <- baseMetadataRust(counts, sizeFactors = c(1, 2, 1))
  nativeMeta <- baseMetadataRust(counts, sizeFactors = c(1, 2, 1), native = TRUE)

  testthat::expect_equal(rownames(meta), rownames(counts))
  testthat::expect_equal(meta$baseMean, unname(rowMeans(normalized)), tolerance = 1e-12)
  testthat::expect_equal(meta$baseVar, as.numeric(apply(normalized, 1L, stats::var)), tolerance = 1e-12)
  testthat::expect_equal(meta$allZero, c(FALSE, FALSE, TRUE))
  testthat::expect_equal(nativeMeta, meta)
})

testthat::test_that("R wrapper base metadata supports normalization factors and weights", {
  counts <- matrix(
    c(
      10L, 20L,
      30L, 60L
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("s1", "s2"))
  )
  normalizationFactors <- matrix(
    c(
      1, 4,
      3, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = dimnames(counts)
  )
  weights <- matrix(
    c(
      1, 0.5,
      0, 1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = dimnames(counts)
  )
  metadataCounts <- counts / normalizationFactors * weights
  meta <- baseMetadataRust(
    counts,
    normalizationFactors = normalizationFactors,
    weights = weights
  )

  testthat::expect_equal(meta$baseMean, unname(rowMeans(metadataCounts)), tolerance = 1e-12)
  testthat::expect_equal(meta$baseVar, as.numeric(apply(metadataCounts, 1L, stats::var)), tolerance = 1e-12)
})

testthat::test_that("R wrapper validates primitive inputs", {
  testthat::expect_error(estimateSizeFactorsRust(matrix(c(0L, 0L), nrow = 1)))
  testthat::expect_error(estimateSizeFactorsRust(matrix(-1L, nrow = 1)))
  testthat::expect_error(estimateSizeFactorsRust(matrix(1L, nrow = 1), geoMeans = NA_real_))
  testthat::expect_error(estimateSizeFactorsRust(matrix(1L, nrow = 1), controlGenes = NA))
  testthat::expect_error(normalizedCountsRust(matrix(1L, nrow = 1), 0))
  testthat::expect_error(normalizedCountsRust(matrix(1L, nrow = 1)))
  testthat::expect_error(normalizedCountsRust(matrix(1L, nrow = 1), normalizationFactors = matrix(0, nrow = 1)))
  testthat::expect_error(normalizedCountsRust(matrix(1L, nrow = 1), normalizationFactors = matrix(1, nrow = 1, ncol = 2)))
  testthat::expect_error(baseMetadataRust(matrix(1L, nrow = 1), sizeFactors = 1, weights = matrix(Inf, nrow = 1)))
  testthat::expect_error(baseMetadataRust(matrix(1L, nrow = 1), sizeFactors = 1, weights = matrix(1, nrow = 1, ncol = 2)))
})
