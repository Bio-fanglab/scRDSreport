test_that("QC columns and thresholds are detected", {
  meta <- data.frame(nCount_RNA = c(0, 100, 1000), nFeature_RNA = c(0, 100, 500))
  result <- scRDSreport:::.qc_keep(meta, min_features = 200, min_counts = 500)
  expect_equal(result$columns$features, "nFeature_RNA")
  expect_equal(result$columns$counts, "nCount_RNA")
  expect_equal(result$keep, c(FALSE, FALSE, TRUE))
})

test_that("QC filtering tolerates missing count columns", {
  meta <- data.frame(nFeature_RNA = c(10, 300))
  result <- scRDSreport:::.qc_keep(meta, min_features = 200)
  expect_equal(result$keep, c(FALSE, TRUE))
  expect_null(result$columns$counts)
})

test_that("low-expression feature filtering validates its threshold", {
  expect_error(
    scRDSreport:::.prefilter_low_expression_features(
      matrix(1, 10, 10), mode = "always", min_cells = 0
    ),
    "positive integer"
  )
})

test_that("typical 67k raw barcode matrices are detected", {
  features <- c(rep(5, 60000), rep(500, 7000))
  expect_true(scRDSreport:::.looks_like_raw_droplets(features, min_features = 200, n_cells = 67000))
  expect_false(scRDSreport:::.looks_like_raw_droplets(features[1:5000], min_features = 200, n_cells = 5000))
})
