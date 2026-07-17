test_that("matrix bundle contains matrix, features, and barcodes", {
  matrix <- Matrix::Matrix(matrix(c(1, 0, 2, 3), nrow = 2,
                                  dimnames = list(c("g1", "g2"), c("c1", "c2"))), sparse = TRUE)
  root <- tempfile("matrix-export-")
  dir.create(root)
  manifest <- scRDSreport:::.write_matrix_bundle(matrix, root, "RNA_counts", root)
  expect_equal(nrow(manifest), 3L)
  expect_true(all(file.exists(file.path(root, manifest$path))))

  con <- gzfile(file.path(root, "RNA_counts.mtx.gz"), open = "rt")
  on.exit(close(con), add = TRUE)
  restored <- Matrix::readMM(con)
  expect_equal(as.matrix(restored), as.matrix(matrix), ignore_attr = TRUE)
  expect_equal(readLines(gzfile(file.path(root, "RNA_counts_features.tsv.gz"))), rownames(matrix))
  expect_equal(readLines(gzfile(file.path(root, "RNA_counts_barcodes.tsv.gz"))), colnames(matrix))
})

test_that("matrix preview uses features as rows and cells as columns", {
  matrix <- Matrix::Matrix(matrix(seq_len(12), nrow = 3,
                                  dimnames = list(paste0("g", 1:3), paste0("cell", 1:4))), sparse = TRUE)
  preview <- scRDSreport:::.matrix_preview_table(matrix, max_features = 2, max_cells = 3)
  expect_equal(preview$feature, c("g1", "g2"))
  expect_equal(names(preview), c("feature", "cell1", "cell2", "cell3"))
})

test_that("relative manifest paths stay below the report root", {
  root <- tempfile("report-root-")
  dir.create(file.path(root, "tables"), recursive = TRUE)
  path <- file.path(root, "tables", "x.csv.gz")
  writeLines("x", path)
  row <- scRDSreport:::.manifest_row("test", "x", path, root)
  expect_equal(row$path, "tables/x.csv.gz")
})
