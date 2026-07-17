test_that("analysis feature limiting preserves Seurat feature metadata identity", {
  skip_if_not_installed("SeuratObject")
  features <- sprintf("feature%02d", seq_len(12L))
  cells <- sprintf("cell%02d", seq_len(20L))
  counts <- matrix(0, nrow = length(features), ncol = length(cells),
                   dimnames = list(features, cells))
  for (index in seq_along(features)) counts[index, seq_len(index)] <- index
  object <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE)
  )
  feature_metadata <- data.frame(
    gene_symbol = paste0("symbol_", features),
    row.names = features
  )
  object[["RNA"]] <- SeuratObject::AddMetaData(object[["RNA"]], feature_metadata)

  limited <- scRDSreport:::.subset_top_detected_features(object, 10L)
  # Character subsetting preserves the assay's original row order; the ten
  # most-detected identities are feature03 through feature12.
  expected <- features[3:12]
  metadata <- as.data.frame(limited[["RNA"]][[]])

  expect_equal(rownames(limited), expected)
  expect_identical(rownames(metadata), expected)
  expect_equal(metadata$gene_symbol, paste0("symbol_", expected))
  expect_equal(
    rownames(scRDSreport:::.layer_data(limited, "RNA", "counts")),
    expected
  )
})
