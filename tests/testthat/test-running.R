test_that("running preserves input annotations and exports original and analysis matrices", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(sample.int(5, 240, replace = TRUE) - 1L, nrow = 20,
           dimnames = list(paste0("gene", seq_len(20)), paste0("cell", seq_len(12)))),
    sparse = TRUE
  )
  meta <- data.frame(
    sample = rep(c("Ctrl_rep1", "Ctrl_rep2"), each = 6),
    custom_labels = rep(c("TypeA", "TypeB"), 6),
    row.names = colnames(counts)
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)
  input <- tempfile(fileext = ".rds")
  output <- tempfile("scrdsreport-output-")
  saveRDS(object, input)

  result <- running(
    input = input, output = output, analyze = "never", render = FALSE,
    annotation_col = "custom_labels", matrix_layers = "counts"
  )
  manifest <- readRDS(result$manifest)
  expect_equal(manifest$annotation$source, "rds_metadata")
  expect_equal(manifest$annotation$primary_column, "custom_labels")
  expect_equal(manifest$download_embedding$mode, "auto")
  expect_true(any(grepl("original_RNA_counts[.]mtx[.]gz$", manifest$files$path)))
  expect_true(any(grepl("analysis_RNA_counts[.]mtx[.]gz$", manifest$files$path)))
  original_annotation_file <- manifest$files$path[
    manifest$files$section == "annotation_original"
  ]
  analysis_annotation_file <- manifest$files$path[
    manifest$files$section == "annotation_analysis"
  ]
  expect_length(original_annotation_file, 1L)
  expect_length(analysis_annotation_file, 1L)
  original_annotations <- utils::read.csv(
    file.path(output, original_annotation_file), check.names = FALSE
  )
  analysis_annotations <- utils::read.csv(
    file.path(output, analysis_annotation_file), check.names = FALSE
  )
  expect_equal(names(original_annotations), c("cell", "custom_labels"))
  expect_equal(names(analysis_annotations), c("cell", "custom_labels"))
  expect_equal(nrow(original_annotations), ncol(object))
  expect_equal(nrow(analysis_annotations), ncol(object))
  expect_false(any(grepl("scRDSreport_celltype", names(original_annotations))))
})

test_that("original metadata remains complete when the analysis object is a subset", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(sample.int(4, 300, replace = TRUE) - 1L, nrow = 20,
           dimnames = list(paste0("gene", seq_len(20)), paste0("cell", seq_len(15)))),
    sparse = TRUE
  )
  original <- SeuratObject::CreateSeuratObject(counts = counts)
  original$celltype_manual <- rep(c("A", "B", "C"), length.out = ncol(original))
  analysis <- original[, seq_len(10)]
  input <- tempfile(fileext = ".rds")
  output <- scRDSreport:::.ensure_output(tempfile("scrdsreport-export-"), FALSE)
  saveRDS(original, input)
  files <- scRDSreport:::.export_object(
    object = analysis, input = input, output = output,
    sample_design = data.frame(sample_id = "sample", group = "sample"),
    raw_matrix_object = original, matrix_object = analysis,
    original_annotation_columns = "celltype_manual",
    analysis_annotation_columns = "celltype_manual",
    matrix_layers = "counts", verbose = FALSE
  )
  original_meta <- utils::read.csv(
    file.path(output, files$path[files$section == "metadata_original"]),
    check.names = FALSE
  )
  analysis_meta <- utils::read.csv(
    file.path(output, files$path[files$section == "metadata_analysis"]),
    check.names = FALSE
  )
  expect_equal(nrow(original_meta), 15L)
  expect_equal(nrow(analysis_meta), 10L)
  expect_equal(original_meta$cell, colnames(original))
  expect_equal(analysis_meta$cell, colnames(analysis))
})

test_that("running validates named argument lists before reading input", {
  expect_error(
    running("missing.rds", tempfile(), scp_args = list(1)),
    "argument names"
  )
  expect_error(
    running("missing.rds", tempfile(), embed_max_mb = 0),
    "embed_max_mb"
  )
})
