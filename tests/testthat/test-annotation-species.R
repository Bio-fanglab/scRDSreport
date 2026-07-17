test_that("existing annotation columns are detected and custom columns are accepted", {
  meta <- data.frame(
    celltype_manual = c("T", "B", "T"),
    custom_labels = c("x", "y", "x"),
    nCount_RNA = c(10, 20, 30)
  )
  expect_equal(scRDSreport:::.annotation_columns(meta), "celltype_manual")
  expect_equal(scRDSreport:::.primary_annotation_column(meta), "celltype_manual")
  expect_equal(scRDSreport:::.primary_annotation_column(meta, "custom_labels"), "custom_labels")
  expect_error(scRDSreport:::.primary_annotation_column(meta, "missing"), "not present")
})

test_that("stable feature prefixes identify multiple species", {
  make_ids <- function(prefix) paste0(prefix, sprintf("%011d", 1:20))
  expect_equal(scRDSreport:::.detect_species(make_ids("ENSG"))$species, "human")
  expect_equal(scRDSreport:::.detect_species(make_ids("ENSMUSG"))$species, "mouse")
  expect_equal(scRDSreport:::.detect_species(make_ids("ENSRNOG"))$species, "rat")
  expect_equal(scRDSreport:::.detect_species(make_ids("ENSDARG"))$species, "zebrafish")
  expect_equal(scRDSreport:::.detect_species(c("Actb", "Gapdh", "Rpl3"))$species, "unknown")
})

test_that("Seurat v5 split layers are joined before analysis", {
  skip_if_not_installed("SeuratObject", minimum_version = "5.0.0")
  counts <- Matrix::Matrix(matrix(seq_len(24), nrow = 6,
                                  dimnames = list(paste0("g", 1:6), paste0("c", 1:4))), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(counts)
  object[["RNA"]] <- split(object[["RNA"]], f = c("a", "a", "b", "b"))
  expect_true(anyDuplicated(sub("[.].*$", "", SeuratObject::Layers(object[["RNA"]]))) > 0L)
  joined <- scRDSreport:::.join_split_layers(object, verbose = FALSE)
  expect_equal(SeuratObject::Layers(joined[["RNA"]]), "counts")
  expect_equal(ncol(SeuratObject::LayerData(joined, assay = "RNA", layer = "counts")), 4L)
})
