test_that("numeric cluster identifiers use a discrete reduction-plot scale", {
  template <- system.file("quarto", "report.qmd", package = "scRDSreport")
  expect_true(nzchar(template))

  lines <- readLines(template, warn = FALSE, encoding = "UTF-8")
  expect_true(any(grepl(
    "embedding[[color_col]] <- factor(embedding[[color_col]])",
    lines,
    fixed = TRUE
  )))
  expect_true(any(grepl("scale_color_manual", lines, fixed = TRUE)))
})

test_that("a report renders with numeric cluster identifiers", {
  skip_on_cran()
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("DT")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggsci")
  skip_if_not_installed("htmltools")
  skip_if_not_installed("knitr")
  skip_if(!nzchar(Sys.which("quarto")), "Quarto CLI is not available")

  counts <- Matrix::Matrix(
    matrix(
      sample.int(5L, 240L, replace = TRUE) - 1L,
      nrow = 20L,
      dimnames = list(paste0("gene", seq_len(20L)), paste0("cell", seq_len(12L)))
    ),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$seurat_clusters <- rep(0:1, each = 6L)
  embeddings <- matrix(
    stats::rnorm(24L),
    ncol = 2L,
    dimnames = list(colnames(object), c("UMAP_1", "UMAP_2"))
  )
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    key = "UMAP_",
    assay = "RNA"
  )

  input <- tempfile(fileext = ".rds")
  output <- tempfile("scrdsreport-render-")
  saveRDS(object, input)

  result <- suppressWarnings(running(
    input = input,
    output = output,
    analyze = "never",
    run_markers = FALSE,
    matrix_layers = "counts",
    embed_downloads = "never",
    verbose = FALSE
  ))

  expect_true(file.exists(result$report))
  expect_gt(file.info(result$report)$size, 10000)
  html <- paste(readLines(result$report, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(html, "seurat_clusters", fixed = TRUE)
})
