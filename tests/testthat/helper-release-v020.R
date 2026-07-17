.release_v020_object <- function(n_features = 48L, n_cells = 24L) {
  skip_if_not_installed("SeuratObject")
  set.seed(20260717)
  counts <- matrix(
    stats::rpois(n_features * n_cells, lambda = 3),
    nrow = n_features,
    dimnames = list(
      c("MT-ND1", "RPS3", "HBA1", paste0("gene", seq_len(n_features - 3L))),
      paste0("cell", seq_len(n_cells))
    )
  )
  counts[counts < 2L] <- 0L
  object <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE)
  )
  object$sample <- rep(c("Ctrl_rep1", "Treat_rep1"), length.out = n_cells)
  object$group <- rep(c("Ctrl", "Treat"), length.out = n_cells)
  object$celltype_manual <- rep(c("TypeA", "TypeB"), each = ceiling(n_cells / 2L))[seq_len(n_cells)]
  object$seurat_clusters <- rep(c("0", "1"), each = ceiling(n_cells / 2L))[seq_len(n_cells)]

  embedding <- matrix(
    stats::rnorm(n_cells * 2L),
    ncol = 2L,
    dimnames = list(colnames(object), c("UMAP_1", "UMAP_2"))
  )
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "UMAP_",
    assay = "RNA"
  )
  object
}
