.read_single_cell_rds <- function(input) {
  if (length(input) != 1L || !file.exists(input)) .sc_stop("Input RDS does not exist: %s", input)
  if (tolower(tools::file_ext(input)) != "rds") .sc_stop("Input must be an .rds file: %s", input)

  object <- tryCatch(readRDS(input), error = function(e) {
    .sc_stop("Could not read RDS '%s': %s", input, conditionMessage(e))
  })

  if (.is_seurat(object) || .is_sce(object) || is.matrix(object) || methods::is(object, "Matrix")) {
    return(object)
  }
  if (is.list(object)) {
    seurat_items <- vapply(object, .is_seurat, logical(1))
    if (sum(seurat_items) == 1L) return(object[[which(seurat_items)]])
    matrix_names <- intersect(names(object), c("counts", "count", "raw_counts", "matrix", "expr", "expression"))
    if (length(matrix_names)) return(object)
  }
  .sc_stop(
    "Unsupported RDS object class: %s. Supported inputs are Seurat, SingleCellExperiment, a matrix, or a list containing counts.",
    paste(class(object), collapse = "/")
  )
}

.as_seurat <- function(object, project = "scRDSreport") {
  if (.is_seurat(object)) {
    .require_optional("SeuratObject", "inspect and export a Seurat object")
    return(object)
  }
  if (.is_sce(object)) {
    .require_optional("Seurat", "convert a SingleCellExperiment object")
    .require_optional("SummarizedExperiment", "read assays from a SingleCellExperiment object")
    assays <- SummarizedExperiment::assayNames(object)
    counts_name <- if ("counts" %in% assays) "counts" else assays[[1L]]
    data_name <- if ("logcounts" %in% assays) "logcounts" else NULL
    return(Seurat::as.Seurat(object, counts = counts_name, data = data_name))
  }

  counts <- object
  metadata <- NULL
  if (is.list(object)) {
    key <- intersect(names(object), c("counts", "count", "raw_counts", "matrix", "expr", "expression"))[[1L]]
    counts <- object[[key]]
    meta_key <- intersect(names(object), c("metadata", "meta.data", "colData", "cell_metadata"))
    if (length(meta_key)) metadata <- as.data.frame(object[[meta_key[[1L]]]])
  }
  if (!(is.matrix(counts) || methods::is(counts, "Matrix"))) {
    .sc_stop("The counts element is not a matrix.")
  }
  .require_optional("SeuratObject", "create a Seurat object from a count matrix")
  SeuratObject::CreateSeuratObject(counts = counts, meta.data = metadata, project = project)
}

.join_split_layers <- function(object, verbose = TRUE) {
  if (!.is_seurat(object) ||
      !exists("JoinLayers", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    return(object)
  }
  for (assay in .seurat_assays(object)) {
    layers <- .assay_layers(object, assay)
    bases <- sub("[.].*$", "", layers)
    if (anyDuplicated(bases)) {
      .sc_message(verbose, "Joining split Seurat v5 layers in assay %s...", assay)
      object <- tryCatch(
        SeuratObject::JoinLayers(object, assay = assay),
        error = function(e) .sc_stop("Could not join split layers in assay '%s': %s", assay, conditionMessage(e))
      )
    }
  }
  object
}
