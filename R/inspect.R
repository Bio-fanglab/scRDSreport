.seurat_metadata <- function(object) {
  if (!.is_seurat(object)) return(data.frame(row.names = colnames(object)))
  meta <- tryCatch(as.data.frame(object[[]]), error = function(e) {
    meta <- .slot_or_null(object, "meta.data")
    if (is.null(meta)) data.frame(row.names = colnames(object)) else as.data.frame(meta)
  })
  active <- .slot_or_null(object, "active.ident")
  if (!is.null(active) && length(active) == nrow(meta) && !"active_ident" %in% names(meta)) {
    meta$active_ident <- as.character(active)
  }
  meta
}

.seurat_assays <- function(object) {
  if (!.is_seurat(object)) return("counts")
  out <- tryCatch(names(object@assays), error = function(e) character())
  if (!length(out) && requireNamespace("SeuratObject", quietly = TRUE)) {
    out <- tryCatch(SeuratObject::Assays(object), error = function(e) character())
  }
  out
}

.seurat_reductions <- function(object) {
  if (!.is_seurat(object)) return(character())
  tryCatch(names(object@reductions), error = function(e) character())
}

.seurat_graphs <- function(object) {
  if (!.is_seurat(object)) return(character())
  tryCatch(names(object@graphs), error = function(e) character())
}

.cluster_columns <- function(meta) {
  if (!ncol(meta)) return(character())
  candidates <- names(meta)[grepl(
    "(^|[._])(cluster|clusters|seurat_clusters)([._]|$)|clusters?$|_snn_res\\.",
    names(meta), ignore.case = TRUE
  )]
  if ("active_ident" %in% names(meta)) {
    is_sample_ident <- "orig.ident" %in% names(meta) && identical(as.character(meta$active_ident), as.character(meta$orig.ident))
    if (!is_sample_ident && length(unique(meta$active_ident)) > 1L) candidates <- c(candidates, "active_ident")
  }
  candidates <- unique(candidates)
  candidates[vapply(candidates, function(name) {
    values <- meta[[name]]
    length(unique(values[!is.na(values)])) >= 2L
  }, logical(1))]
}

.annotation_columns <- function(meta) {
  if (!ncol(meta)) return(character())
  candidates <- names(meta)[grepl(
    "cell.?type|annotation|annotated|predicted[._]?(id|label|celltype)|cell[._]?ontology|^labels?$",
    names(meta), ignore.case = TRUE
  )]
  candidates <- candidates[!grepl("^\\.scRDSreport_", candidates)]
  candidates[vapply(candidates, function(name) {
    values <- meta[[name]]
    n <- length(unique(as.character(values[!is.na(values) & nzchar(as.character(values))])))
    n >= 1L && n <= max(200L, ceiling(nrow(meta) * 0.5))
  }, logical(1))]
}

.primary_annotation_column <- function(meta, annotation_col = NULL) {
  if (!is.null(annotation_col)) {
    if (length(annotation_col) != 1L || !annotation_col %in% names(meta)) {
      .sc_stop("annotation_col '%s' is not present in cell metadata.", annotation_col)
    }
    return(annotation_col)
  }
  candidates <- .annotation_columns(meta)
  if (!length(candidates)) return(NULL)
  score <- vapply(candidates, function(name) {
    lname <- tolower(name)
    value <- 0
    if (grepl("manual|curated", lname)) value <- value + 50
    if (lname %in% c("celltype", "cell_type", "cell.type")) value <- value + 40
    if (grepl("cell.?type", lname)) value <- value + 30
    if (grepl("annotation", lname)) value <- value + 20
    if (grepl("fine|subtype", lname)) value <- value + 5
    value
  }, numeric(1))
  candidates[[which.max(score)]]
}

.has_normalized_data <- function(object) {
  if (!.is_seurat(object)) return(FALSE)
  assays <- .seurat_assays(object)
  any(vapply(assays, function(assay) {
    layers <- .assay_layers(object, assay)
    normalized <- layers[grepl("^(data|logcounts|normalized)([.]|$)", layers, ignore.case = TRUE)]
    any(vapply(normalized, function(layer) {
      value <- .layer_data(object, assay, layer)
      !is.null(value) && length(value) > 0L && nrow(value) > 0L && ncol(value) > 0L
    }, logical(1)))
  }, logical(1)))
}

.object_inventory <- function(object) {
  meta <- .seurat_metadata(object)
  assays <- .seurat_assays(object)
  reductions <- .seurat_reductions(object)
  graphs <- .seurat_graphs(object)
  cluster_cols <- .cluster_columns(meta)
  annotation_cols <- .annotation_columns(meta)
  commands <- if (.is_seurat(object)) names(.slot_or_null(object, "commands") %||% list()) else character()
  tools <- if (.is_seurat(object)) names(.slot_or_null(object, "tools") %||% list()) else character()
  misc <- if (.is_seurat(object)) names(.slot_or_null(object, "misc") %||% list()) else character()
  normalized <- .has_normalized_data(object)
  has_reduction <- length(reductions) > 0L
  has_cluster <- length(cluster_cols) > 0L || length(annotation_cols) > 0L

  status <- if (has_reduction && has_cluster) {
    "analyzed"
  } else if (normalized || has_reduction || has_cluster || length(graphs)) {
    "partial"
  } else {
    "raw"
  }

  data.frame(
    item = c("class", "cells", "features", "assays", "reductions", "graphs", "cluster_columns", "annotation_columns", "commands", "tools", "misc", "normalized_data", "analysis_status"),
    value = c(
      paste(class(object), collapse = "/"),
      if (!is.null(ncol(object))) ncol(object) else NA,
      if (!is.null(nrow(object))) nrow(object) else NA,
      paste(assays, collapse = ", "),
      paste(reductions, collapse = ", "),
      paste(graphs, collapse = ", "),
      paste(cluster_cols, collapse = ", "),
      paste(annotation_cols, collapse = ", "),
      paste(commands, collapse = ", "),
      paste(tools, collapse = ", "),
      paste(misc, collapse = ", "),
      normalized,
      status
    ),
    stringsAsFactors = FALSE
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Inspect a single-cell RDS file without changing it
#'
#' @param input Path to an RDS file.
#' @return A two-column data frame describing the object and its analysis state.
#' @export
inspect_rds <- function(input) {
  object <- .read_single_cell_rds(input)
  object <- .as_seurat(object)
  object <- .join_split_layers(object, verbose = FALSE)
  .object_inventory(object)
}
