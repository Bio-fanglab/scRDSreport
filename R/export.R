.assay_layers <- function(object, assay) {
  if (!.is_seurat(object)) return("counts")
  assay_object <- tryCatch(object[[assay]], error = function(e) NULL)
  if (is.null(assay_object)) return(character())
  if (requireNamespace("SeuratObject", quietly = TRUE) && exists("Layers", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    layers <- tryCatch(SeuratObject::Layers(assay_object), error = function(e) character())
    if (length(layers)) return(layers)
  }
  intersect(c("counts", "data", "scale.data"), methods::slotNames(assay_object))
}

.layer_data <- function(object, assay, layer) {
  if (!.is_seurat(object)) return(object)
  if (requireNamespace("SeuratObject", quietly = TRUE) && exists("LayerData", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    value <- tryCatch(SeuratObject::LayerData(object, assay = assay, layer = layer), error = function(e) NULL)
    if (!is.null(value)) return(value)
  }
  tryCatch(SeuratObject::GetAssayData(object, assay = assay, slot = layer), error = function(e) NULL)
}

.write_lines_gz <- function(x, path) {
  con <- gzfile(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(as.character(x), con, useBytes = TRUE)
  invisible(path)
}

.write_matrix_bundle <- function(matrix, directory, prefix, root) {
  if (is.null(matrix) || !length(matrix)) return(data.frame())
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  sparse <- if (methods::is(matrix, "sparseMatrix")) matrix else Matrix::Matrix(matrix, sparse = TRUE)
  prefix <- .safe_name(prefix)
  temporary <- tempfile(fileext = ".mtx")
  on.exit(unlink(temporary), add = TRUE)
  Matrix::writeMM(sparse, temporary)
  matrix_path <- file.path(directory, paste0(prefix, ".mtx.gz"))
  feature_path <- file.path(directory, paste0(prefix, "_features.tsv.gz"))
  barcode_path <- file.path(directory, paste0(prefix, "_barcodes.tsv.gz"))
  .gzip_file(temporary, matrix_path)
  .write_lines_gz(rownames(sparse) %||% seq_len(nrow(sparse)), feature_path)
  .write_lines_gz(colnames(sparse) %||% seq_len(ncol(sparse)), barcode_path)
  rbind(
    .manifest_row("expression_matrix", paste0(prefix, " matrix (features x cells)"), matrix_path, root,
                  nrow(sparse), ncol(sparse), "Rows=features; columns=cells; Matrix Market coordinate format"),
    .manifest_row("expression_matrix", paste0(prefix, " row names (features)"), feature_path, root,
                  nrow(sparse), 1L, "One feature name per matrix row, in identical order"),
    .manifest_row("expression_matrix", paste0(prefix, " column names (cells)"), barcode_path, root,
                  ncol(sparse), 1L, "One cell name per matrix column, in identical order")
  )
}

.matrix_preview_table <- function(matrix, max_features = 12L, max_cells = 12L) {
  if (is.null(matrix) || !nrow(matrix) || !ncol(matrix)) return(NULL)
  nr <- min(as.integer(max_features), nrow(matrix))
  nc <- min(as.integer(max_cells), ncol(matrix))
  block <- as.matrix(matrix[seq_len(nr), seq_len(nc), drop = FALSE])
  data.frame(feature = rownames(block) %||% as.character(seq_len(nr)),
             block, row.names = NULL, check.names = FALSE)
}

.write_table <- function(x, path, section, label, root, note = "") {
  x <- as.data.frame(x)
  rn <- rownames(x)
  default_rn <- is.null(rn) || identical(rn, as.character(seq_len(nrow(x))))
  if (!default_rn && !"rowname" %in% names(x)) {
    x <- data.frame(rowname = rn, x, row.names = NULL, check.names = FALSE)
  }
  .write_csv_gz(x, path, row.names = FALSE)
  .manifest_row(section, label, path, root, nrow(x), ncol(x), note)
}

.export_metadata_tables <- function(object, output, prefix,
                                    annotation_columns = NULL) {
  meta <- .seurat_metadata(object)
  if (!nrow(meta)) return(data.frame())
  prefix <- .safe_name(prefix)
  cell_names <- rownames(meta)
  annotation_columns <- unique(annotation_columns %||% character())
  annotation_columns <- annotation_columns[annotation_columns %in% names(meta)]
  exported_meta <- data.frame(cell = cell_names, meta, check.names = FALSE, row.names = NULL)
  meta_path <- file.path(output, "tables", paste0(prefix, "_cell_metadata.csv.gz"))
  rows <- .write_table(
    exported_meta, meta_path, paste0("metadata_", prefix),
    paste(tools::toTitleCase(prefix), "cell metadata"), output,
    note = paste0("One row per ", prefix, " object cell; all metadata columns are retained.")
  )
  if (length(annotation_columns)) {
    annotations <- data.frame(
      cell = cell_names, meta[annotation_columns],
      check.names = FALSE, row.names = NULL
    )
    annotation_path <- file.path(
      output, "tables", paste0(prefix, "_cell_annotations.csv.gz")
    )
    rows <- rbind(rows, .write_table(
      annotations, annotation_path, paste0("annotation_", prefix),
      paste(tools::toTitleCase(prefix), "cell annotations"), output,
      note = if (identical(prefix, "original")) {
        "One row per original RDS cell; annotation values are preserved from the input RDS."
      } else {
        paste0(
          "One row per analysis object cell; columns may contain preserved input annotations ",
          "and explicitly generated reference/manual annotations. See the annotation source in manifest."
        )
      }
    ))
  }
  rows
}

.embedding_table <- function(object, reduction) {
  embedding <- tryCatch(SeuratObject::Embeddings(object, reduction = reduction), error = function(e) NULL)
  if (is.null(embedding)) return(NULL)
  data.frame(cell = rownames(embedding), embedding, check.names = FALSE, row.names = NULL)
}

.variable_feature_table <- function(object, assay) {
  features <- tryCatch(SeuratObject::VariableFeatures(object[[assay]]), error = function(e) character())
  if (!length(features)) return(NULL)
  data.frame(feature = features, rank = seq_along(features), stringsAsFactors = FALSE)
}

.collect_result_tables <- function(x, prefix = "result", depth = 0L, max_depth = 3L) {
  if (depth > max_depth || is.null(x)) return(list())
  if (is.data.frame(x)) return(stats::setNames(list(x), prefix))
  if (is.matrix(x) && length(x) <= 5e6) return(stats::setNames(list(as.data.frame(x)), prefix))
  if (isS4(x)) {
    output <- list()
    for (slot in methods::slotNames(x)) {
      child <- .collect_result_tables(
        methods::slot(x, slot), paste(prefix, .safe_name(slot), sep = "__"), depth + 1L, max_depth
      )
      output <- c(output, child)
    }
    return(output)
  }
  if (!is.list(x)) return(list())
  output <- list()
  nms <- names(x) %||% paste0("item", seq_along(x))
  for (i in seq_along(x)) {
    child <- .collect_result_tables(x[[i]], paste(prefix, .safe_name(nms[i]), sep = "__"), depth + 1L, max_depth)
    output <- c(output, child)
  }
  output
}

.export_object <- function(object, input, output, sample_design,
                           raw_matrix_object = object, matrix_object = object,
                           original_annotation_columns = NULL,
                           analysis_annotation_columns = NULL,
                           matrix_layers = c("counts", "data"), verbose = TRUE) {
  manifest <- data.frame()
  original_path <- file.path(output, "downloads", paste0("original_", basename(input)))
  if (!file.copy(input, original_path, overwrite = TRUE)) .sc_stop("Could not copy the original RDS to the report directory.")
  manifest <- rbind(manifest, .manifest_row("rds", "Original input RDS", original_path, output, note = "Unmodified input"))

  analyzed_path <- file.path(output, "downloads", "analysis_object.rds")
  saveRDS(object, analyzed_path, compress = "gzip")
  manifest <- rbind(manifest, .manifest_row("rds", "Report analysis RDS", analyzed_path, output, note = "Includes inferred design and any SCP results"))

  design_path <- file.path(output, "tables", "sample_design.csv.gz")
  manifest <- rbind(manifest, .write_table(
    sample_design, design_path, "design", "Sample design", output,
    note = paste0(
      "One row per sample candidate. n_cells_post_qc is measured after barcode QC and before ",
      "the optional analysis subset; n_cells_analysis is the number represented in analysis results."
    )
  ))

  manifest <- rbind(manifest, .export_metadata_tables(
    raw_matrix_object, output, "original",
    annotation_columns = original_annotation_columns
  ))
  manifest <- rbind(manifest, .export_metadata_tables(
    object, output, "analysis",
    annotation_columns = analysis_annotation_columns
  ))

  inventory <- .object_inventory(object)
  inventory_path <- file.path(output, "tables", "object_inventory.csv.gz")
  manifest <- rbind(manifest, .write_table(inventory, inventory_path, "inventory", "Object inventory", output))

  for (assay in .seurat_assays(raw_matrix_object)) {
    available <- .assay_layers(raw_matrix_object, assay)
    raw_layers <- available[sub("[.].*$", "", tolower(available)) == "counts"]
    for (layer in raw_layers) {
      .sc_message(verbose, "Exporting original assay %s layer %s...", assay, layer)
      matrix <- .layer_data(raw_matrix_object, assay, layer)
      prefix <- paste("original", assay, layer, sep = "_")
      manifest <- rbind(manifest, .write_matrix_bundle(
        matrix, file.path(output, "matrices"), prefix, output
      ))
      preview <- .matrix_preview_table(matrix)
      if (!is.null(preview)) {
        preview_path <- file.path(output, "tables", paste0("matrix_preview_", .safe_name(prefix), ".csv.gz"))
        manifest <- rbind(manifest, .write_table(
          preview, preview_path, "matrix_preview", paste("Original matrix preview:", assay, layer), output,
          note = "Preview only. Rows=features; columns after 'feature'=original cell names."
        ))
      }
    }
  }

  for (assay in .seurat_assays(matrix_object)) {
    available <- .assay_layers(matrix_object, assay)
    layer_base <- sub("[.].*$", "", tolower(available))
    wanted <- available[layer_base %in% tolower(matrix_layers)]
    for (layer in wanted) {
      .sc_message(verbose, "Exporting assay %s layer %s...", assay, layer)
      matrix <- .layer_data(matrix_object, assay, layer)
      prefix <- paste("analysis", assay, layer, sep = "_")
      manifest <- rbind(manifest, .write_matrix_bundle(
        matrix, file.path(output, "matrices"), prefix, output
      ))
      preview <- .matrix_preview_table(matrix)
      if (!is.null(preview)) {
        preview_path <- file.path(
          output, "tables", paste0("matrix_preview_", .safe_name(prefix), ".csv.gz")
        )
        manifest <- rbind(manifest, .write_table(
          preview, preview_path, "matrix_preview",
          paste("Matrix preview:", assay, layer), output,
          note = "Preview only. Rows=features; columns after 'feature'=cell names."
        ))
      }
    }
    vf <- .variable_feature_table(object, assay)
    if (!is.null(vf)) {
      path <- file.path(output, "tables", paste0("variable_features_", .safe_name(assay), ".csv.gz"))
      manifest <- rbind(manifest, .write_table(vf, path, "features", paste("Variable features:", assay), output))
    }
    feature_meta <- tryCatch(as.data.frame(object[[assay]][[]]), error = function(e) NULL)
    if (!is.null(feature_meta) && nrow(feature_meta) && ncol(feature_meta)) {
      feature_meta <- data.frame(feature = rownames(feature_meta), feature_meta, row.names = NULL, check.names = FALSE)
      path <- file.path(output, "tables", paste0("feature_metadata_", .safe_name(assay), ".csv.gz"))
      manifest <- rbind(manifest, .write_table(feature_meta, path, "features", paste("Feature metadata:", assay), output))
    }
  }

  for (reduction in .seurat_reductions(object)) {
    embedding <- .embedding_table(object, reduction)
    if (!is.null(embedding)) {
      path <- file.path(output, "tables", paste0("embedding_", .safe_name(reduction), ".csv.gz"))
      manifest <- rbind(manifest, .write_table(embedding, path, "embedding", paste("Embedding:", reduction), output))
    }
  }

  if (.is_seurat(object)) {
    containers <- list(misc = .slot_or_null(object, "misc"), tools = .slot_or_null(object, "tools"))
    result_tables <- .collect_result_tables(containers)
    for (name in names(result_tables)) {
      path <- file.path(output, "tables", paste0(.safe_name(name), ".csv.gz"))
      manifest <- rbind(manifest, .write_table(result_tables[[name]], path, "analysis_result", name, output))
    }
    commands <- .slot_or_null(object, "commands")
    if (length(commands)) {
      command_table <- data.frame(command = names(commands), stringsAsFactors = FALSE)
      path <- file.path(output, "tables", "seurat_commands.csv.gz")
      manifest <- rbind(manifest, .write_table(command_table, path, "provenance", "Seurat command history", output))
    }
  }
  manifest
}
