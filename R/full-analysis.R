# Optional full-analysis layer -------------------------------------------------
#
# This file deliberately keeps every advanced dependency optional.  The main
# package can therefore inspect and export an RDS with its core dependencies,
# while this layer records a precise skipped/needs_input state when a requested
# scientific module cannot be run safely.

utils::globalVariables(c(
  "minus_log10", "Description", "plot_size", "dimension_1", "dimension_2",
  "Phase", "score_value", "fraction", "cluster", "n_cells",
  "mean_absolute_cnv_signal", "target", "interaction_weight", "value",
  "nCount", "nFeature", "pass", "color_value", "group", "symbol",
  "z_expression"
))

.fa_module_aliases <- function(id) {
  aliases <- list(
    qc = c("qc", "quality_control"),
    reduction = c("reduction", "dimensionality", "dimensional_reduction"),
    cluster = c("cluster", "clustering"),
    celltype = c("celltype", "cell_type", "annotation", "composition"),
    differential = c("differential", "de", "differential_expression"),
    enrichment = c("enrichment", "ora", "gsea", "gsva"),
    pseudotime = c("pseudotime", "trajectory", "monocle3"),
    communication = c("communication", "cellchat", "cell_communication"),
    cell_cycle = c("cell_cycle", "cellcycle"),
    tf = c("tf", "tf_expression", "transcription_factor"),
    cnv = c("cnv", "infercnv", "infer_cnv"),
    downloads = c("downloads", "artifacts")
  )
  aliases[[id]] %||% id
}

.fa_merge_list <- function(x, y) {
  if (!is.list(x)) x <- list()
  if (is.logical(y) && length(y) == 1L && !is.na(y)) y <- list(enabled = y)
  if (!is.list(y)) return(x)
  utils::modifyList(x, y, keep.null = TRUE)
}

.fa_module_config <- function(config, id) {
  if (!is.list(config)) config <- list()
  out <- list()
  containers <- list(config, config$modules, config$module_options, config$analysis)
  aliases <- .fa_module_aliases(id)
  for (container in containers) {
    if (!is.list(container)) next
    for (alias in aliases) {
      if (!is.null(container[[alias]])) out <- .fa_merge_list(out, container[[alias]])
    }
  }
  if (identical(id, "celltype") && is.null(out$mapping) && is.null(out$manual_markers) &&
      !is.null(config$resource_overrides$manual_markers)) {
    out$manual_markers <- config$resource_overrides$manual_markers
  }
  overrides <- config$resource_overrides %||% list()
  resource_fields <- switch(
    id,
    qc = c(
      "orgdb", "feature_keytype", "symbol_column", "mitochondrial_pattern",
      "ribosomal_pattern", "hemoglobin_pattern", "pattern_ignore_case"
    ),
    celltype = c(
      "manual_markers", "auto_annotation_reference", "orgdb",
      "feature_keytype", "symbol_column"
    ),
    enrichment = c(
      "orgdb", "feature_keytype", "symbol_column", "kegg_code",
      "gene_sets", "gmt", "gmt_files", "scientific_name",
      "msigdbr_species", "msigdbr_db_species", "msigdbr_default_collection",
      "msigdbr_ortholog_projection", "gene_sets_strategy"
    ),
    communication = c(
      "cellchat_db", "orgdb", "feature_keytype", "symbol_column"
    ),
    cell_cycle = c(
      "s_genes", "g2m_genes", "cell_cycle_genes", "cell_cycle_strategy",
      "scientific_name", "orgdb", "feature_keytype", "symbol_column"
    ),
    tf = c(
      "tf_genes", "tf_catalog", "tf_catalog_strategy", "scientific_name",
      "orgdb", "feature_keytype", "symbol_column"
    ),
    cnv = c(
      "gene_order", "gtf", "txdb", "orgdb", "feature_keytype",
      "symbol_column"
    ),
    character()
  )
  for (field in resource_fields) {
    if (is.null(out[[field]]) && !is.null(overrides[[field]])) out[[field]] <- overrides[[field]]
  }
  limits <- config$limits %||% list()
  if (is.null(out$max_plot_cells) && !is.null(limits$plot_max_cells)) out$max_plot_cells <- limits$plot_max_cells
  if (is.null(out$cores) && !is.null(limits$workers)) out$cores <- limits$workers
  if (is.null(out$num_threads) && !is.null(limits$workers)) out$num_threads <- limits$workers
  if (is.null(out$min_cells) && !is.null(limits$min_cells_per_group)) out$min_cells <- limits$min_cells_per_group
  if (identical(id, "differential") && is.null(out$marker_max_cells_per_ident) &&
      !is.null(limits$marker_max_cells_per_ident)) {
    out$marker_max_cells_per_ident <- limits$marker_max_cells_per_ident
  }
  out
}

.fa_config_value <- function(config, name, default = NULL) {
  value <- config[[name]]
  if (is.null(value)) default else value
}

.fa_requested <- function(config, id, default = TRUE) {
  cfg <- .fa_module_config(config, id)
  value <- NULL
  modules <- config$modules
  if (!is.null(modules) && !is.list(modules) && !is.null(names(modules))) {
    aliases <- .fa_module_aliases(id)
    hit <- aliases[aliases %in% names(modules)]
    if (length(hit)) value <- modules[[hit[[1L]]]]
  }
  value <- value %||% cfg$enabled %||% cfg$requested %||% default
  isTRUE(value)
}

.fa_pkg_available <- function(package) {
  requireNamespace(package, quietly = TRUE)
}

.fa_pkg_fun <- function(package, name, exported = TRUE) {
  if (!.fa_pkg_available(package)) return(NULL)
  if (isTRUE(exported)) {
    tryCatch(getExportedValue(package, name), error = function(e) NULL)
  } else {
    tryCatch(get(name, envir = asNamespace(package), inherits = FALSE), error = function(e) NULL)
  }
}

.fa_pkg_object <- function(package, name) {
  if (!.fa_pkg_available(package)) return(NULL)
  object <- tryCatch(get(name, envir = asNamespace(package), inherits = FALSE), error = function(e) NULL)
  if (!is.null(object)) return(object)
  environment <- new.env(parent = emptyenv())
  suppressWarnings(try(utils::data(list = name, package = package, envir = environment), silent = TRUE))
  if (exists(name, envir = environment, inherits = FALSE)) get(name, envir = environment) else NULL
}

.fa_package_version <- function(package) {
  if (!.fa_pkg_available(package)) return(NA_character_)
  as.character(utils::packageVersion(package))
}

.fa_default_assay <- function(object) {
  value <- tryCatch(SeuratObject::DefaultAssay(object), error = function(e) NULL)
  if (is.null(value) || !length(value) || !nzchar(value[[1L]])) {
    assays <- .seurat_assays(object)
    if (length(assays)) assays[[1L]] else NULL
  } else {
    value[[1L]]
  }
}

.fa_matrix <- function(object, layer = "counts", assay = NULL) {
  assay <- assay %||% .fa_default_assay(object)
  if (is.null(assay)) return(NULL)
  value <- .layer_data(object, assay = assay, layer = layer)
  if (is.null(value) || !nrow(value) || !ncol(value)) NULL else value
}

.fa_species <- function(species_info, config = list()) {
  override <- config$resource_overrides$species %||% config$species %||% NULL
  if (is.list(species_info)) {
    selected <- species_info$selected %||% species_info$species %||% species_info$detected
    selected_value <- tolower(as.character(selected %||% "unknown")[[1L]])
    argument_was_auto <- identical(tolower(as.character(species_info$requested %||% "auto")[[1L]]), "auto")
    if (!is.null(override) && (argument_was_auto || selected_value %in% c("auto", "unknown"))) selected <- override
  } else {
    selected <- species_info
    selected_value <- tolower(as.character(selected %||% "unknown")[[1L]])
    if (!is.null(override) && selected_value %in% c("auto", "unknown")) selected <- override
  }
  selected <- selected %||% "unknown"
  .normalize_species_name(as.character(selected[[1L]]))
}

.fa_species_resources <- function(species, config = list()) {
  public <- species_resources(species)
  base <- list(
    orgdb = public$orgdb,
    feature_keytype = public$feature_keytype,
    symbol_column = public$symbol_column,
    kegg = public$kegg_code,
    cellchat = public$cellchat_db,
    txdb = public$txdb,
    gtf = public$gtf,
    gene_sets = public$gene_sets,
    gene_sets_strategy = public$gene_sets_strategy,
    msigdbr_species = public$msigdbr_species,
    msigdbr_db_species = public$msigdbr_db_species,
    msigdbr_default_collection = public$msigdbr_default_collection,
    scientific_name = public$scientific_name,
    taxonomy_id = public$taxonomy_id,
    auto_annotation_reference = public$auto_annotation_reference,
    cell_cycle_strategy = public$cell_cycle_strategy,
    tf_catalog_strategy = public$tf_catalog_strategy,
    mitochondrial_pattern = public$mitochondrial_pattern,
    ribosomal_pattern = public$ribosomal_pattern,
    hemoglobin_pattern = public$hemoglobin_pattern,
    pattern_ignore_case = public$pattern_ignore_case
  )
  supplied <- config$species_resources %||% list()
  if (is.list(supplied[[species]])) supplied <- supplied[[species]]
  if (is.list(supplied)) {
    if (is.null(supplied$kegg) && !is.null(supplied$kegg_code)) supplied$kegg <- supplied$kegg_code
    if (is.null(supplied$cellchat) && !is.null(supplied$cellchat_db)) supplied$cellchat <- supplied$cellchat_db
  }
  base <- .fa_merge_list(base, supplied)
  direct <- list(
    orgdb = config$orgdb,
    feature_keytype = config$feature_keytype,
    symbol_column = config$symbol_column,
    kegg = config$kegg %||% config$kegg_code,
    cellchat = config$cellchat %||% config$cellchat_db,
    txdb = config$txdb,
    gtf = config$gtf,
    gene_sets = config$gene_sets,
    gene_sets_strategy = config$gene_sets_strategy,
    scientific_name = config$scientific_name,
    taxonomy_id = config$taxonomy_id,
    msigdbr_species = config$msigdbr_species,
    msigdbr_db_species = config$msigdbr_db_species,
    msigdbr_default_collection = config$msigdbr_default_collection,
    msigdbr_ortholog_projection = config$msigdbr_ortholog_projection,
    auto_annotation_reference = config$auto_annotation_reference,
    cell_cycle_strategy = config$cell_cycle_strategy,
    tf_catalog_strategy = config$tf_catalog_strategy,
    mitochondrial_pattern = config$mitochondrial_pattern,
    ribosomal_pattern = config$ribosomal_pattern,
    hemoglobin_pattern = config$hemoglobin_pattern,
    pattern_ignore_case = config$pattern_ignore_case
  )
  output <- .fa_merge_list(base, direct[!vapply(direct, is.null, logical(1))])
  if (!is.null(output$msigdbr_db_species) &&
      is.null(config$msigdbr_ortholog_projection)) {
    database <- toupper(as.character(output$msigdbr_db_species[[1L]]))
    if (is.null(config$msigdbr_default_collection)) {
      output$msigdbr_default_collection <- if (identical(database, "MM")) "MH" else "H"
    }
    output$msigdbr_ortholog_projection <- !(
      (identical(species, "human") && identical(database, "HS")) ||
        (identical(species, "mouse") && identical(database, "MM"))
    )
    if (is.null(config$gene_sets_strategy)) {
      target <- output$scientific_name %||% species
      output$gene_sets_strategy <- paste0(
        "msigdbr_", tolower(database), "_to_",
        gsub("[^a-z0-9]+", "_", tolower(target))
      )
    }
  }
  output
}

.fa_feature_symbols <- function(object, assay = NULL) {
  assay <- assay %||% .fa_default_assay(object)
  features <- rownames(object)
  if (is.null(features)) return(character())
  feature_meta <- tryCatch(as.data.frame(object[[assay]][[]]), error = function(e) NULL)
  symbol_columns <- c("gene_symbols", "gene_symbol", "symbol", "gene_name", "gene_names")
  column <- if (!is.null(feature_meta)) intersect(symbol_columns, names(feature_meta)) else character()
  symbols <- if (length(column)) as.character(feature_meta[[column[[1L]]]]) else features
  missing <- is.na(symbols) | !nzchar(symbols)
  symbols[missing] <- features[missing]
  stats::setNames(symbols, features)
}

.fa_match_features <- function(requested, feature_symbols) {
  if (!length(requested) || !length(feature_symbols)) return(character())
  map <- match(toupper(as.character(requested)), toupper(as.character(feature_symbols)))
  unique(names(feature_symbols)[stats::na.omit(map)])
}

.fa_cluster_column <- function(object, config = list()) {
  meta <- .seurat_metadata(object)
  explicit <- config$cluster_column %||% config$column
  if (!is.null(explicit) && explicit %in% names(meta)) return(explicit)
  candidates <- .cluster_columns(meta)
  preferred <- c(".scRDSreport_cluster", "seurat_clusters")
  hit <- intersect(preferred, candidates)
  if (length(hit)) hit[[1L]] else if (length(candidates)) candidates[[1L]] else NULL
}

.fa_annotation_column <- function(object, config = list()) {
  meta <- .seurat_metadata(object)
  explicit <- config$annotation_column %||% config$column
  if (!is.null(explicit) && explicit %in% names(meta)) return(explicit)
  misc <- .slot_or_null(object, "misc") %||% list()
  active <- misc$scRDSreport_full_analysis$active_annotation_column %||% NULL
  if (!is.null(active) && length(active) == 1L && active %in% names(meta)) {
    active_values <- as.character(meta[[active]])
    active_levels <- unique(active_values[!is.na(active_values) & nzchar(active_values)])
    if (length(active_levels) >= 1L && length(active_levels) <= max(200L, ceiling(nrow(meta) * 0.5))) {
      return(active)
    }
  }
  generated <- names(meta)[grepl("^\\.scRDSreport_celltype_", names(meta), ignore.case = TRUE)]
  generated <- generated[vapply(generated, function(name) {
    values <- as.character(meta[[name]])
    n <- length(unique(values[!is.na(values) & nzchar(values)]))
    n >= 1L && n <= max(200L, ceiling(nrow(meta) * 0.5))
  }, logical(1))]
  if (length(generated)) return(generated[[length(generated)]])
  .primary_annotation_column(meta)
}

.fa_sample_column <- function(object, config = list()) {
  meta <- .seurat_metadata(object)
  explicit <- config$sample_column
  candidates <- unique(c(explicit, ".scRDSreport_sample", "sample", "sample_id", "orig.ident"))
  hit <- candidates[!is.na(candidates) & candidates %in% names(meta)]
  if (length(hit)) hit[[1L]] else NULL
}

.fa_group_column <- function(object, config = list()) {
  meta <- .seurat_metadata(object)
  explicit <- config$group_column
  explicit_hit <- unique(explicit[!is.na(explicit) & explicit %in% names(meta)])
  if (length(explicit_hit)) return(explicit_hit[[1L]])

  inferred_column <- ".scRDSreport_group"
  has_inferred <- inferred_column %in% names(meta)
  user_mapped <- FALSE
  if (has_inferred) {
    provenance <- list(
      c(".scRDSreport_grouping_rule", "user_map"),
      c(".scRDSreport_design_confidence", "user")
    )
    user_mapped <- any(vapply(provenance, function(spec) {
      if (!spec[[1L]] %in% names(meta)) return(FALSE)
      values <- tolower(trimws(as.character(meta[[spec[[1L]]]])))
      values <- values[!is.na(values) & nzchar(values)]
      length(values) > 0L && all(values == spec[[2L]])
    }, logical(1)))
  }
  if (has_inferred && user_mapped) return(inferred_column)

  original_hit <- intersect(c("group", "condition", "treatment"), names(meta))
  if (length(original_hit)) return(original_hit[[1L]])
  if (has_inferred) inferred_column else NULL
}

.fa_palette <- function(n) {
  colors <- c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
    "#F0E442", "#000000", "#332288", "#88CCEE", "#44AA99", "#117733",
    "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#661100",
    "#6699CC", "#AA4466", "#4477AA", "#228833", "#EE6677", "#BBBBBB"
  )
  if (n <= length(colors)) return(colors[seq_len(max(0L, n))])
  grDevices::hcl.colors(n, palette = "Dark 3")
}

.fa_module_directory <- function(output, module) {
  path <- file.path(output, "analysis", .safe_name(module))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

.fa_sha256 <- function(path) {
  digest_fun <- .fa_pkg_fun("digest", "digest")
  if (is.null(digest_fun) || !file.exists(path)) return(NA_character_)
  tryCatch(
    digest_fun(object = path, algo = "sha256", file = TRUE, serialize = FALSE),
    error = function(e) NA_character_
  )
}

.fa_column_dictionary <- function(x, supplied = NULL) {
  defaults <- vapply(names(x), function(name) {
    lname <- tolower(name)
    if (lname %in% c("cell", "cell_id")) return("Unique cell barcode or cell name.")
    if (lname %in% c("feature", "gene", "gene_id", "symbol")) return("Feature or gene identifier.")
    if (grepl("sample", lname)) return("Biological sample identifier.")
    if (grepl("group|condition", lname)) return("Experimental group or condition.")
    if (grepl("cluster", lname)) return("Cluster identity from the analysis object.")
    if (grepl("cell.?type|annotation|label", lname)) return("Cell annotation retained from the RDS or explicitly requested annotation step.")
    if (grepl("p_val_adj|padj|fdr", lname)) return("Multiple-testing adjusted P value.")
    if (grepl("p_val|pvalue|p.value", lname)) return("Unadjusted P value.")
    if (grepl("logfc|log2fc|fold", lname)) return("Log2 fold change for the stated contrast.")
    if (grepl("fraction|proportion|percent|pct", lname)) return("Fraction, proportion, or percentage defined by the column name.")
    if (grepl("count|cells|n_", lname)) return("Count defined by the column name.")
    if (grepl("reason", lname)) return("Machine-readable explanation for this row or decision.")
    "Value defined by the source object or analysis module; see the artifact description."
  }, character(1))
  if (is.null(supplied)) return(as.list(defaults))
  supplied <- unlist(supplied, use.names = TRUE)
  defaults[names(supplied)] <- as.character(supplied)
  as.list(defaults)
}

.fa_artifact_record <- function(module, type, path, output, label, description,
                                rows = NA_integer_, columns = NA_integer_,
                                row_unit = NA_character_, column_unit = NA_character_,
                                column_dictionary = list(), units = list(), complete = TRUE) {
  relative <- .relative_path(path, output)
  extension <- if (grepl("[.]csv[.]gz$", path, ignore.case = TRUE)) {
    "csv.gz"
  } else if (grepl("[.]mtx[.]gz$", path, ignore.case = TRUE)) {
    "mtx.gz"
  } else {
    tolower(tools::file_ext(path))
  }
  list(
    artifact_id = paste(.safe_name(module), .safe_name(basename(path)), sep = "__"),
    module = module,
    type = type,
    format = extension,
    path = relative,
    path_is_relative = TRUE,
    label = label,
    description = description,
    rows = as.integer(rows),
    columns = as.integer(columns),
    row_unit = row_unit,
    column_unit = column_unit,
    column_dictionary = column_dictionary,
    units = units,
    bytes = if (file.exists(path)) as.numeric(file.info(path)$size) else NA_real_,
    sha256 = .fa_sha256(path),
    complete = isTRUE(complete)
  )
}

.fa_write_table_artifact <- function(x, output, module, name, label, description,
                                     row_unit, column_dictionary = NULL,
                                     units = list()) {
  directory <- .fa_module_directory(output, module)
  path <- file.path(directory, paste0(.safe_name(name), ".csv.gz"))
  x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  .write_csv_gz(x, path, row.names = FALSE)
  .fa_artifact_record(
    module, "table", path, output, label, description,
    rows = nrow(x), columns = ncol(x), row_unit = row_unit,
    column_unit = "field", column_dictionary = .fa_column_dictionary(x, column_dictionary),
    units = units
  )
}

.fa_write_rds_artifact <- function(x, output, module, name, label, description) {
  directory <- .fa_module_directory(output, module)
  path <- file.path(directory, paste0(.safe_name(name), ".rds"))
  saveRDS(x, path, compress = "gzip")
  .fa_artifact_record(
    module, "rds", path, output, label, description,
    row_unit = "R object", column_unit = NA_character_, column_dictionary = list()
  )
}

.fa_write_plot_artifacts <- function(plot, output, module, name, label, description,
                                     width = 8, height = 6) {
  save_plot <- .fa_pkg_fun("ggplot2", "ggsave")
  if (is.null(save_plot)) return(list())
  directory <- .fa_module_directory(output, module)
  artifacts <- list()
  for (extension in c("png", "pdf")) {
    path <- file.path(directory, paste0(.safe_name(name), ".", extension))
    ok <- tryCatch({
      save_plot(filename = path, plot = plot, width = width, height = height,
                units = "in", dpi = if (extension == "png") 180 else 300,
                limitsize = FALSE)
      TRUE
    }, error = function(e) FALSE)
    if (ok && file.exists(path)) {
      artifacts[[length(artifacts) + 1L]] <- .fa_artifact_record(
        module, "figure", path, output, paste(label, toupper(extension)), description,
        row_unit = "figure", column_unit = NA_character_, column_dictionary = list()
      )
    }
  }
  artifacts
}

.fa_write_sparse_bundle <- function(x, output, module, name, label, description) {
  if (is.null(x) || !nrow(x) || !ncol(x)) return(list())
  directory <- .fa_module_directory(output, module)
  prefix <- .safe_name(name)
  temporary <- tempfile(fileext = ".mtx")
  on.exit(unlink(temporary), add = TRUE)
  sparse <- if (methods::is(x, "sparseMatrix")) x else Matrix::Matrix(x, sparse = TRUE)
  Matrix::writeMM(sparse, temporary)
  matrix_path <- file.path(directory, paste0(prefix, ".mtx.gz"))
  feature_path <- file.path(directory, paste0(prefix, "_features.tsv.gz"))
  column_path <- file.path(directory, paste0(prefix, "_columns.tsv.gz"))
  .gzip_file(temporary, matrix_path)
  .write_lines_gz(rownames(sparse) %||% seq_len(nrow(sparse)), feature_path)
  .write_lines_gz(colnames(sparse) %||% seq_len(ncol(sparse)), column_path)
  list(
    .fa_artifact_record(
      module, "matrix", matrix_path, output, label, description,
      rows = nrow(sparse), columns = ncol(sparse), row_unit = "feature",
      column_unit = "sample or group", column_dictionary = list(
        matrix = "Matrix Market sparse coordinate matrix; row and column names are stored in companion files."
      )
    ),
    .fa_artifact_record(
      module, "matrix_row_names", feature_path, output, paste(label, "row names"),
      "One feature identifier for every matrix row, in the same order.",
      rows = nrow(sparse), columns = 1L, row_unit = "feature", column_unit = "feature_id",
      column_dictionary = list(feature_id = "Feature identifier for the corresponding matrix row.")
    ),
    .fa_artifact_record(
      module, "matrix_column_names", column_path, output, paste(label, "column names"),
      "One sample or group identifier for every matrix column, in the same order.",
      rows = ncol(sparse), columns = 1L, row_unit = "matrix column", column_unit = "column_id",
      column_dictionary = list(column_id = "Identifier for the corresponding matrix column.")
    )
  )
}

.fa_sanitize_config <- function(x, depth = 0L) {
  if (depth > 4L) return("<nested configuration omitted>")
  if (is.null(x) || is.atomic(x)) {
    if (length(x) > 100L) return(sprintf("<%s values>", length(x)))
    return(x)
  }
  if (is.data.frame(x) || is.matrix(x) || methods::is(x, "Matrix")) {
    return(sprintf("<%s x %s %s>", nrow(x), ncol(x), paste(class(x), collapse = "/")))
  }
  if (isS4(x)) return(sprintf("<%s object>", paste(class(x), collapse = "/")))
  if (is.list(x)) return(lapply(x, .fa_sanitize_config, depth = depth + 1L))
  sprintf("<%s>", paste(class(x), collapse = "/"))
}

.fa_result <- function(object, status = "completed", reason_code = "ok", message = "",
                       engine = NULL, artifacts = list(), details = list()) {
  list(
    object = object, status = status, reason_code = reason_code, message = message,
    engine = engine, artifacts = artifacts, details = details
  )
}

.fa_run_module <- function(id, object, output, config, seed, verbose, fun) {
  cfg <- .fa_module_config(config, id)
  started <- Sys.time()
  warnings <- character()
  requested <- .fa_requested(config, id, default = TRUE)
  if (!requested) {
    result <- .fa_result(object, "skipped", "disabled", "Module disabled by configuration.")
  } else {
    result <- tryCatch(
      withCallingHandlers(
        fun(object, output, cfg, seed, verbose),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) .fa_result(
        object, "failed", "module_error", conditionMessage(e),
        details = list(error_class = class(e))
      )
    )
  }
  if (is.null(result$object)) result$object <- object
  if (is.null(result$artifacts)) result$artifacts <- list()
  finished <- Sys.time()
  engine_package <- if (!is.null(result$engine)) sub("::.*$", "", result$engine) else NULL
  module <- list(
    id = id,
    requested = requested,
    eligible = result$status %in% c("completed", "partial"),
    status = result$status,
    reason_code = result$reason_code,
    message = result$message,
    engine = result$engine,
    engine_version = if (!is.null(engine_package)) .fa_package_version(engine_package) else NA_character_,
    parameters = .fa_sanitize_config(cfg),
    seed = as.integer(seed),
    started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"),
    finished_at = format(finished, "%Y-%m-%dT%H:%M:%S%z"),
    elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
    cells_used = ncol(result$object),
    features_used = nrow(result$object),
    warnings = unique(warnings),
    error = if (identical(result$status, "failed")) result$message else NULL,
    artifact_ids = vapply(result$artifacts, function(x) x$artifact_id %||% "", character(1)),
    details = result$details %||% list()
  )
  .sc_message(verbose, "Full analysis [%s]: %s%s", id, module$status,
              if (nzchar(module$message)) paste0(" - ", module$message) else "")
  list(object = result$object, module = module, artifacts = result$artifacts,
       warnings = unique(warnings))
}

.fa_sparse_group_sum <- function(matrix, groups) {
  groups <- factor(as.character(groups))
  keep <- !is.na(groups)
  if (!all(keep)) {
    matrix <- matrix[, keep, drop = FALSE]
    groups <- droplevels(groups[keep])
  }
  indicator <- Matrix::sparseMatrix(
    i = seq_along(groups), j = as.integer(groups), x = 1,
    dims = c(length(groups), nlevels(groups)),
    dimnames = list(colnames(matrix), levels(groups))
  )
  matrix %*% indicator
}

.fa_scale_columns <- function(matrix, factors) {
  factors <- as.numeric(factors)
  if (length(factors) != ncol(matrix)) stop("Column scaling factors do not match matrix columns.")
  scaled <- matrix %*% Matrix::Diagonal(x = factors)
  dimnames(scaled) <- dimnames(matrix)
  scaled
}

.fa_long_average <- function(matrix, groups, feature_symbols = NULL, scale_factor = 1e4) {
  sums <- .fa_sparse_group_sum(matrix, groups)
  library_sizes <- Matrix::colSums(sums)
  normalized <- .fa_scale_columns(sums, scale_factor / pmax(library_sizes, 1))
  normalized <- log1p(normalized)
  values <- as.vector(as.matrix(normalized))
  output <- data.frame(
    feature = rep(rownames(normalized), times = ncol(normalized)),
    group = rep(colnames(normalized), each = nrow(normalized)),
    average_log_normalized_expression = values,
    stringsAsFactors = FALSE
  )
  if (!is.null(feature_symbols)) {
    output$symbol <- unname(feature_symbols[output$feature])
    output <- output[c("feature", "symbol", "group", "average_log_normalized_expression")]
  }
  output
}

# Quality control -------------------------------------------------------------

.fa_qc_summary <- function(qc, sample, stage, keep = rep(TRUE, nrow(qc))) {
  keep <- as.logical(keep) & !is.na(sample)
  if (!any(keep)) return(data.frame())
  split_rows <- split(which(keep), as.character(sample[keep]))
  rows <- lapply(names(split_rows), function(id) {
    index <- split_rows[[id]]
    data.frame(
      stage = stage,
      sample = id,
      n_cells = length(index),
      median_nCount = stats::median(qc$nCount[index], na.rm = TRUE),
      median_nFeature = stats::median(qc$nFeature[index], na.rm = TRUE),
      median_percent_mt = stats::median(qc$percent_mt[index], na.rm = TRUE),
      median_percent_ribo = stats::median(qc$percent_ribo[index], na.rm = TRUE),
      median_percent_hb = stats::median(qc$percent_hb[index], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.fa_module_qc <- function(object, output, cfg, seed, verbose) {
  assay <- cfg$assay %||% .fa_default_assay(object)
  counts <- .fa_matrix(object, "counts", assay)
  if (is.null(counts)) {
    return(.fa_result(object, "skipped", "counts_missing",
                      "QC requires a counts layer; no counts layer was found."))
  }
  symbols <- .fa_feature_symbols(object, assay)
  orgdb <- .fa_orgdb(list(orgdb = cfg$orgdb), cfg)
  mapping <- .fa_feature_mapping(
    rownames(counts), symbols, orgdb,
    preferred_keytype = cfg$feature_keytype,
    symbol_column = cfg$symbol_column %||% "SYMBOL"
  )
  symbol_values <- as.character(mapping$SYMBOL[match(rownames(counts), mapping$feature)])
  ignore_case <- isTRUE(cfg$pattern_ignore_case %||% TRUE)
  mt_pattern <- cfg$mitochondrial_pattern %||% "^MT-"
  ribo_pattern <- cfg$ribosomal_pattern %||% "^RP[SL]"
  hb_pattern <- cfg$hemoglobin_pattern %||% "^HB[ABDEGMQZ]"
  mt_features <- rownames(counts)[grepl(mt_pattern, symbol_values, ignore.case = ignore_case)]
  ribo_features <- rownames(counts)[grepl(ribo_pattern, symbol_values, ignore.case = ignore_case)]
  hb_features <- rownames(counts)[grepl(hb_pattern, symbol_values, ignore.case = ignore_case)]
  totals <- as.numeric(Matrix::colSums(counts))
  detected <- as.numeric(Matrix::colSums(counts != 0))
  percent_for <- function(features) {
    if (!length(features)) return(rep(NA_real_, ncol(counts)))
    100 * as.numeric(Matrix::colSums(counts[features, , drop = FALSE])) / pmax(totals, 1)
  }
  qc <- data.frame(
    cell = colnames(counts),
    nCount = totals,
    nFeature = detected,
    percent_mt = percent_for(mt_features),
    percent_ribo = percent_for(ribo_features),
    percent_hb = percent_for(hb_features),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  thresholds <- list(
    min_features = as.numeric(cfg$min_features %||% 200),
    max_features = as.numeric(cfg$max_features %||% 7500),
    min_counts = as.numeric(cfg$min_counts %||% 500),
    max_counts = as.numeric(cfg$max_counts %||% Inf),
    max_percent_mt = as.numeric(cfg$max_percent_mt %||% 20),
    max_percent_ribo = as.numeric(cfg$max_percent_ribo %||% Inf),
    max_percent_hb = as.numeric(cfg$max_percent_hb %||% Inf)
  )
  failures <- list(
    low_features = qc$nFeature < thresholds$min_features,
    high_features = qc$nFeature > thresholds$max_features,
    low_counts = qc$nCount < thresholds$min_counts,
    high_counts = qc$nCount > thresholds$max_counts,
    high_mitochondrial_fraction = !is.na(qc$percent_mt) & qc$percent_mt > thresholds$max_percent_mt,
    high_ribosomal_fraction = !is.na(qc$percent_ribo) & qc$percent_ribo > thresholds$max_percent_ribo,
    high_hemoglobin_fraction = !is.na(qc$percent_hb) & qc$percent_hb > thresholds$max_percent_hb
  )
  failure_matrix <- do.call(cbind, failures)
  qc$pass <- rowSums(failure_matrix) == 0L
  qc$reason <- vapply(seq_len(nrow(qc)), function(i) {
    reasons <- names(failures)[failure_matrix[i, ]]
    if (length(reasons)) paste(reasons, collapse = ";") else "pass"
  }, character(1))
  meta <- .seurat_metadata(object)
  sample_column <- .fa_sample_column(object, cfg)
  group_column <- .fa_group_column(object, cfg)
  qc$sample <- if (!is.null(sample_column)) as.character(meta[qc$cell, sample_column]) else "sample1"
  if (!is.null(group_column)) qc$group <- as.character(meta[qc$cell, group_column])

  added <- data.frame(
    .scRDSreport_nCount = qc$nCount,
    .scRDSreport_nFeature = qc$nFeature,
    .scRDSreport_percent_mt = qc$percent_mt,
    .scRDSreport_percent_ribo = qc$percent_ribo,
    .scRDSreport_percent_hb = qc$percent_hb,
    .scRDSreport_qc_pass = qc$pass,
    .scRDSreport_qc_reason = qc$reason,
    row.names = qc$cell,
    check.names = FALSE
  )
  object <- SeuratObject::AddMetaData(object, metadata = added)
  before_summary <- .fa_qc_summary(qc, qc$sample, "before", rep(TRUE, nrow(qc)))
  after_summary <- .fa_qc_summary(qc, qc$sample, "after_pass", qc$pass)
  summary_table <- rbind(before_summary, after_summary)
  threshold_table <- data.frame(
    metric = names(thresholds),
    threshold = unlist(thresholds, use.names = FALSE),
    finite = is.finite(unlist(thresholds, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
  feature_set_table <- data.frame(
    feature_set = c("mitochondrial", "ribosomal", "hemoglobin"),
    matched_features = c(length(mt_features), length(ribo_features), length(hb_features)),
    symbol_source = if (!is.null(orgdb)) "species_orgdb_or_feature_metadata" else if (identical(unname(symbols), names(symbols))) "feature_id" else "feature_metadata",
    matching_pattern = c(mt_pattern, ribo_pattern, hb_pattern),
    stringsAsFactors = FALSE
  )
  artifacts <- list(
    .fa_write_table_artifact(
      qc, output, "qc", "cell_qc", "Per-cell QC metrics",
      "One row per analyzed cell. Metrics are calculated from the sparse counts layer; reason is a compact semicolon-delimited failure code.",
      "cell", list(
        nCount = "Total counts/UMIs for the cell.",
        nFeature = "Number of detected features for the cell.",
        percent_mt = "Percentage of counts assigned to matched mitochondrial genes; NA means no mitochondrial symbols matched.",
        percent_ribo = "Percentage of counts assigned to matched ribosomal protein genes; NA means no ribosomal symbols matched.",
        percent_hb = "Percentage of counts assigned to matched hemoglobin genes; NA means no hemoglobin symbols matched.",
        pass = "TRUE only when the cell passes every finite configured threshold.",
        reason = "Compact failure codes; pass means no configured threshold failed."
      ), units = list(percent_mt = "percent", percent_ribo = "percent", percent_hb = "percent")
    ),
    .fa_write_table_artifact(
      summary_table, output, "qc", "sample_before_after", "QC before/after sample summary",
      "One row per sample and stage. The after_pass stage summarizes cells satisfying every configured QC threshold.",
      "sample-stage", list(stage = "before includes all input cells; after_pass includes cells with pass=TRUE."),
      units = list(median_percent_mt = "percent", median_percent_ribo = "percent", median_percent_hb = "percent")
    ),
    .fa_write_table_artifact(
      threshold_table, output, "qc", "thresholds", "QC thresholds",
      "One row per threshold. Infinite thresholds are recorded but do not exclude cells.",
      "threshold", list(metric = "QC threshold name.", threshold = "Configured cutoff.", finite = "Whether this cutoff actively filters cells.")
    ),
    .fa_write_table_artifact(
      feature_set_table, output, "qc", "feature_set_matching", "QC feature-set matching",
      "One row per QC gene family; zero matches explains an NA percentage column.",
      "gene family", list(matched_features = "Number of matrix rows matched to this family.")
    )
  )

  max_plot_cells <- as.integer(cfg$max_plot_cells %||% 50000L)
  plot_rows <- seq_len(nrow(qc))
  if (length(plot_rows) > max_plot_cells) {
    set.seed(as.integer(seed))
    plot_rows <- sort(sample(plot_rows, max_plot_cells))
  }
  plot_qc <- qc[plot_rows, , drop = FALSE]
  long <- rbind(
    data.frame(sample = plot_qc$sample, metric = "nCount", value = plot_qc$nCount),
    data.frame(sample = plot_qc$sample, metric = "nFeature", value = plot_qc$nFeature),
    data.frame(sample = plot_qc$sample, metric = "percent_mt", value = plot_qc$percent_mt),
    data.frame(sample = plot_qc$sample, metric = "percent_ribo", value = plot_qc$percent_ribo),
    data.frame(sample = plot_qc$sample, metric = "percent_hb", value = plot_qc$percent_hb)
  )
  long <- long[is.finite(long$value), , drop = FALSE]
  if (nrow(long)) {
    violin <- ggplot2::ggplot(long, ggplot2::aes(x = sample, y = value, fill = sample)) +
      ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
      ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.25) +
      ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 2) +
      ggplot2::scale_fill_manual(values = .fa_palette(length(unique(long$sample)))) +
      ggplot2::labs(x = "Sample", y = "QC value", title = "Per-sample QC distributions") +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    artifacts <- c(artifacts, .fa_write_plot_artifacts(
      violin, output, "qc", "qc_distributions", "QC distributions",
      "Violin and embedded box plots show the cell-level distributions for each sample; plotting may use a deterministic cell subset."
    ))
  }
  scatter <- ggplot2::ggplot(plot_qc, ggplot2::aes(x = nCount, y = nFeature, color = pass)) +
    ggplot2::geom_point(size = 0.45, alpha = 0.55) +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#D55E00")) +
    ggplot2::labs(x = "Total counts", y = "Detected features", color = "QC pass",
                  title = "Cell complexity and QC status") +
    ggplot2::theme_bw(base_size = 12)
  artifacts <- c(artifacts, .fa_write_plot_artifacts(
    scatter, output, "qc", "counts_features", "Counts versus features",
    "Each point is one cell; color indicates whether all configured QC thresholds were passed."
  ))

  status <- "completed"
  reason <- "qc_metrics_exported"
  message <- sprintf("QC metrics were calculated for %s cells; %s passed.", nrow(qc), sum(qc$pass))
  if (isTRUE(cfg$filter)) {
    if (sum(qc$pass) < 10L) {
      status <- "partial"
      reason <- "filter_would_retain_too_few_cells"
      message <- paste0(message, " Filtering was not applied because fewer than 10 cells would remain.")
    } else {
      object <- object[, qc$cell[qc$pass], drop = FALSE]
      reason <- "qc_filter_applied"
      message <- paste0(message, " The analysis object was subset to passing cells by explicit configuration.")
    }
  }
  .fa_result(
    object, status, reason, message, "Matrix sparse QC", artifacts,
    list(assay = assay, matched_feature_sets = feature_set_table, filtering_applied = isTRUE(cfg$filter) && sum(qc$pass) >= 10L)
  )
}

# Dimensional reduction -------------------------------------------------------

.fa_module_reduction <- function(object, output, cfg, seed, verbose) {
  artifacts <- list()
  assay <- cfg$assay %||% .fa_default_assay(object)
  variable_features <- tryCatch(SeuratObject::VariableFeatures(object[[assay]]), error = function(e) character())
  if (length(variable_features)) {
    vf <- data.frame(feature = variable_features, rank = seq_along(variable_features), stringsAsFactors = FALSE)
    artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
      vf, output, "reduction", "variable_features", "Highly variable features",
      "One row per variable feature, in the order stored in the Seurat assay.", "feature",
      list(rank = "One-based rank in SeuratObject::VariableFeatures order.")
    )
  }
  reductions <- .seurat_reductions(object)
  if (!length(reductions) && !length(variable_features)) {
    return(.fa_result(object, "skipped", "no_reduction",
                      "The object contains no dimensional reduction or variable-feature result."))
  }
  meta <- .seurat_metadata(object)
  color_column <- cfg$color_column %||% .fa_annotation_column(object, cfg) %||%
    .fa_cluster_column(object, cfg) %||% .fa_sample_column(object, cfg)
  max_plot_cells <- as.integer(cfg$max_plot_cells %||% 100000L)
  for (reduction in reductions) {
    embedding <- .embedding_table(object, reduction)
    if (is.null(embedding) || ncol(embedding) < 2L) next
    artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
      embedding, output, "reduction", paste0("embedding_", reduction), paste("Embedding", reduction),
      "One row per cell. The first column is the exact cell name; remaining columns are stored reduction coordinates.",
      "cell", list(cell = "Exact cell name corresponding to a column of the expression matrix.")
    )
    if (ncol(embedding) < 3L) next
    reduction_object <- tryCatch(object[[reduction]], error = function(e) NULL)
    stdev <- tryCatch(SeuratObject::Stdev(reduction_object), error = function(e) NULL)
    if (!is.null(stdev) && length(stdev)) {
      variance <- stdev^2
      pca_variance <- data.frame(
        reduction = reduction,
        component = seq_along(stdev),
        standard_deviation = as.numeric(stdev),
        variance = as.numeric(variance),
        variance_fraction = as.numeric(variance / sum(variance)),
        cumulative_variance_fraction = as.numeric(cumsum(variance) / sum(variance)),
        stringsAsFactors = FALSE
      )
      artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
        pca_variance, output, "reduction", paste0("variance_", reduction), paste("Variance", reduction),
        "One row per stored linear-reduction component.", "component",
        list(variance_fraction = "Fraction of total stored component variance explained by this component.",
             cumulative_variance_fraction = "Cumulative fraction through this component."),
        units = list(variance_fraction = "fraction", cumulative_variance_fraction = "fraction")
      )
    }
    coordinates <- embedding[, 2:3, drop = FALSE]
    names(coordinates) <- c("dimension_1", "dimension_2")
    coordinates$cell <- embedding$cell
    if (!is.null(color_column) && color_column %in% names(meta)) {
      coordinates$color_value <- as.character(meta[coordinates$cell, color_column])
    } else {
      coordinates$color_value <- "all_cells"
      color_column <- "all_cells"
    }
    if (nrow(coordinates) > max_plot_cells) {
      set.seed(as.integer(seed))
      coordinates <- coordinates[sort(sample(seq_len(nrow(coordinates)), max_plot_cells)), , drop = FALSE]
    }
    coordinates$color_value <- factor(coordinates$color_value)
    plot <- ggplot2::ggplot(coordinates, ggplot2::aes(x = dimension_1, y = dimension_2, color = color_value)) +
      ggplot2::geom_point(size = 0.42, alpha = 0.78) +
      ggplot2::scale_color_manual(values = .fa_palette(nlevels(coordinates$color_value))) +
      ggplot2::coord_equal() +
      ggplot2::labs(x = names(embedding)[2L], y = names(embedding)[3L], color = color_column,
                    title = paste("Embedding:", reduction)) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    artifacts <- c(artifacts, .fa_write_plot_artifacts(
      plot, output, "reduction", paste0("embedding_", reduction), paste("Embedding", reduction),
      paste0("Two-dimensional view of ", reduction, " colored by ", color_column,
             "; large objects use a deterministic plotting subset only.")
    ))
  }
  .fa_result(
    object, if (length(artifacts)) "completed" else "partial",
    if (length(artifacts)) "reductions_exported" else "reduction_export_failed",
    sprintf("Exported %s stored reductions and the available variable-feature/PCA variance results.", length(reductions)),
    "SeuratObject reductions", artifacts,
    list(reductions = reductions, color_column = color_column)
  )
}

# Clustering ------------------------------------------------------------------

.fa_module_cluster <- function(object, output, cfg, seed, verbose) {
  cluster_column <- .fa_cluster_column(object, cfg)
  if (is.null(cluster_column)) {
    return(.fa_result(object, "skipped", "cluster_missing",
                      "No multi-level cluster identity is present; no cluster was invented."))
  }
  if (isTRUE(cfg$run_markers)) {
    existing_misc <- .slot_or_null(object, "misc") %||% list()
    existing_names <- names(.collect_result_tables(existing_misc))
    if (!any(grepl("marker|detest|differential", existing_names, ignore.case = TRUE))) {
      marker_result <- .run_cluster_markers(
        object, marker_args = cfg$marker_args %||% list(),
        cluster_prefix = cfg$cluster_prefix %||% "scRDSreport", verbose = verbose
      )
      object <- marker_result$object
    }
  }
  meta <- .seurat_metadata(object)
  assignments <- data.frame(
    cell = rownames(meta),
    cluster = as.character(meta[[cluster_column]]),
    stringsAsFactors = FALSE
  )
  sample_column <- .fa_sample_column(object, cfg)
  if (!is.null(sample_column)) assignments$sample <- as.character(meta[[sample_column]])
  sizes <- as.data.frame(table(assignments$cluster), stringsAsFactors = FALSE)
  names(sizes) <- c("cluster", "n_cells")
  sizes$fraction <- sizes$n_cells / sum(sizes$n_cells)
  artifacts <- list(
    .fa_write_table_artifact(
      assignments, output, "cluster", "cluster_assignments", "Cluster assignments",
      "One row per cell. Cluster values are copied from the selected existing metadata column.",
      "cell", list(cluster = paste0("Value of metadata column '", cluster_column, "'."))
    ),
    .fa_write_table_artifact(
      sizes, output, "cluster", "cluster_sizes", "Cluster sizes",
      "One row per cluster, with cell count and fraction of analyzed cells.",
      "cluster", list(fraction = "Cluster cell count divided by the total analyzed cell count."),
      units = list(fraction = "fraction")
    )
  )
  misc_tables <- .collect_result_tables(.slot_or_null(object, "misc") %||% list())
  marker_tables <- misc_tables[grepl("marker|detest|differential", names(misc_tables), ignore.case = TRUE)]
  for (name in names(marker_tables)) {
    table <- marker_tables[[name]]
    if (!is.data.frame(table) || !nrow(table)) next
    artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
      table, output, "cluster", paste0("markers_", name), paste("Cluster marker result", name),
      "Rows are marker-test records retained in the Seurat object; exact fields are preserved.",
      "marker-test record"
    )
  }
  size_plot <- ggplot2::ggplot(sizes, ggplot2::aes(x = stats::reorder(cluster, n_cells), y = n_cells, fill = cluster)) +
    ggplot2::geom_col(width = 0.78) +
    ggplot2::scale_fill_manual(values = .fa_palette(nrow(sizes))) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Cluster", y = "Cells", title = "Cluster sizes") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "none")
  artifacts <- c(artifacts, .fa_write_plot_artifacts(
    size_plot, output, "cluster", "cluster_sizes", "Cluster sizes",
    "Bar length is the number of cells in each existing cluster."
  ))
  .fa_result(
    object, "completed", "clusters_exported",
    sprintf("Cluster assignments were exported from '%s' (%s clusters).", cluster_column, nrow(sizes)),
    "Seurat cluster metadata", artifacts,
    list(cluster_column = cluster_column, n_clusters = nrow(sizes))
  )
}

# Cell annotation and composition --------------------------------------------

.fa_reference_labels <- function(reference, cfg) {
  if (!is.null(cfg$reference_labels)) return(as.character(cfg$reference_labels))
  col_data_fun <- .fa_pkg_fun("SummarizedExperiment", "colData")
  if (is.null(col_data_fun)) return(NULL)
  col_data <- tryCatch(as.data.frame(col_data_fun(reference)), error = function(e) NULL)
  if (is.null(col_data)) return(NULL)
  candidates <- unique(c(cfg$reference_label_column, "label.main", "label.fine", "label", "cell_type"))
  candidates <- candidates[!is.na(candidates) & candidates %in% names(col_data)]
  if (length(candidates)) as.character(col_data[[candidates[[1L]]]]) else NULL
}

.fa_load_annotation_reference <- function(species, cfg) {
  reference <- cfg$reference %||% cfg$auto_annotation_reference
  if (is.character(reference) && length(reference) == 1L && file.exists(reference)) {
    return(readRDS(reference))
  }
  if (is.character(reference) && length(reference) == 1L && grepl("::", reference, fixed = TRUE)) {
    fields <- strsplit(reference, "::", fixed = TRUE)[[1L]]
    loader <- if (length(fields) == 2L) .fa_pkg_fun(fields[[1L]], fields[[2L]]) else NULL
    if (is.null(loader)) return(NULL)
    return(tryCatch(loader(), error = function(e) NULL))
  }
  if (!is.null(reference) && !is.character(reference)) return(reference)
  if (!is.null(reference)) return(NULL)
  if (isFALSE(cfg$allow_reference_download)) return(NULL)
  function_name <- switch(
    species,
    mouse = "MouseRNAseqData",
    human = "HumanPrimaryCellAtlasData",
    NULL
  )
  if (is.null(function_name)) return(NULL)
  loader <- .fa_pkg_fun("celldex", function_name)
  if (is.null(loader)) return(NULL)
  tryCatch(loader(), error = function(e) NULL)
}

.fa_apply_manual_annotation <- function(object, cfg, cluster_column) {
  mapping <- cfg$mapping
  marker_scores <- NULL
  if (is.null(mapping) && !is.null(cfg$manual_markers) && !is.list(cfg$manual_markers) &&
      !is.null(names(cfg$manual_markers))) {
    mapping <- cfg$manual_markers
  }
  if (is.null(mapping) && is.list(cfg$manual_markers) && length(cfg$manual_markers) &&
      !is.null(names(cfg$manual_markers)) && !is.null(cluster_column)) {
    assay <- cfg$assay %||% .fa_default_assay(object)
    expression <- .fa_matrix(object, "data", assay)
    if (is.null(expression)) expression <- .fa_matrix(object, "counts", assay)
    meta <- .seurat_metadata(object)
    if (!is.null(expression)) {
      symbols <- .fa_feature_symbols(object, assay)
      clusters <- as.character(meta[colnames(expression), cluster_column])
      sums <- .fa_sparse_group_sum(expression, clusters)
      group_sizes <- as.numeric(table(factor(clusters, levels = colnames(sums))))
      averages <- .fa_scale_columns(sums, 1 / pmax(group_sizes, 1))
      score_rows <- list()
      for (label in names(cfg$manual_markers)) {
        features <- .fa_match_features(as.character(cfg$manual_markers[[label]]), symbols)
        if (length(features) < as.integer(cfg$min_manual_marker_genes %||% 1L)) next
        score_rows[[label]] <- data.frame(
          annotation = label,
          cluster = colnames(averages),
          marker_score = as.numeric(Matrix::colMeans(averages[features, , drop = FALSE])),
          n_markers_matched = length(features),
          stringsAsFactors = FALSE
        )
      }
      if (length(score_rows)) {
        marker_scores <- do.call(rbind, score_rows)
        split_scores <- split(marker_scores, marker_scores$cluster)
        winners <- lapply(split_scores, function(x) x[which.max(x$marker_score), , drop = FALSE])
        winners <- do.call(rbind, winners)
        mapping <- stats::setNames(as.character(winners$annotation), as.character(winners$cluster))
      }
    }
  }
  if (is.null(mapping) || is.null(cluster_column)) {
    return(list(object = object, column = NULL, scores = marker_scores,
                message = "Manual annotation requires both a cluster column and a named mapping."))
  }
  if (is.data.frame(mapping)) {
    cluster_name <- intersect(c("cluster", "cluster_id", cluster_column), names(mapping))
    label_name <- intersect(c("annotation", "celltype", "cell_type", "label"), names(mapping))
    if (!length(cluster_name) || !length(label_name)) {
      return(list(object = object, column = NULL, scores = marker_scores,
                  message = "Manual mapping data frame needs cluster and annotation columns."))
    }
    values <- as.character(mapping[[label_name[[1L]]]])
    names(values) <- as.character(mapping[[cluster_name[[1L]]]])
    mapping <- values
  }
  if (is.null(names(mapping)) || any(!nzchar(names(mapping)))) {
    return(list(object = object, column = NULL, scores = marker_scores,
                message = "Manual annotation mapping must be named by cluster."))
  }
  meta <- .seurat_metadata(object)
  output_column <- cfg$output_column %||% ".scRDSreport_celltype_manual"
  if (output_column %in% names(meta)) {
    return(list(object = object, column = NULL, scores = marker_scores,
                message = paste0("Refusing to overwrite existing metadata column '", output_column, "'.")))
  }
  labels <- unname(as.character(mapping)[match(as.character(meta[[cluster_column]]), names(mapping))])
  if (!any(!is.na(labels) & nzchar(labels))) {
    return(list(
      object = object, column = NULL, scores = marker_scores,
      message = "The explicit manual mapping did not match any analyzed cluster; no all-missing annotation column was added."
    ))
  }
  annotation <- data.frame(labels, row.names = rownames(meta), check.names = FALSE)
  names(annotation) <- output_column
  object <- SeuratObject::AddMetaData(object, annotation)
  list(
    object = object, column = output_column, scores = marker_scores,
    message = sprintf("Manual mapping annotated %s of %s cells in a new metadata column.", sum(!is.na(labels)), length(labels))
  )
}

.fa_apply_singler_annotation <- function(object, species, cfg, cluster_column) {
  single_r <- .fa_pkg_fun("SingleR", "SingleR")
  if (is.null(single_r)) {
    return(list(object = object, column = NULL, predictions = NULL,
                reason = "dependency_missing", message = "SingleR is not installed."))
  }
  reference <- .fa_load_annotation_reference(species, cfg)
  if (is.null(reference)) {
    return(list(
      object = object, column = NULL, predictions = NULL, reason = "reference_missing",
      message = "Automatic annotation was explicitly requested, but no reference object was supplied. Set reference/reference_labels, or explicitly allow a celldex download for human or mouse."
    ))
  }
  labels <- .fa_reference_labels(reference, cfg)
  if (is.null(labels) || length(labels) == 0L) {
    return(list(object = object, column = NULL, predictions = NULL,
                reason = "reference_labels_missing", message = "The annotation reference has no usable label vector."))
  }
  assay <- cfg$assay %||% .fa_default_assay(object)
  test <- .fa_matrix(object, "data", assay)
  if (is.null(test)) {
    normalize <- .fa_pkg_fun("Seurat", "NormalizeData")
    if (is.null(normalize)) {
      return(list(object = object, column = NULL, predictions = NULL,
                  reason = "normalized_data_missing", message = "SingleR requires normalized expression data."))
    }
    object <- normalize(object, assay = assay, verbose = FALSE)
    test <- .fa_matrix(object, "data", assay)
  }
  symbols_from_object <- .fa_feature_symbols(object, assay)
  resources <- .fa_species_resources(species, cfg)
  orgdb <- .fa_orgdb(resources, cfg)
  mapping <- .fa_feature_mapping(
    rownames(test), symbols_from_object, orgdb,
    preferred_keytype = cfg$feature_keytype %||% resources$feature_keytype,
    symbol_column = cfg$symbol_column %||% resources$symbol_column %||% "SYMBOL"
  )
  symbols <- as.character(mapping$SYMBOL[match(rownames(test), mapping$feature)])
  valid <- !is.na(symbols) & nzchar(symbols) & !duplicated(toupper(symbols))
  test <- test[valid, , drop = FALSE]
  rownames(test) <- symbols[valid]
  if (nrow(test) < 100L) {
    return(list(object = object, column = NULL, predictions = NULL,
                reason = "insufficient_feature_overlap", message = "Fewer than 100 unique feature symbols were available for reference annotation."))
  }
  meta <- .seurat_metadata(object)
  clusters <- if (!is.null(cluster_column)) as.character(meta[[cluster_column]]) else NULL
  if (is.null(clusters) && ncol(test) > as.integer(cfg$max_cell_level_cells %||% 50000L)) {
    return(list(object = object, column = NULL, predictions = NULL,
                reason = "cluster_required_for_large_object",
                message = "Cell-level SingleR was not run on a large object without clusters; supply a cluster column."))
  }
  prediction_mode <- if (is.null(clusters)) "cell_level" else "cluster_level"
  predictions <- tryCatch(
    single_r(test = test, ref = reference, labels = labels, clusters = clusters),
    error = function(e) e
  )
  if (inherits(predictions, "error") && !is.null(clusters) &&
      grepl("scrapper", conditionMessage(predictions), ignore.case = TRUE) &&
      ncol(test) <= as.integer(cfg$max_cell_level_cells %||% 50000L) &&
      !isFALSE(cfg$allow_cell_level_fallback)) {
    clusters <- NULL
    prediction_mode <- "cell_level_fallback_missing_cluster_aggregation_dependency"
    predictions <- tryCatch(
      single_r(test = test, ref = reference, labels = labels, clusters = NULL),
      error = function(e) e
    )
  }
  if (inherits(predictions, "error")) {
    return(list(object = object, column = NULL, predictions = NULL,
                reason = "singler_failed", message = conditionMessage(predictions)))
  }
  prediction_table <- as.data.frame(predictions)
  prediction_table$prediction_id <- rownames(prediction_table)
  use_pruned <- !isFALSE(cfg$use_pruned_labels)
  selected <- if (use_pruned && "pruned.labels" %in% names(prediction_table)) {
    as.character(prediction_table$pruned.labels)
  } else if ("labels" %in% names(prediction_table)) {
    as.character(prediction_table$labels)
  } else {
    rep(NA_character_, nrow(prediction_table))
  }
  selected[is.na(selected) | !nzchar(selected)] <- NA_character_
  prediction_table$selected_label <- selected
  prediction_table$selected_label_source <- if (use_pruned && "pruned.labels" %in% names(prediction_table)) {
    "pruned.labels"
  } else {
    "labels"
  }
  if (!is.null(clusters)) {
    cell_labels <- selected[match(clusters, rownames(prediction_table))]
  } else {
    cell_labels <- selected[match(rownames(meta), rownames(prediction_table))]
  }
  reference_value <- cfg$reference %||% cfg$auto_annotation_reference
  reference_label <- if (is.null(reference_value)) {
    "registered celldex reference"
  } else if (is.character(reference_value)) {
    paste(reference_value, collapse = ", ")
  } else {
    paste0("<", paste(class(reference_value), collapse = "/"), " object>")
  }
  confident_cells <- sum(!is.na(cell_labels) & nzchar(cell_labels))
  provenance <- list(
    species = species,
    reference = reference_label,
    label_source = if (use_pruned && "pruned.labels" %in% names(prediction_table)) "pruned.labels" else "labels",
    prediction_mode = prediction_mode,
    confident_cells = confident_cells,
    total_cells = length(cell_labels)
  )
  if (confident_cells == 0L) {
    return(list(
      object = object, column = NULL, predictions = prediction_table,
      reason = "no_confident_labels",
      message = "SingleR returned no confident labels after pruning; predictions were exported, but no all-missing annotation column was added.",
      provenance = provenance
    ))
  }
  output_column <- cfg$output_column %||% ".scRDSreport_celltype_SingleR"
  if (output_column %in% names(meta)) {
    return(list(object = object, column = NULL, predictions = prediction_table,
                reason = "output_column_exists",
                message = paste0("Refusing to overwrite existing metadata column '", output_column, "'.")))
  }
  annotation <- data.frame(cell_labels, row.names = rownames(meta), check.names = FALSE)
  names(annotation) <- output_column
  object <- SeuratObject::AddMetaData(object, annotation)
  list(
    object = object, column = output_column, predictions = prediction_table,
    reason = "singler_completed",
    message = sprintf(
      "SingleR added a new species-matched reference annotation column (%s): %s of %s cells received a confident label; pruned/uncertain labels remain missing.",
      prediction_mode, confident_cells, length(cell_labels)
    ),
    provenance = provenance
  )
}

.fa_module_celltype <- function(object, output, cfg, seed, verbose, species) {
  requested_mode <- tolower(as.character(cfg$mode %||% "auto_if_missing"))
  if (!requested_mode %in% c("auto_if_missing", "preserve", "manual", "auto")) {
    requested_mode <- "auto_if_missing"
  }
  original_annotation <- .fa_annotation_column(object, cfg)
  mode <- if (identical(requested_mode, "auto_if_missing")) {
    if (is.null(original_annotation)) "auto" else "preserve"
  } else {
    requested_mode
  }
  cluster_column <- .fa_cluster_column(object, cfg)
  annotation_column <- original_annotation
  annotation_engine <- "RDS metadata"
  annotation_message <- "Existing annotations were preserved without modification."
  prediction_table <- NULL
  annotation_provenance <- NULL
  annotation_issue <- NULL
  if (identical(mode, "manual")) {
    manual <- .fa_apply_manual_annotation(object, cfg, cluster_column)
    object <- manual$object
    annotation_column <- manual$column
    annotation_engine <- "explicit manual mapping"
    annotation_message <- manual$message
    prediction_table <- manual$scores
    if (is.null(annotation_column)) annotation_issue <- "manual_mapping_invalid"
  } else if (identical(mode, "auto")) {
    automatic <- .fa_apply_singler_annotation(object, species, cfg, cluster_column)
    object <- automatic$object
    annotation_column <- automatic$column
    prediction_table <- automatic$predictions
    annotation_provenance <- automatic$provenance %||% NULL
    annotation_engine <- "SingleR"
    annotation_message <- automatic$message
    if (is.null(annotation_column)) annotation_issue <- automatic$reason
  }
  if (!is.null(annotation_column) && !identical(mode, "preserve")) {
    object <- .fa_store_analysis_misc(object, "active_annotation_column", annotation_column)
  }
  meta <- .seurat_metadata(object)
  artifacts <- list()
  if (!is.null(annotation_column) && annotation_column %in% names(meta)) {
    annotation_table <- data.frame(
      cell = rownames(meta),
      annotation = as.character(meta[[annotation_column]]),
      annotation_source_column = annotation_column,
      stringsAsFactors = FALSE
    )
    if (!is.null(cluster_column)) annotation_table$cluster <- as.character(meta[[cluster_column]])
    artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
      annotation_table, output, "celltype", "cell_annotations", "Cell annotations",
      "One row per cell. Original RDS annotations are copied exactly; manual/SingleR columns exist only when explicitly requested and never overwrite an existing column.",
      "cell", list(annotation = paste0("Value of metadata column '", annotation_column, "'."),
                   annotation_source_column = "Metadata column from which annotation was read.")
    )
  }
  if (!is.null(prediction_table) && nrow(prediction_table)) {
    artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
      prediction_table, output, "celltype",
      if (identical(mode, "manual")) "manual_marker_scores" else "singler_predictions",
      if (identical(mode, "manual")) "Manual marker scores" else "SingleR predictions",
      if (identical(mode, "manual")) {
        "One row per supplied cell-type marker set and cluster. The maximum mean marker score assigns the explicitly requested manual label; these scores are not a reference-atlas annotation."
      } else {
        "One row per tested cell or cluster. Labels and scores are retained exactly from SingleR."
      },
      if (!is.null(cluster_column)) "cluster" else "cell"
    )
  }

  label_column <- annotation_column %||% cluster_column
  label_kind <- if (!is.null(annotation_column)) "cell annotation" else "cluster"
  sample_column <- .fa_sample_column(object, cfg)
  group_column <- .fa_group_column(object, cfg)
  if (!is.null(label_column) && !is.null(sample_column)) {
    sample <- as.character(meta[[sample_column]])
    label <- as.character(meta[[label_column]])
    valid <- !is.na(sample) & nzchar(sample) & !is.na(label) & nzchar(label)
    counts <- as.data.frame(table(sample = sample[valid], label = label[valid]), stringsAsFactors = FALSE)
    counts <- counts[counts$Freq > 0L, , drop = FALSE]
    names(counts)[names(counts) == "Freq"] <- "n_cells"
    sample_totals <- stats::setNames(as.numeric(tapply(counts$n_cells, counts$sample, sum)),
                                     names(tapply(counts$n_cells, counts$sample, sum)))
    counts$fraction <- counts$n_cells / sample_totals[counts$sample]
    if (!is.null(group_column)) {
      sample_group <- tapply(as.character(meta[[group_column]]), sample, function(x) {
        values <- unique(x[!is.na(x) & nzchar(x)])
        if (length(values)) values[[1L]] else NA_character_
      })
      counts$group <- unname(sample_group[counts$sample])
    }
    artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
      counts, output, "celltype", "sample_composition", paste("Sample", label_kind, "composition"),
      paste0("One row per observed sample-by-", label_kind,
             " combination. Fractions sum to one within each sample after excluding missing labels."),
      paste("sample-by", label_kind),
      list(label = paste0("Value of '", label_column, "'."),
           n_cells = "Number of cells in this sample and label.",
           fraction = "n_cells divided by all labeled cells in the sample."),
      units = list(fraction = "fraction")
    )
    counts$label <- factor(counts$label)
    composition_plot <- ggplot2::ggplot(counts, ggplot2::aes(x = sample, y = fraction, fill = label)) +
      ggplot2::geom_col(width = 0.82) +
      ggplot2::scale_fill_manual(values = .fa_palette(nlevels(counts$label))) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
      ggplot2::labs(x = "Sample", y = "Cell fraction", fill = label_kind,
                    title = paste("Per-sample", label_kind, "composition")) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    artifacts <- c(artifacts, .fa_write_plot_artifacts(
      composition_plot, output, "celltype", "sample_composition", paste("Sample", label_kind, "composition"),
      "Stacked bars show cell fractions within each sample. This is descriptive composition, not an inferential replicate-level test."
    ))
  }
  reduction <- .seurat_reductions(object)
  reduction <- c(reduction[grepl("umap", reduction, ignore.case = TRUE)], reduction)
  reduction <- unique(reduction)
  if (!is.null(label_column) && length(reduction)) {
    embedding <- .embedding_table(object, reduction[[1L]])
    if (!is.null(embedding) && ncol(embedding) >= 3L) {
      plotted <- data.frame(
        dimension_1 = embedding[[2L]], dimension_2 = embedding[[3L]],
        label = factor(as.character(meta[embedding$cell, label_column])),
        stringsAsFactors = FALSE
      )
      max_plot <- as.integer(cfg$max_plot_cells %||% 100000L)
      if (nrow(plotted) > max_plot) {
        set.seed(as.integer(seed))
        plotted <- plotted[sort(sample(seq_len(nrow(plotted)), max_plot)), , drop = FALSE]
      }
      annotation_plot <- ggplot2::ggplot(plotted, ggplot2::aes(x = dimension_1, y = dimension_2, color = label)) +
        ggplot2::geom_point(size = 0.42, alpha = 0.78) +
        ggplot2::scale_color_manual(values = .fa_palette(nlevels(plotted$label))) +
        ggplot2::coord_equal() +
        ggplot2::labs(x = names(embedding)[2L], y = names(embedding)[3L], color = label_kind,
                      title = paste("Embedding colored by", label_kind)) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(panel.grid = ggplot2::element_blank())
      embedding_stem <- if (!is.null(annotation_column)) "annotation_embedding" else "cluster_context_embedding"
      embedding_label <- if (!is.null(annotation_column)) {
        "Cell annotation embedding"
      } else {
        "Cluster context embedding (not cell annotation)"
      }
      embedding_description <- if (!is.null(annotation_column)) {
        paste0("Stored ", reduction[[1L]], " coordinates colored by annotation column ", label_column, ".")
      } else {
        paste0(
          "Stored ", reduction[[1L]], " coordinates colored by cluster column ", label_column,
          ". This is clustering context only and must not be interpreted as a cell-type annotation."
        )
      }
      artifacts <- c(artifacts, .fa_write_plot_artifacts(
        annotation_plot, output, "celltype", embedding_stem, embedding_label,
        embedding_description
      ))
    }
  }
  if (!is.null(annotation_issue)) {
    return(.fa_result(
      object, "needs_input", annotation_issue, annotation_message, annotation_engine,
      artifacts, list(mode = mode, original_annotation_column = original_annotation,
                      requested_mode = requested_mode, annotation_column = annotation_column,
                      cluster_column = cluster_column,
                      annotation_provenance = annotation_provenance)
    ))
  }
  if (is.null(annotation_column)) {
    return(.fa_result(
      object, if (length(artifacts)) "partial" else "skipped", "annotation_missing",
      "No original cell annotation was present. Cluster composition may be shown, but no cell type was invented.",
      "RDS metadata", artifacts,
      list(mode = mode, requested_mode = requested_mode,
           annotation_column = NULL, cluster_column = cluster_column)
    ))
  }
  .fa_result(
    object, "completed", if (identical(mode, "preserve")) "annotation_preserved" else "annotation_reference_added",
    annotation_message, annotation_engine, artifacts,
    list(mode = mode, requested_mode = requested_mode,
         original_annotation_column = original_annotation,
         annotation_column = annotation_column, cluster_column = cluster_column,
         annotation_provenance = annotation_provenance)
  )
}

# Differential expression ----------------------------------------------------

.fa_store_analysis_misc <- function(object, name, value) {
  if (!.is_seurat(object) || !"misc" %in% methods::slotNames(object)) return(object)
  misc <- .slot_or_null(object, "misc") %||% list()
  container <- misc$scRDSreport_full_analysis %||% list()
  container[[name]] <- value
  misc$scRDSreport_full_analysis <- container
  methods::slot(object, "misc") <- misc
  object
}

.fa_get_analysis_misc <- function(object, name) {
  misc <- .slot_or_null(object, "misc") %||% list()
  container <- misc$scRDSreport_full_analysis %||% list()
  container[[name]]
}

.fa_contrasts <- function(groups, cfg) {
  groups <- sort(unique(as.character(groups[!is.na(groups) & nzchar(groups)])))
  supplied <- cfg$contrasts
  output <- data.frame(name = character(), group_a = character(), group_b = character(), stringsAsFactors = FALSE)
  if (is.data.frame(supplied) && all(c("group_a", "group_b") %in% names(supplied))) {
    output <- data.frame(
      name = if ("name" %in% names(supplied)) as.character(supplied$name) else paste(supplied$group_b, "vs", supplied$group_a),
      group_a = as.character(supplied$group_a), group_b = as.character(supplied$group_b),
      stringsAsFactors = FALSE
    )
  } else if (is.list(supplied) && length(supplied)) {
    rows <- lapply(seq_along(supplied), function(i) {
      contrast <- as.character(supplied[[i]])
      if (length(contrast) < 2L) return(NULL)
      data.frame(
        name = names(supplied)[i] %||% paste(contrast[[2L]], "vs", contrast[[1L]]),
        group_a = contrast[[1L]], group_b = contrast[[2L]], stringsAsFactors = FALSE
      )
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows)) output <- do.call(rbind, rows)
  } else if (length(groups) >= 2L) {
    pairs <- utils::combn(groups, 2L)
    output <- data.frame(
      name = paste(pairs[2L, ], "vs", pairs[1L, ]),
      group_a = pairs[1L, ], group_b = pairs[2L, ], stringsAsFactors = FALSE
    )
  }
  output <- output[output$group_a %in% groups & output$group_b %in% groups & output$group_a != output$group_b, , drop = FALSE]
  max_contrasts <- as.integer(cfg$max_contrasts %||% 20L)
  if (nrow(output) > max_contrasts) output <- output[seq_len(max_contrasts), , drop = FALSE]
  rownames(output) <- NULL
  output
}

.fa_edgeR_contrast <- function(pseudobulk, sample_groups, group_a, group_b, cfg) {
  functions <- lapply(
    c("DGEList", "calcNormFactors", "filterByExpr", "estimateDisp", "glmQLFit", "glmQLFTest", "topTags"),
    function(name) .fa_pkg_fun("edgeR", name)
  )
  names(functions) <- c("DGEList", "calcNormFactors", "filterByExpr", "estimateDisp", "glmQLFit", "glmQLFTest", "topTags")
  if (any(vapply(functions, is.null, logical(1)))) stop("The installed edgeR does not expose the required quasi-likelihood functions.")
  selected <- sample_groups %in% c(group_a, group_b)
  counts <- pseudobulk[, selected, drop = FALSE]
  group <- droplevels(factor(sample_groups[selected], levels = c(group_a, group_b)))
  if (any(table(group) < 2L)) stop("Both contrast groups need at least two biological samples for inferential pseudobulk DE.")
  design <- stats::model.matrix(~0 + group)
  colnames(design) <- levels(group)
  y <- functions$DGEList(counts = counts, group = group)
  keep <- functions$filterByExpr(y, design = design, min.count = cfg$min_count %||% 10)
  if (sum(keep) < 10L) stop("Fewer than 10 features passed edgeR filtering for this contrast.")
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- functions$calcNormFactors(y)
  y <- functions$estimateDisp(y, design, robust = TRUE)
  fit <- functions$glmQLFit(y, design, robust = TRUE)
  contrast <- rep(0, ncol(design))
  contrast[match(group_b, colnames(design))] <- 1
  contrast[match(group_a, colnames(design))] <- -1
  test <- functions$glmQLFTest(fit, contrast = contrast)
  table <- as.data.frame(functions$topTags(test, n = Inf, sort.by = "none")$table)
  table$feature <- rownames(table)
  rownames(table) <- NULL
  table
}

.fa_exploratory_contrast <- function(pseudobulk, sample_groups, group_a, group_b) {
  group_a_sum <- Matrix::rowSums(pseudobulk[, sample_groups == group_a, drop = FALSE])
  group_b_sum <- Matrix::rowSums(pseudobulk[, sample_groups == group_b, drop = FALSE])
  cpm_a <- 1e6 * group_a_sum / max(sum(group_a_sum), 1)
  cpm_b <- 1e6 * group_b_sum / max(sum(group_b_sum), 1)
  data.frame(
    feature = rownames(pseudobulk),
    mean_cpm_group_a = as.numeric(cpm_a),
    mean_cpm_group_b = as.numeric(cpm_b),
    logFC = log2((as.numeric(cpm_b) + 0.5) / (as.numeric(cpm_a) + 0.5)),
    stringsAsFactors = FALSE
  )
}

.fa_wilcox_contrast <- function(object, cells, groups, group_a, group_b, cfg, seed) {
  find_markers <- .fa_pkg_fun("Seurat", "FindMarkers")
  normalize <- .fa_pkg_fun("Seurat", "NormalizeData")
  if (is.null(find_markers)) {
    stop("Seurat::FindMarkers is required for the explicitly requested Wilcoxon analysis.")
  }
  cells <- intersect(as.character(cells), colnames(object))
  groups <- as.character(groups[cells])
  selected <- !is.na(groups) & nzchar(groups) & groups %in% c(group_a, group_b)
  cells <- cells[selected]
  groups <- groups[selected]
  if (length(unique(groups)) != 2L || any(table(factor(groups, levels = c(group_a, group_b))) < 3L)) {
    stop("Both Wilcoxon groups need at least three cells in the current stratum.")
  }
  local_object <- tryCatch(object[, cells], error = function(e) e)
  if (inherits(local_object, "error")) stop(conditionMessage(local_object))
  assay <- cfg$assay %||% .fa_default_assay(local_object)
  if (is.null(.fa_matrix(local_object, "data", assay))) {
    if (is.null(normalize)) stop("The selected assay has no normalized data and Seurat::NormalizeData is unavailable.")
    local_object <- normalize(local_object, assay = assay, verbose = FALSE)
  }
  metadata <- data.frame(
    .scRDSreport_wilcox_group = factor(groups, levels = c(group_a, group_b)),
    row.names = cells,
    check.names = FALSE
  )
  local_object <- SeuratObject::AddMetaData(local_object, metadata = metadata)
  max_cells <- cfg$max_cells_per_ident %||% cfg$marker_max_cells_per_ident %||% Inf
  if (is.finite(max_cells)) max_cells <- max(3L, as.integer(max_cells))
  result <- find_markers(
    local_object,
    assay = assay,
    ident.1 = group_b,
    ident.2 = group_a,
    group.by = ".scRDSreport_wilcox_group",
    test.use = "wilcox",
    logfc.threshold = as.numeric(cfg$logfc_threshold %||% 0),
    min.pct = as.numeric(cfg$min_pct %||% 0.1),
    only.pos = FALSE,
    max.cells.per.ident = max_cells,
    random.seed = as.integer(seed),
    verbose = FALSE
  )
  result <- as.data.frame(result, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(result)) stop("Seurat::FindMarkers returned no features for this contrast.")
  result$feature <- rownames(result)
  rownames(result) <- NULL
  logfc_column <- intersect(c("avg_log2FC", "avg_logFC", "log2FC", "logFC"), names(result))
  p_column <- intersect(c("p_val", "PValue", "p.value"), names(result))
  fdr_column <- intersect(c("p_val_adj", "FDR", "padj"), names(result))
  result$logFC <- if (length(logfc_column)) as.numeric(result[[logfc_column[[1L]]]]) else NA_real_
  result$PValue <- if (length(p_column)) as.numeric(result[[p_column[[1L]]]]) else NA_real_
  result$FDR <- if (length(fdr_column)) as.numeric(result[[fdr_column[[1L]]]]) else {
    stats::p.adjust(result$PValue, method = "BH")
  }
  result
}

.fa_module_differential <- function(object, output, cfg, seed, verbose) {
  strategy <- tolower(as.character(cfg$strategy %||% "auto"))
  wilcox_requested <- identical(strategy, "wilcox")
  if (identical(strategy, "none")) {
    return(.fa_result(object, "skipped", "differential_disabled",
                      "Differential analysis was disabled explicitly in report_config()."))
  }
  assay <- cfg$assay %||% .fa_default_assay(object)
  counts <- .fa_matrix(object, "counts", assay)
  if (is.null(counts)) {
    return(.fa_result(object, "skipped", "counts_missing", "Differential analysis requires a counts layer."))
  }
  meta <- .seurat_metadata(object)
  sample_column <- .fa_sample_column(object, cfg)
  group_column <- .fa_group_column(object, cfg)
  if (is.null(group_column) || (!wilcox_requested && is.null(sample_column))) {
    return(.fa_result(
      object, "needs_input", "sample_or_group_missing",
      if (wilcox_requested) {
        "The explicitly requested Wilcoxon analysis needs a per-cell group column. A sample column is recommended but is not required for this exploratory cell-level test."
      } else {
        "Differential analysis needs explicit per-cell sample and group columns. Automatic names are not sufficient when these fields are absent."
      }
    ))
  }
  sample_id <- if (!is.null(sample_column)) {
    as.character(meta[colnames(counts), sample_column])
  } else {
    rep(NA_character_, ncol(counts))
  }
  group_id <- as.character(meta[colnames(counts), group_column])
  valid_group <- !is.na(group_id) & nzchar(group_id)
  valid_sample <- !is.na(sample_id) & nzchar(sample_id)
  valid <- valid_group & (wilcox_requested | valid_sample)
  review_flags <- rep(FALSE, ncol(counts))
  if (identical(group_column, ".scRDSreport_group")) {
    review_column <- ".scRDSreport_design_needs_review"
    review_flags <- if (review_column %in% names(meta)) {
      suppressWarnings(as.logical(meta[colnames(counts), review_column]))
    } else {
      rep(TRUE, ncol(counts))
    }
    review_flags[is.na(review_flags)] <- TRUE
  }
  design_review_required <- !wilcox_requested && any(review_flags[valid])
  if (!wilcox_requested) {
    sample_group_sets <- tapply(group_id[valid], sample_id[valid], function(x) unique(x))
    ambiguous <- names(sample_group_sets)[vapply(sample_group_sets, length, integer(1)) != 1L]
    if (length(ambiguous)) {
      return(.fa_result(
        object, "needs_input", "sample_maps_to_multiple_groups",
        paste0("Each biological sample must map to one group. Ambiguous samples: ", paste(ambiguous, collapse = ", "), ".")
      ))
    }
  }
  strata_column <- cfg$strata_column
  if (isTRUE(cfg$by_annotation) && is.null(strata_column)) strata_column <- .fa_annotation_column(object, cfg)
  strata <- if (!is.null(strata_column) && strata_column %in% names(meta)) {
    as.character(meta[colnames(counts), strata_column])
  } else {
    rep("all_cells", ncol(counts))
  }
  strata[is.na(strata) | !nzchar(strata)] <- "unlabeled"
  stratum_levels <- sort(unique(strata[valid]))
  max_strata <- as.integer(cfg$max_strata %||% 50L)
  if (length(stratum_levels) > max_strata) stratum_levels <- stratum_levels[seq_len(max_strata)]
  artifacts <- list()
  all_results <- list()
  contrast_registry <- list()
  any_inferential <- FALSE
  any_exploratory <- FALSE
  any_wilcox <- FALSE
  any_descriptive <- FALSE
  for (stratum in stratum_levels) {
    cell_keep <- valid & strata == stratum
    if (sum(cell_keep) < as.integer(cfg$min_cells_per_stratum %||% 20L)) next
    matrix <- counts[, cell_keep, drop = FALSE]
    cell_groups <- stats::setNames(group_id[cell_keep], colnames(matrix))
    if (wilcox_requested) {
      pseudobulk <- NULL
      sample_groups <- NULL
      contrasts <- .fa_contrasts(cell_groups, cfg)
    } else {
      samples <- sample_id[cell_keep]
      pseudobulk <- .fa_sparse_group_sum(matrix, samples)
      sample_groups <- vapply(colnames(pseudobulk), function(id) {
        values <- unique(group_id[cell_keep & sample_id == id])
        values[[1L]]
      }, character(1))
      names(sample_groups) <- colnames(pseudobulk)
      contrasts <- .fa_contrasts(sample_groups, cfg)
    }
    if (!nrow(contrasts)) next
    if (!wilcox_requested) {
      bundle_name <- paste0("pseudobulk_counts_", stratum)
      artifacts <- c(artifacts, .fa_write_sparse_bundle(
        pseudobulk, output, "differential", bundle_name,
        paste("Pseudobulk counts", stratum),
        paste0("Sparse feature-by-sample counts aggregated within stratum '", stratum,
               "'. Biological samples, not cells, are the statistical units.")
      ))
      sample_table <- data.frame(
        stratum = stratum, sample = colnames(pseudobulk), group = unname(sample_groups),
        needs_review = vapply(colnames(pseudobulk), function(id) {
          any(review_flags[cell_keep & sample_id == id])
        }, logical(1)),
        library_size = as.numeric(Matrix::colSums(pseudobulk)), stringsAsFactors = FALSE
      )
      artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
        sample_table, output, "differential", paste0("pseudobulk_samples_", stratum),
        paste("Pseudobulk sample design", stratum),
        "One row per biological sample used in this stratum; library_size is the summed feature count.",
        "biological sample"
      )
    }
    for (i in seq_len(nrow(contrasts))) {
      contrast <- contrasts[i, , drop = FALSE]
      if (wilcox_requested) {
        sample_counts <- if (!is.null(sample_column)) {
          vapply(c(contrast$group_a, contrast$group_b), function(group) {
            length(unique(sample_id[cell_keep & group_id == group & valid_sample]))
          }, integer(1))
        } else {
          stats::setNames(c(NA_integer_, NA_integer_), c(contrast$group_a, contrast$group_b))
        }
        names(sample_counts) <- c(contrast$group_a, contrast$group_b)
        result <- tryCatch(
          .fa_wilcox_contrast(
            object, colnames(matrix), cell_groups,
            contrast$group_a, contrast$group_b, cfg, seed + i
          ),
          error = function(e) e
        )
        inferential <- FALSE
        if (inherits(result, "error")) {
          wilcox_error <- conditionMessage(result)
          result <- .fa_exploratory_contrast(matrix, cell_groups, contrast$group_a, contrast$group_b)
          result$analysis_note <- paste0(
            "The explicitly requested Wilcoxon test could not complete; descriptive effect sizes only: ",
            wilcox_error
          )
          analysis_type <- "exploratory_effect_size_only"
          any_descriptive <- TRUE
        } else {
          analysis_type <- "exploratory_cell_level_wilcox"
          any_wilcox <- TRUE
        }
        group_counts <- sample_counts
      } else {
        group_counts <- table(sample_groups[sample_groups %in% c(contrast$group_a, contrast$group_b)])
        has_replicates <- length(group_counts) == 2L && all(group_counts >= 2L)
        inferential <- has_replicates && !design_review_required && .fa_pkg_available("edgeR")
        result <- if (inferential) {
          tryCatch(
            .fa_edgeR_contrast(pseudobulk, sample_groups, contrast$group_a, contrast$group_b, cfg),
            error = function(e) e
          )
        } else {
          .fa_exploratory_contrast(pseudobulk, sample_groups, contrast$group_a, contrast$group_b)
        }
        if (inherits(result, "error")) {
          edge_error <- conditionMessage(result)
          inferential <- FALSE
          result <- .fa_exploratory_contrast(pseudobulk, sample_groups, contrast$group_a, contrast$group_b)
          result$analysis_note <- paste0("edgeR could not complete; descriptive effect sizes only: ", edge_error)
        }
        analysis_type <- if (inferential) "pseudobulk_edgeR_QL" else "exploratory_effect_size_only"
        if (!inferential && design_review_required) {
          result$analysis_note <- paste0(
            "Automatic sample grouping is marked needs_review; only descriptive effect sizes were exported. ",
            "Provide a verified sample_map before formal replicate-aware inference."
          )
        }
        if (!inferential) any_descriptive <- TRUE
      }
      result$stratum <- stratum
      result$contrast <- contrast$name
      result$group_a <- contrast$group_a
      result$group_b <- contrast$group_b
      result$statistical_unit <- if (identical(analysis_type, "exploratory_cell_level_wilcox")) {
        "cell"
      } else if (identical(analysis_type, "pseudobulk_edgeR_QL")) {
        "biological_sample"
      } else {
        "none_descriptive"
      }
      result$analysis_type <- analysis_type
      result$n_samples_group_a <- as.integer(group_counts[contrast$group_a] %||% 0L)
      result$n_samples_group_b <- as.integer(group_counts[contrast$group_b] %||% 0L)
      result$n_cells_group_a <- as.integer(sum(cell_groups == contrast$group_a, na.rm = TRUE))
      result$n_cells_group_b <- as.integer(sum(cell_groups == contrast$group_b, na.rm = TRUE))
      if (inferential) any_inferential <- TRUE else any_exploratory <- TRUE
      key <- paste(.safe_name(stratum), .safe_name(contrast$name), sep = "__")
      all_results[[key]] <- result
      contrast_registry[[length(contrast_registry) + 1L]] <- data.frame(
        stratum = stratum, contrast = contrast$name, group_a = contrast$group_a,
        group_b = contrast$group_b,
        n_samples_group_a = as.integer(group_counts[contrast$group_a] %||% 0L),
        n_samples_group_b = as.integer(group_counts[contrast$group_b] %||% 0L),
        n_cells_group_a = as.integer(sum(cell_groups == contrast$group_a, na.rm = TRUE)),
        n_cells_group_b = as.integer(sum(cell_groups == contrast$group_b, na.rm = TRUE)),
        analysis_type = analysis_type,
        stringsAsFactors = FALSE
      )
      artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
        result, output, "differential", paste0("de_", key),
        paste("Differential result", stratum, contrast$name),
        if (inferential) {
          "One row per tested feature. edgeR quasi-likelihood statistics use biological pseudobulk samples as replicates."
        } else if (identical(analysis_type, "exploratory_cell_level_wilcox")) {
          "One row per tested feature. Seurat's Wilcoxon test was explicitly requested and treats cells as statistical units; its P values are exploratory because cells do not replace biological replicates."
        } else {
          "One row per feature. No valid biological replication was available, so only descriptive CPM and log2 fold change are reported; no P value was manufactured."
        },
        "feature-contrast record",
        list(
          analysis_type = "pseudobulk_edgeR_QL is replicate-aware inference; exploratory_cell_level_wilcox is an explicitly requested cell-level test; exploratory_effect_size_only has no formal test.",
          statistical_unit = "Unit used by the method: biological_sample for pseudobulk, cell for explicitly requested Wilcoxon, or none_descriptive when no significance test was run."
        )
      )
    }
  }
  if (!length(all_results)) {
    return(.fa_result(
      object, "needs_input", "no_valid_contrast",
      "No stratum contained enough cells and at least two valid groups for a contrast.",
      artifacts = artifacts,
      details = list(sample_column = sample_column, group_column = group_column, strata_column = strata_column)
    ))
  }
  registry <- do.call(rbind, contrast_registry)
  artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
    registry, output, "differential", "contrast_registry", "Differential contrast registry",
    "One row per executed stratum and contrast. The analysis_type field distinguishes replicate-aware inference, explicitly requested exploratory cell-level Wilcoxon tests, and descriptive output.",
    "contrast", list(
      n_samples_group_a = "Number of biological samples in reference group A.",
      n_samples_group_b = "Number of biological samples in comparison group B.",
      n_cells_group_a = "Number of cells in reference group A within this stratum.",
      n_cells_group_b = "Number of cells in comparison group B within this stratum.",
      analysis_type = "Method actually used. Cell-level Wilcoxon P values are exploratory; effect-size-only results contain no P values."
    )
  )
  rankings <- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
      feature = x$feature,
      logFC = x$logFC,
      PValue = if ("PValue" %in% names(x)) x$PValue else NA_real_,
      FDR = if ("FDR" %in% names(x)) x$FDR else NA_real_,
      stratum = x$stratum,
      contrast = x$contrast,
      analysis_type = x$analysis_type,
      stringsAsFactors = FALSE
    )
  }))
  for (column in c("PValue", "FDR")) {
    if (column %in% names(rankings) && !any(is.finite(rankings[[column]]))) {
      rankings[[column]] <- NULL
    }
  }
  object <- .fa_store_analysis_misc(object, "differential_rankings", rankings)
  status <- if (design_review_required && !wilcox_requested && !any_inferential) {
    "needs_input"
  } else if ((any_inferential && any_exploratory) ||
                (any_wilcox && any_descriptive) ||
                (wilcox_requested && !any_wilcox)) "partial" else "completed"
  reason <- if (wilcox_requested && !any_wilcox) {
    "wilcox_failed_descriptive_only"
  } else if (any_wilcox && any_descriptive) {
    "mixed_wilcox_and_descriptive"
  } else if (any_wilcox) {
    "explicit_exploratory_wilcox_completed"
  } else if (any_inferential && any_exploratory) {
    "mixed_inferential_and_exploratory"
  } else if (any_inferential) {
    "replicate_aware_pseudobulk_completed"
  } else if (design_review_required) {
    "sample_design_needs_review"
  } else {
    "exploratory_no_valid_replication"
  }
  message <- if (wilcox_requested && !any_wilcox) {
    "The explicitly requested Wilcoxon tests could not complete; only descriptive group effect sizes were exported."
  } else if (any_wilcox) {
    "The explicitly requested Seurat Wilcoxon tests were run at cell level. Their P values are exploratory and do not substitute for biological replication."
  } else if (any_inferential) {
    "Replicate-aware edgeR pseudobulk results were generated where both groups had at least two biological samples."
  } else if (design_review_required) {
    paste0(
      "Automatic sample grouping requires review; descriptive effect sizes were exported without formal P values. ",
      "Provide a verified sample_map to enable replicate-aware inference."
    )
  } else {
    "Biological replication was insufficient; only descriptive effect sizes were exported and no formal P values were created."
  }
  .fa_result(
    object, status, reason, message,
    if (any_wilcox) {
      "Seurat::FindMarkers(test.use = 'wilcox')"
    } else if (any_inferential) {
      "edgeR::glmQLFTest"
    } else if (wilcox_requested) {
      "descriptive group aggregate"
    } else {
      "descriptive pseudobulk"
    },
    artifacts,
    list(sample_column = sample_column, group_column = group_column,
         strata_column = strata_column %||% "all_cells", inferential = any_inferential,
         exploratory = any_exploratory, wilcox = any_wilcox,
         design_needs_review = design_review_required,
         requested_strategy = strategy)
  )
}

# Functional enrichment ------------------------------------------------------

.fa_read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  sets <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = FALSE)[[1L]]
    if (length(fields) < 3L) return(NULL)
    genes <- unique(fields[-c(1L, 2L)])
    genes <- genes[nzchar(genes)]
    if (!length(genes)) return(NULL)
    stats::setNames(list(genes), fields[[1L]])
  })
  sets <- Filter(Negate(is.null), sets)
  if (!length(sets)) return(list())
  unlist(sets, recursive = FALSE)
}

.fa_gene_sets <- function(cfg, species = NULL, resources = list()) {
  sets <- cfg[["gene_sets", exact = TRUE]] %||% list()
  if (!is.list(sets)) sets <- list()
  paths <- unique(unlist(c(
    cfg[["gmt", exact = TRUE]], cfg[["gmt_files", exact = TRUE]]
  ), use.names = FALSE))
  if (length(paths)) {
    paths <- as.character(paths)
    paths <- paths[!is.na(paths) & file.exists(paths)]
  } else {
    paths <- character()
  }
  for (path in paths) sets <- c(sets, .fa_read_gmt(path))
  has_user_sets <- length(sets) > 0L
  source <- if (has_user_sets) "user_gene_sets_or_gmt" else "none"
  auto_error <- NULL
  auto_membership <- NULL
  database_versions <- character()
  db_species <- toupper(as.character(
    cfg$msigdbr_db_species %||% resources$msigdbr_db_species %||%
      if (identical(species, "mouse")) "MM" else "HS"
  )[[1L]])
  default_collection <- resources$msigdbr_default_collection %||%
    if (identical(db_species, "MM")) "MH" else "H"
  collections <- if (has_user_sets) {
    unique(as.character(cfg$gene_set_collections %||% "user_supplied"))
  } else {
    unique(as.character(cfg$msigdb_collections %||% default_collection))
  }
  collections <- collections[!is.na(collections) & nzchar(collections)]
  if (!length(sets) && !isFALSE(cfg$auto_gene_sets) && length(collections) &&
      .fa_pkg_available("msigdbr")) {
    loader <- .fa_pkg_fun("msigdbr", "msigdbr")
    target_species <- cfg$msigdbr_species %||% resources$msigdbr_species %||%
      resources$scientific_name %||% species
    loaded <- lapply(collections, function(collection) {
      tryCatch(
        suppressMessages(loader(
          db_species = db_species,
          species = as.character(target_species)[[1L]],
          collection = collection
        )),
        error = function(e) {
          auto_error <<- paste0(collection, ": ", conditionMessage(e))
          NULL
        }
      )
    })
    loaded <- Filter(function(x) !is.null(x) && nrow(x) &&
                       all(c("gs_name", "gene_symbol") %in% names(x)), loaded)
    if (length(loaded)) {
      table <- do.call(rbind, loaded)
      table <- table[!is.na(table$gs_name) & nzchar(table$gs_name) &
                       !is.na(table$gene_symbol) & nzchar(table$gene_symbol), , drop = FALSE]
      sets <- split(as.character(table$gene_symbol), as.character(table$gs_name))
      sets <- lapply(sets, unique)
      membership_columns <- intersect(
        c(
          "gs_name", "gene_symbol", "gs_collection", "gs_subcollection",
          "db_version", "db_gene_symbol", "db_ensembl_gene",
          "ortholog_sources", "num_ortholog_sources"
        ),
        names(table)
      )
      auto_membership <- unique(as.data.frame(table[membership_columns], stringsAsFactors = FALSE))
      if ("db_version" %in% names(table)) {
        database_versions <- unique(as.character(table$db_version[!is.na(table$db_version)]))
      }
      source <- "msigdbr"
    }
    attr(sets, "msigdbr_species") <- target_species
    attr(sets, "msigdbr_db_species") <- db_species
  }
  sets <- sets[vapply(sets, function(x) length(unique(x[!is.na(x) & nzchar(x)])) >= 2L, logical(1))]
  sets <- lapply(sets, unique)
  attr(sets, "source") <- source
  attr(sets, "collections") <- collections
  attr(sets, "auto_error") <- auto_error
  attr(sets, "membership") <- auto_membership
  attr(sets, "database_versions") <- database_versions
  attr(sets, "gmt_files") <- paths
  if (is.null(attr(sets, "msigdbr_species"))) {
    attr(sets, "msigdbr_species") <- cfg$gene_set_species %||%
      resources$scientific_name %||% resources$msigdbr_species
  }
  if (is.null(attr(sets, "msigdbr_db_species"))) {
    attr(sets, "msigdbr_db_species") <- if (has_user_sets) {
      cfg$gene_set_database_species %||% NA_character_
    } else {
      cfg$msigdbr_db_species %||% resources$msigdbr_db_species
    }
  }
  sets
}

.fa_orgdb <- function(resources, cfg) {
  if (!is.null(cfg$orgdb) && !is.character(cfg$orgdb)) return(cfg$orgdb)
  package <- cfg$orgdb %||% resources$orgdb
  if (is.null(package) || !.fa_pkg_available(package)) return(NULL)
  object_name <- sub("[.]db$", ".db", package)
  object <- .fa_pkg_object(package, object_name)
  if (is.null(object)) object <- .fa_pkg_object(package, package)
  object
}

.fa_feature_mapping <- function(features, feature_symbols, orgdb,
                                preferred_keytype = NULL,
                                symbol_column = "SYMBOL") {
  base <- data.frame(
    feature = as.character(features),
    symbol_from_object = unname(feature_symbols[features]),
    stringsAsFactors = FALSE
  )
  base$symbol_from_object[is.na(base$symbol_from_object) | !nzchar(base$symbol_from_object)] <-
    base$feature[is.na(base$symbol_from_object) | !nzchar(base$symbol_from_object)]
  if (is.null(orgdb) || !.fa_pkg_available("AnnotationDbi")) {
    base$SYMBOL <- base$symbol_from_object
    base$ENTREZID <- NA_character_
    base$mapping_keytype <- NA_character_
    return(base)
  }
  keytypes_fun <- .fa_pkg_fun("AnnotationDbi", "keytypes")
  columns_fun <- .fa_pkg_fun("AnnotationDbi", "columns")
  select_fun <- .fa_pkg_fun("AnnotationDbi", "select")
  available <- tryCatch(keytypes_fun(orgdb), error = function(e) character())
  available_columns <- tryCatch(columns_fun(orgdb), error = function(e) available)
  cleaned <- sub("[.][0-9]+$", "", base$feature)
  ensembl_fraction <- mean(grepl("^ENS[A-Z]*G[0-9]+", cleaned, ignore.case = TRUE))
  entrez_fraction <- mean(grepl("^[0-9]+$", cleaned))
  heuristic <- if (ensembl_fraction >= 0.2) {
    "ENSEMBL"
  } else if (entrez_fraction >= 0.8) {
    "ENTREZID"
  } else {
    "SYMBOL"
  }
  preferred_keytype <- toupper(as.character(preferred_keytype %||% character()))
  candidates <- unique(c(heuristic, preferred_keytype, "SYMBOL", "ENSEMBL", "ENTREZID"))
  candidates <- candidates[candidates %in% available]
  requested_symbol_column <- toupper(as.character(symbol_column %||% "SYMBOL")[[1L]])
  query_columns <- unique(c(requested_symbol_column, "SYMBOL", "ENTREZID"))
  query_columns <- query_columns[query_columns %in% available_columns]
  if (!length(candidates) || !length(query_columns)) {
    base$SYMBOL <- base$symbol_from_object
    base$ENTREZID <- NA_character_
    base$mapping_keytype <- NA_character_
    return(base)
  }
  best <- NULL
  best_score <- -1L
  for (keytype in candidates) {
    keys <- if (identical(keytype, "SYMBOL")) base$symbol_from_object else cleaned
    mapped <- tryCatch(
      suppressMessages(select_fun(
        orgdb, keys = unique(keys), keytype = keytype, columns = query_columns
      )),
      error = function(e) NULL
    )
    if (is.null(mapped) || !nrow(mapped) || !keytype %in% names(mapped)) next
    mapped <- as.data.frame(mapped, stringsAsFactors = FALSE)
    mapped <- mapped[!duplicated(mapped[[keytype]]), , drop = FALSE]
    index <- match(keys, mapped[[keytype]])
    symbol_name <- intersect(c(requested_symbol_column, "SYMBOL"), names(mapped))
    mapped_symbols <- if (length(symbol_name)) {
      as.character(mapped[[symbol_name[[1L]]]][index])
    } else {
      rep(NA_character_, length(index))
    }
    score <- sum(!is.na(mapped_symbols) & nzchar(mapped_symbols))
    if (score > best_score) {
      best <- list(
        keytype = keytype,
        symbols = mapped_symbols,
        entrez = if ("ENTREZID" %in% names(mapped)) {
          as.character(mapped$ENTREZID[index])
        } else if (identical(keytype, "ENTREZID")) {
          as.character(keys)
        } else {
          rep(NA_character_, length(index))
        }
      )
      best_score <- score
    }
    if (score >= ceiling(0.8 * nrow(base))) break
  }
  if (is.null(best) || best_score < 1L) {
    base$SYMBOL <- base$symbol_from_object
    base$ENTREZID <- NA_character_
    base$mapping_keytype <- NA_character_
    return(base)
  }
  base$SYMBOL <- best$symbols
  base$SYMBOL[is.na(base$SYMBOL) | !nzchar(base$SYMBOL)] <- base$symbol_from_object[is.na(base$SYMBOL) | !nzchar(base$SYMBOL)]
  base$ENTREZID <- best$entrez
  base$mapping_keytype <- best$keytype
  base
}

.fa_enrichment_plot <- function(table, title) {
  if (!nrow(table) || !"Description" %in% names(table)) return(NULL)
  p_column <- intersect(c("p.adjust", "pvalue", "qvalue"), names(table))
  size_column <- intersect(c("Count", "setSize", "core_enrichment"), names(table))
  if (!length(p_column)) return(NULL)
  shown <- table[order(table[[p_column[[1L]]]], na.last = NA), , drop = FALSE]
  shown <- utils::head(shown, 20L)
  if (!nrow(shown)) return(NULL)
  shown$minus_log10 <- -log10(pmax(as.numeric(shown[[p_column[[1L]]]]), .Machine$double.xmin))
  shown$plot_size <- if (length(size_column) && is.numeric(shown[[size_column[[1L]]]])) {
    as.numeric(shown[[size_column[[1L]]]])
  } else {
    seq_len(nrow(shown))
  }
  ggplot2::ggplot(shown, ggplot2::aes(x = minus_log10, y = stats::reorder(Description, minus_log10),
                                      size = plot_size, color = minus_log10)) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_color_gradient(low = "#56B4E9", high = "#D55E00") +
    ggplot2::labs(x = paste0("-log10(", p_column[[1L]], ")"), y = NULL,
                  size = "Gene/set size", color = "Evidence", title = title) +
    ggplot2::theme_bw(base_size = 11)
}

.fa_enrichment_artifacts <- function(result, output, name, label, description) {
  if (is.null(result)) return(list())
  table <- tryCatch(as.data.frame(result), error = function(e) data.frame())
  if (!nrow(table)) return(list())
  artifacts <- list(.fa_write_table_artifact(
    table, output, "enrichment", name, label, description,
    "gene-set result", list(
      Description = "Gene-set or pathway description.",
      p.adjust = "Multiple-testing adjusted enrichment P value.",
      geneID = "Input genes contributing to this result; delimiter follows the enrichment engine.",
      core_enrichment = "Leading-edge genes for a GSEA result."
    )
  ))
  plot <- .fa_enrichment_plot(table, label)
  if (!is.null(plot)) artifacts <- c(artifacts, .fa_write_plot_artifacts(
    plot, output, "enrichment", name, label,
    paste0(description, " The plot shows at most the 20 most statistically supported rows.")
  ))
  artifacts
}

.fa_run_gsva <- function(object, cfg, gene_sets, resources = list()) {
  if (!length(gene_sets) || !.fa_pkg_available("GSVA")) return(NULL)
  counts <- .fa_matrix(object, "counts", cfg$assay %||% .fa_default_assay(object))
  meta <- .seurat_metadata(object)
  sample_column <- .fa_sample_column(object, cfg)
  if (is.null(counts) || is.null(sample_column)) return(NULL)
  samples <- as.character(meta[colnames(counts), sample_column])
  valid <- !is.na(samples) & nzchar(samples)
  if (length(unique(samples[valid])) < 2L) return(NULL)
  pseudobulk <- .fa_sparse_group_sum(counts[, valid, drop = FALSE], samples[valid])
  symbols <- .fa_feature_symbols(object, cfg$assay %||% .fa_default_assay(object))
  orgdb <- .fa_orgdb(resources, cfg)
  mapping <- .fa_feature_mapping(
    rownames(pseudobulk), symbols, orgdb,
    preferred_keytype = cfg$feature_keytype %||% resources$feature_keytype,
    symbol_column = cfg$symbol_column %||% resources$symbol_column %||% "SYMBOL"
  )
  mapped_symbols <- as.character(mapping$SYMBOL[match(rownames(pseudobulk), mapping$feature)])
  keep <- !is.na(mapped_symbols) & nzchar(mapped_symbols) & !duplicated(toupper(mapped_symbols))
  pseudobulk <- pseudobulk[keep, , drop = FALSE]
  rownames(pseudobulk) <- mapped_symbols[keep]
  library_sizes <- Matrix::colSums(pseudobulk)
  expression <- as.matrix(log1p(.fa_scale_columns(pseudobulk, 1e6 / pmax(library_sizes, 1))))
  gsva_fun <- .fa_pkg_fun("GSVA", "gsva")
  param_fun <- .fa_pkg_fun("GSVA", "gsvaParam")
  scores <- tryCatch({
    if (!is.null(param_fun)) {
      parameter <- param_fun(exprData = expression, geneSets = gene_sets, kcdf = "Gaussian")
      gsva_fun(parameter, verbose = FALSE)
    } else {
      gsva_fun(expr = expression, gset.idx.list = gene_sets, method = "gsva", kcdf = "Gaussian", verbose = FALSE)
    }
  }, error = function(e) NULL)
  if (is.null(scores)) return(NULL)
  list(scores = scores, samples = colnames(scores))
}

.fa_module_enrichment <- function(object, output, cfg, seed, verbose, species) {
  rankings <- .fa_get_analysis_misc(object, "differential_rankings")
  resources <- .fa_species_resources(species, cfg)
  gene_sets <- .fa_gene_sets(cfg, species = species, resources = resources)
  gene_set_source <- attr(gene_sets, "source") %||% "none"
  gene_set_collections <- attr(gene_sets, "collections") %||% character()
  gene_set_error <- attr(gene_sets, "auto_error") %||% NULL
  gene_set_membership <- attr(gene_sets, "membership") %||% NULL
  gene_set_db_versions <- attr(gene_sets, "database_versions") %||% character()
  gsva_result <- .fa_run_gsva(object, cfg, gene_sets, resources = resources)
  orgdb <- .fa_orgdb(resources, cfg)
  cluster_profiler_available <- .fa_pkg_available("clusterProfiler")
  artifacts <- list()
  gene_set_provenance <- data.frame(
    source = gene_set_source,
    species = as.character(attr(gene_sets, "msigdbr_species") %||% resources$scientific_name %||% species),
    database_species = as.character(attr(gene_sets, "msigdbr_db_species") %||% NA_character_),
    collections = paste(gene_set_collections, collapse = ","),
    resource_version = if (length(gene_set_db_versions)) paste(gene_set_db_versions, collapse = ",") else NA_character_,
    package_version = if (.fa_pkg_available("msigdbr")) .fa_package_version("msigdbr") else NA_character_,
    ortholog_projection = if (identical(gene_set_source, "msigdbr")) {
      isTRUE(resources$msigdbr_ortholog_projection)
    } else {
      NA
    },
    n_gene_sets = length(gene_sets),
    n_unique_genes = length(unique(unlist(gene_sets, use.names = FALSE))),
    loading_note = as.character(gene_set_error %||% "ok"),
    stringsAsFactors = FALSE
  )
  artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
    gene_set_provenance, output, "enrichment", "gene_set_provenance",
    "Gene-set resource provenance",
    "One row describing the gene-set source used by ORA, GSEA, or GSVA. For non-human species, the database_species and target species distinguish native collections from ortholog-mapped collections.",
    "gene-set resource",
    list(n_gene_sets = "Number of loaded gene sets.",
         n_unique_genes = "Number of unique target-species gene symbols across the loaded sets.",
         loading_note = "ok, or the captured resource-loading error when no automatic set could be loaded.")
  )
  if (length(gene_sets)) {
    if (is.null(gene_set_membership) || !nrow(gene_set_membership)) {
      gene_set_membership <- do.call(rbind, lapply(names(gene_sets), function(term) {
        data.frame(gs_name = term, gene_symbol = as.character(gene_sets[[term]]), stringsAsFactors = FALSE)
      }))
    }
    gene_set_membership$resource_source <- gene_set_source
    gene_set_membership$target_species <- as.character(
      attr(gene_sets, "msigdbr_species") %||% resources$scientific_name %||% species
    )
    gene_set_membership$database_species <- as.character(
      attr(gene_sets, "msigdbr_db_species") %||% NA_character_
    )
    artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
      gene_set_membership, output, "enrichment", "gene_set_membership",
      "Gene-set membership",
      paste0(
        "One row per gene-set membership used by enrichment and GSVA. ",
        if (identical(gene_set_source, "msigdbr") &&
            isTRUE(resources$msigdbr_ortholog_projection)) {
          "The source database is human MSigDB and target symbols are ortholog-projected; ortholog evidence columns are retained when provided by msigdbr."
        } else if (identical(gene_set_source, "msigdbr")) {
          "The selected MSigDB database matches the registered database species."
        } else {
          "Rows came from explicit user gene sets or GMT files; no MSigDB database origin or ortholog projection is inferred."
        }
      ),
      "gene-set membership",
      list(gs_name = "Gene-set identifier.", gene_symbol = "Target-species symbol used for matching.",
           resource_source = "msigdbr or explicit user gene-set/GMT source.",
           target_species = "Species whose symbols appear in gene_symbol.",
           database_species = "MSigDB database species code when applicable.")
    )
  }
  completed_types <- character()
  skipped_reasons <- character()
  if (!is.null(rankings) && nrow(rankings) && cluster_profiler_available) {
    symbols <- .fa_feature_symbols(object, cfg$assay %||% .fa_default_assay(object))
    mapping <- .fa_feature_mapping(
      unique(rankings$feature), symbols, orgdb,
      preferred_keytype = cfg$feature_keytype %||% resources$feature_keytype,
      symbol_column = cfg$symbol_column %||% resources$symbol_column %||% "SYMBOL"
    )
    artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
      mapping, output, "enrichment", "feature_identifier_mapping", "Feature identifier mapping",
      "One row per differential feature. SYMBOL/ENTREZID are resolved with the selected species OrgDb when available; no other species is substituted.",
      "feature", list(symbol_from_object = "Symbol retained in assay feature metadata, or the feature ID when unavailable.",
                      ENTREZID = "Species-matched Entrez identifier; NA when no validated mapping exists.")
    )
    rankings <- merge(rankings, mapping, by = "feature", all.x = TRUE, sort = FALSE)
    keys <- unique(rankings[c("stratum", "contrast")])
    max_jobs <- as.integer(cfg$max_enrichment_contrasts %||% 8L)
    if (nrow(keys) > max_jobs) keys <- keys[seq_len(max_jobs), , drop = FALSE]
    enrich_go <- .fa_pkg_fun("clusterProfiler", "enrichGO")
    enrich_kegg <- .fa_pkg_fun("clusterProfiler", "enrichKEGG")
    gse_go <- .fa_pkg_fun("clusterProfiler", "gseGO")
    gse_kegg <- .fa_pkg_fun("clusterProfiler", "gseKEGG")
    enricher <- .fa_pkg_fun("clusterProfiler", "enricher")
    gsea <- .fa_pkg_fun("clusterProfiler", "GSEA")
    term2gene <- if (length(gene_sets)) {
      do.call(rbind, lapply(names(gene_sets), function(name) {
        data.frame(term = name, gene = as.character(gene_sets[[name]]), stringsAsFactors = FALSE)
      }))
    } else NULL
    for (i in seq_len(nrow(keys))) {
      subset <- rankings[rankings$stratum == keys$stratum[i] & rankings$contrast == keys$contrast[i], , drop = FALSE]
      if (!any(subset$analysis_type == "pseudobulk_edgeR_QL")) {
        skipped_reasons <- c(skipped_reasons, paste(keys$stratum[i], keys$contrast[i], "has exploratory DE only"))
        next
      }
      suffix <- paste(.safe_name(keys$stratum[i]), .safe_name(keys$contrast[i]), sep = "__")
      significant <- subset[is.finite(subset$FDR) & subset$FDR <= as.numeric(cfg$fdr %||% 0.05) &
                              abs(subset$logFC) >= as.numeric(cfg$min_abs_logfc %||% 0), , drop = FALSE]
      entrez <- unique(significant$ENTREZID[!is.na(significant$ENTREZID) & nzchar(significant$ENTREZID)])
      symbols_sig <- unique(significant$SYMBOL[!is.na(significant$SYMBOL) & nzchar(significant$SYMBOL)])
      if (length(entrez) >= as.integer(cfg$min_genes %||% 5L) && !is.null(orgdb) && !is.null(enrich_go)) {
        result <- tryCatch(enrich_go(
          gene = entrez, OrgDb = orgdb, keyType = "ENTREZID", ont = cfg$ontology %||% "ALL",
          pAdjustMethod = "BH", pvalueCutoff = cfg$pvalue_cutoff %||% 0.05,
          qvalueCutoff = cfg$qvalue_cutoff %||% 0.2, readable = TRUE
        ), error = function(e) NULL)
        artifacts <- c(artifacts, .fa_enrichment_artifacts(
          result, output, paste0("go_ora_", suffix), paste("GO ORA", keys$stratum[i], keys$contrast[i]),
          "Over-representation analysis of replicate-aware DE genes using the species-matched OrgDb."
        ))
        if (!is.null(result) && nrow(as.data.frame(result))) completed_types <- c(completed_types, "GO ORA")
      }
      if (length(entrez) >= as.integer(cfg$min_genes %||% 5L) && !is.null(resources$kegg) && !is.null(enrich_kegg)) {
        result <- tryCatch(enrich_kegg(
          gene = entrez, organism = resources$kegg, pvalueCutoff = cfg$pvalue_cutoff %||% 0.05,
          pAdjustMethod = "BH", qvalueCutoff = cfg$qvalue_cutoff %||% 0.2
        ), error = function(e) NULL)
        artifacts <- c(artifacts, .fa_enrichment_artifacts(
          result, output, paste0("kegg_ora_", suffix), paste("KEGG ORA", keys$stratum[i], keys$contrast[i]),
          "KEGG over-representation analysis using the selected species code; online KEGG failures are recorded as missing output, never replaced with another species."
        ))
        if (!is.null(result) && nrow(as.data.frame(result))) completed_types <- c(completed_types, "KEGG ORA")
      }
      if (length(symbols_sig) >= as.integer(cfg$min_genes %||% 5L) && !is.null(term2gene) && !is.null(enricher)) {
        result <- tryCatch(enricher(
          gene = symbols_sig, TERM2GENE = term2gene,
          pvalueCutoff = cfg$pvalue_cutoff %||% 0.05, pAdjustMethod = "BH"
        ), error = function(e) NULL)
        artifacts <- c(artifacts, .fa_enrichment_artifacts(
          result, output, paste0("custom_ora_", suffix), paste("Custom gene-set ORA", keys$stratum[i], keys$contrast[i]),
          "Over-representation analysis against explicitly supplied GMT/gene sets."
        ))
        if (!is.null(result) && nrow(as.data.frame(result))) completed_types <- c(completed_types, "custom ORA")
      }
      ranking_entrez <- subset[!is.na(subset$ENTREZID) & nzchar(subset$ENTREZID) & is.finite(subset$logFC), , drop = FALSE]
      ranking_entrez <- ranking_entrez[order(abs(ranking_entrez$logFC), decreasing = TRUE), , drop = FALSE]
      ranking_entrez <- ranking_entrez[!duplicated(ranking_entrez$ENTREZID), , drop = FALSE]
      gene_list_entrez <- stats::setNames(ranking_entrez$logFC, ranking_entrez$ENTREZID)
      gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
      if (length(gene_list_entrez) >= 50L && !is.null(orgdb) && !is.null(gse_go)) {
        result <- tryCatch(gse_go(
          geneList = gene_list_entrez, OrgDb = orgdb, keyType = "ENTREZID",
          ont = cfg$ontology %||% "ALL", pAdjustMethod = "BH",
          minGSSize = cfg$min_gs_size %||% 10L, maxGSSize = cfg$max_gs_size %||% 500L,
          seed = TRUE, verbose = FALSE
        ), error = function(e) NULL)
        artifacts <- c(artifacts, .fa_enrichment_artifacts(
          result, output, paste0("go_gsea_", suffix), paste("GO GSEA", keys$stratum[i], keys$contrast[i]),
          "Rank-based GO enrichment from replicate-aware pseudobulk log fold changes."
        ))
        if (!is.null(result) && nrow(as.data.frame(result))) completed_types <- c(completed_types, "GO GSEA")
      }
      if (length(gene_list_entrez) >= 50L && !is.null(resources$kegg) && !is.null(gse_kegg)) {
        result <- tryCatch(gse_kegg(
          geneList = gene_list_entrez, organism = resources$kegg,
          minGSSize = cfg$min_gs_size %||% 10L, maxGSSize = cfg$max_gs_size %||% 500L,
          pAdjustMethod = "BH", seed = TRUE, verbose = FALSE
        ), error = function(e) NULL)
        artifacts <- c(artifacts, .fa_enrichment_artifacts(
          result, output, paste0("kegg_gsea_", suffix), paste("KEGG GSEA", keys$stratum[i], keys$contrast[i]),
          "Rank-based KEGG enrichment using the selected species code."
        ))
        if (!is.null(result) && nrow(as.data.frame(result))) completed_types <- c(completed_types, "KEGG GSEA")
      }
      ranking_symbol <- subset[!is.na(subset$SYMBOL) & nzchar(subset$SYMBOL) & is.finite(subset$logFC), , drop = FALSE]
      ranking_symbol <- ranking_symbol[order(abs(ranking_symbol$logFC), decreasing = TRUE), , drop = FALSE]
      ranking_symbol <- ranking_symbol[!duplicated(ranking_symbol$SYMBOL), , drop = FALSE]
      gene_list_symbol <- sort(stats::setNames(ranking_symbol$logFC, ranking_symbol$SYMBOL), decreasing = TRUE)
      if (length(gene_list_symbol) >= 50L && !is.null(term2gene) && !is.null(gsea)) {
        result <- tryCatch(gsea(
          geneList = gene_list_symbol, TERM2GENE = term2gene,
          minGSSize = cfg$min_gs_size %||% 10L, maxGSSize = cfg$max_gs_size %||% 500L,
          pAdjustMethod = "BH", seed = TRUE, verbose = FALSE
        ), error = function(e) NULL)
        artifacts <- c(artifacts, .fa_enrichment_artifacts(
          result, output, paste0("custom_gsea_", suffix), paste("Custom GSEA", keys$stratum[i], keys$contrast[i]),
          "Rank-based enrichment against explicitly supplied GMT/gene sets."
        ))
        if (!is.null(result) && nrow(as.data.frame(result))) completed_types <- c(completed_types, "custom GSEA")
      }
    }
  } else {
    if (is.null(rankings) || !nrow(rankings)) skipped_reasons <- c(skipped_reasons, "no differential ranking")
    if (!cluster_profiler_available) skipped_reasons <- c(skipped_reasons, "clusterProfiler is not installed")
  }
  if (!is.null(gsva_result)) {
    scores <- as.data.frame(gsva_result$scores, check.names = FALSE)
    scores <- data.frame(gene_set = rownames(scores), scores, row.names = NULL, check.names = FALSE)
    artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
      scores, output, "enrichment", "gsva_scores", "Pseudobulk GSVA scores",
      "One row per supplied gene set and one column per biological sample. Scores use log-normalized pseudobulk expression; cells are not treated as replicates.",
      "gene set", list(gene_set = "Name from the supplied GMT or gene_sets configuration.")
    )
    completed_types <- c(completed_types, "GSVA")
  } else if (!length(gene_sets)) {
    skipped_reasons <- c(
      skipped_reasons,
      if (is.null(gene_set_error)) {
        "GSVA needs a species-matched msigdbr resource or explicitly supplied GMT/gene_sets"
      } else {
        paste0("automatic gene-set loading failed: ", gene_set_error)
      }
    )
  }
  completed_types <- unique(completed_types)
  if (!length(completed_types)) {
    status <- if (length(skipped_reasons)) "needs_input" else "skipped"
    return(.fa_result(
      object, status, "no_enrichment_result",
      paste("No enrichment result was generated:", paste(unique(skipped_reasons), collapse = "; ")),
      if (cluster_profiler_available) "clusterProfiler/GSVA" else NULL,
      artifacts,
      list(species = species, resources = .fa_sanitize_config(resources), skipped = unique(skipped_reasons))
    ))
  }
  .fa_result(
    object, if (length(skipped_reasons)) "partial" else "completed", "enrichment_completed",
    paste("Generated:", paste(completed_types, collapse = ", "), "."),
    "clusterProfiler/GSVA", artifacts,
    list(species = species, completed_types = completed_types, skipped = unique(skipped_reasons),
         gene_set_source = gene_set_source, gene_sets = length(gene_sets))
  )
}

# Monocle3 trajectory geometry and pseudotime --------------------------------

.fa_build_cds <- function(object, cfg, seed) {
  max_cells <- as.integer(cfg$max_cells %||% 20000L)
  cells <- colnames(object)
  if (length(cells) > max_cells) {
    meta <- .seurat_metadata(object)
    stratify_column <- .fa_cluster_column(object, cfg) %||% .fa_annotation_column(object, cfg)
    set.seed(as.integer(seed))
    if (!is.null(stratify_column)) {
      strata <- split(cells, as.character(meta[cells, stratify_column]))
      target <- pmax(1L, floor(max_cells * lengths(strata) / length(cells)))
      selected <- unlist(Map(function(x, n) sample(x, min(length(x), n)), strata, target), use.names = FALSE)
      if (length(selected) < max_cells) {
        remaining <- setdiff(cells, selected)
        selected <- c(selected, sample(remaining, min(length(remaining), max_cells - length(selected))))
      }
      cells <- sort(unique(selected[seq_len(min(length(selected), max_cells))]))
    } else {
      cells <- sort(sample(cells, max_cells))
    }
  }
  assay <- cfg$assay %||% .fa_default_assay(object)
  counts <- .fa_matrix(object, "counts", assay)
  if (is.null(counts)) return(NULL)
  counts <- counts[, cells, drop = FALSE]
  max_features <- as.integer(cfg$max_features %||% 3000L)
  if (nrow(counts) > max_features) {
    variable <- tryCatch(SeuratObject::VariableFeatures(object[[assay]]), error = function(e) character())
    variable <- intersect(variable, rownames(counts))
    if (length(variable) >= min(500L, max_features)) {
      features <- utils::head(variable, max_features)
    } else {
      detected <- Matrix::rowSums(counts != 0)
      features <- rownames(counts)[order(detected, decreasing = TRUE)[seq_len(max_features)]]
    }
    counts <- counts[features, , drop = FALSE]
  }
  new_cds <- .fa_pkg_fun("monocle3", "new_cell_data_set")
  if (is.null(new_cds)) return(NULL)
  meta <- .seurat_metadata(object)[colnames(counts), , drop = FALSE]
  symbols <- .fa_feature_symbols(object, assay)
  gene_meta <- data.frame(
    gene_short_name = unname(symbols[rownames(counts)]),
    row.names = rownames(counts), stringsAsFactors = FALSE
  )
  cds <- new_cds(counts, cell_metadata = meta, gene_metadata = gene_meta)
  list(cds = cds, cells = colnames(counts), features = rownames(counts))
}

.fa_reduced_dimension <- function(cds, name = "UMAP") {
  reduced_dims <- .fa_pkg_fun("SingleCellExperiment", "reducedDims")
  if (is.null(reduced_dims)) return(NULL)
  dimensions <- tryCatch(reduced_dims(cds), error = function(e) NULL)
  if (is.null(dimensions) || !name %in% names(dimensions)) return(NULL)
  dimensions[[name]]
}

.fa_resolve_trajectory_root <- function(object, built, meta, cfg, graph = NULL, root_column = NULL) {
  output <- list(
    supplied = FALSE, valid = FALSE, type = NULL,
    root_cells = character(), root_pr_nodes = character(),
    message = NULL
  )
  fail <- function(message, type = NULL) {
    output$supplied <- TRUE
    output$valid <- FALSE
    output$type <- type
    output$message <- message
    output
  }
  succeed <- function(type, cells = character(), nodes = character(), message = NULL) {
    output$supplied <- TRUE
    output$valid <- TRUE
    output$type <- type
    output$root_cells <- unique(as.character(cells))
    output$root_pr_nodes <- unique(as.character(nodes))
    output$message <- message
    output
  }

  root <- cfg$root
  if (is.null(root)) {
    legacy_cells <- as.character(cfg$root_cells %||% character())
    if (length(legacy_cells)) root <- list(type = "cells", cells = legacy_cells)
  }
  if (is.null(root) && !is.null(cfg$root_cluster)) {
    root <- list(
      type = "metadata_group",
      column = cfg$root_column %||% root_column,
      values = cfg$root_cluster
    )
  }
  if (is.null(root)) return(output)
  output$supplied <- TRUE
  if ((is.character(root) || is.numeric(root)) && length(root)) {
    root <- list(type = "principal_node", value = root)
  }
  if (!is.list(root) || !length(root)) {
    return(fail("The trajectory root definition must be a principal node or a named list describing cells, markers, or a metadata group."))
  }

  nested_metadata <- root$metadata
  if (is.list(nested_metadata)) root <- .fa_merge_list(root, nested_metadata)
  type <- as.character(root$type %||% "")
  type <- if (length(type) && !is.na(type[[1L]])) type[[1L]] else ""
  type <- tolower(gsub("[ -]+", "_", type))
  if (!nzchar(type)) {
    if (!is.null(root$cells) || !is.null(root$cell_ids) || !is.null(root$root_cells)) {
      type <- "cells"
    } else if (!is.null(root$markers) || !is.null(root$marker_genes)) {
      type <- "markers"
    } else if (!is.null(root$column) || !is.null(root$metadata_column) || !is.null(root$group)) {
      type <- "metadata_group"
    } else if (!is.null(root$principal_node) || !is.null(root$node) || !is.null(root$value)) {
      type <- "principal_node"
    }
  }

  if (type %in% c("principal_node", "principal_nodes", "node", "nodes", "root_pr_nodes")) {
    nodes <- root$principal_node %||% root$nodes %||% root$node %||% root$value
    nodes <- unique(trimws(as.character(nodes %||% character())))
    nodes <- nodes[!is.na(nodes) & nzchar(nodes)]
    if (!length(nodes)) return(fail("No principal-graph node was supplied in trajectory_root.", type))
    graph_nodes <- character()
    if (!is.null(graph) && .fa_pkg_available("igraph")) {
      graph_to_table <- .fa_pkg_fun("igraph", "as_data_frame")
      vertices <- if (!is.null(graph_to_table)) {
        tryCatch(graph_to_table(graph, what = "vertices"), error = function(e) NULL)
      } else {
        NULL
      }
      if (!is.null(vertices) && "name" %in% names(vertices)) graph_nodes <- as.character(vertices$name)
    }
    if (length(graph_nodes)) {
      missing <- setdiff(nodes, graph_nodes)
      if (length(missing)) {
        return(fail(
          paste0("Unknown principal-graph root node(s): ", paste(missing, collapse = ", "),
                 ". Select node IDs from the exported principal_graph_vertices table."),
          type
        ))
      }
    }
    return(succeed(
      "principal_node", nodes = nodes,
      message = paste0("Principal-graph node(s) explicitly selected: ", paste(nodes, collapse = ", "), ".")
    ))
  }

  if (type %in% c("cell", "cells", "cell_ids", "root_cells")) {
    cells <- root$cells %||% root$cell_ids %||% root$root_cells %||% root$value
    cells <- unique(as.character(cells %||% character()))
    cells <- intersect(cells[!is.na(cells) & nzchar(cells)], built$cells)
    if (!length(cells)) {
      return(fail("None of the explicitly supplied trajectory root cells exists in the analyzed cell set.", type))
    }
    return(succeed(
      "cells", cells = cells,
      message = sprintf("%s explicitly supplied cell(s) define the trajectory root.", length(cells))
    ))
  }

  if (type %in% c("metadata", "metadata_group", "group", "label", "cluster")) {
    column <- root$column %||% root$metadata_column %||% root$root_column %||% root_column
    values <- root$values %||% root$groups %||% root$group %||% root$labels %||% root$label %||% root$value
    column <- as.character(column %||% "")
    values <- unique(as.character(values %||% character()))
    values <- values[!is.na(values) & nzchar(values)]
    if (length(column) != 1L || !nzchar(column) || !column %in% names(meta)) {
      return(fail("The metadata-group trajectory root refers to a missing metadata column.", type))
    }
    if (!length(values)) return(fail("The metadata-group trajectory root has no selected value.", type))
    labels <- as.character(meta[[column]])
    cells <- rownames(meta)[!is.na(labels) & labels %in% values]
    cells <- intersect(cells, built$cells)
    if (!length(cells)) {
      return(fail(
        paste0("No analyzed cells matched trajectory root value(s) ", paste(values, collapse = ", "),
               " in metadata column '", column, "'."),
        type
      ))
    }
    return(succeed(
      "metadata_group", cells = cells,
      message = paste0("Metadata column '", column, "' and value(s) ", paste(values, collapse = ", "),
                       " explicitly selected ", length(cells), " root cell(s).")
    ))
  }

  if (type %in% c("marker", "markers", "marker_genes", "genes")) {
    requested <- root$markers %||% root$marker_genes %||% root$genes %||% root$value
    requested <- unique(as.character(requested %||% character()))
    requested <- requested[!is.na(requested) & nzchar(requested)]
    if (!length(requested)) return(fail("The marker-based trajectory root has no marker genes.", type))
    assay <- root$assay %||% cfg$assay %||% .fa_default_assay(object)
    expression <- .fa_matrix(object, "data", assay)
    layer <- "normalized data"
    if (is.null(expression)) {
      expression <- .fa_matrix(object, "counts", assay)
      layer <- "counts"
    }
    if (is.null(expression)) return(fail("Marker-based trajectory rooting requires an expression matrix.", type))
    cells <- intersect(built$cells, colnames(expression))
    symbols <- .fa_feature_symbols(object, assay)
    features <- .fa_match_features(requested, symbols)
    features <- intersect(features, rownames(expression))
    if (!length(features)) {
      return(fail("None of the explicitly supplied trajectory-root markers matched assay features.", type))
    }
    scores <- Matrix::colMeans(expression[features, cells, drop = FALSE], na.rm = TRUE)
    finite <- is.finite(scores)
    if (!any(finite) || length(unique(as.numeric(scores[finite]))) < 2L) {
      return(fail("The supplied root markers did not produce a variable cell score in the analyzed data.", type))
    }
    fraction <- suppressWarnings(as.numeric(root$top_fraction %||% root$fraction %||% 0.05))
    if (!is.finite(fraction) || fraction <= 0 || fraction > 1) fraction <- 0.05
    n_root <- suppressWarnings(as.integer(root$n_cells %||% root$top_n %||%
                                           max(10L, ceiling(length(scores) * fraction))))
    if (!is.finite(n_root) || n_root < 1L) n_root <- max(1L, ceiling(length(scores) * fraction))
    n_root <- min(n_root, sum(finite))
    direction <- as.character(root$direction %||% "high")
    direction <- if (length(direction) && !is.na(direction[[1L]])) tolower(direction[[1L]]) else "high"
    decreasing <- !direction %in% c("low", "lowest", "decreasing")
    ranked <- names(sort(scores[finite], decreasing = decreasing, na.last = NA))
    cells <- utils::head(ranked, n_root)
    return(succeed(
      "markers", cells = cells,
      message = paste0(
        length(cells), " cell(s) with ", if (decreasing) "highest" else "lowest",
        " mean ", layer, " across ", length(features),
        " explicitly supplied marker(s) define the root."
      )
    ))
  }

  fail(
    paste0("Unsupported trajectory root type '", type,
           "'. Use principal_node, cells, markers, or metadata_group."),
    type
  )
}

.fa_module_pseudotime <- function(object, output, cfg, seed, verbose) {
  if (!.fa_pkg_available("monocle3")) {
    return(.fa_result(object, "skipped", "dependency_missing",
                      "monocle3 is not installed; trajectory geometry was not computed."))
  }
  built <- .fa_build_cds(object, cfg, seed)
  if (is.null(built)) {
    return(.fa_result(object, "skipped", "counts_or_constructor_missing",
                      "A sparse counts matrix and monocle3::new_cell_data_set are required."))
  }
  cds <- built$cds
  preprocess <- .fa_pkg_fun("monocle3", "preprocess_cds")
  reduce_dimension <- .fa_pkg_fun("monocle3", "reduce_dimension")
  cluster_cells <- .fa_pkg_fun("monocle3", "cluster_cells")
  learn_graph <- .fa_pkg_fun("monocle3", "learn_graph")
  order_cells <- .fa_pkg_fun("monocle3", "order_cells")
  pseudotime_fun <- .fa_pkg_fun("monocle3", "pseudotime")
  plot_cells <- .fa_pkg_fun("monocle3", "plot_cells")
  graph_test <- .fa_pkg_fun("monocle3", "graph_test")
  required <- list(preprocess, reduce_dimension, cluster_cells, learn_graph)
  if (any(vapply(required, is.null, logical(1)))) {
    return(.fa_result(object, "failed", "monocle3_api_missing",
                      "The installed monocle3 does not expose the required trajectory functions."))
  }
  num_dim <- min(as.integer(cfg$num_dim %||% 30L), nrow(cds) - 1L, ncol(cds) - 1L)
  cds <- preprocess(cds, num_dim = max(2L, num_dim), method = "PCA", verbose = FALSE)
  cds <- reduce_dimension(
    cds, reduction_method = "UMAP",
    preprocess_method = "PCA", umap.fast_sgd = isTRUE(cfg$umap_fast_sgd),
    cores = as.integer(cfg$cores %||% 1L), verbose = FALSE
  )
  cds <- cluster_cells(cds, reduction_method = "UMAP", verbose = FALSE)
  cds <- learn_graph(cds, use_partition = isTRUE(cfg$use_partition), close_loop = isTRUE(cfg$close_loop))
  artifacts <- list()
  coordinates <- .fa_reduced_dimension(cds, "UMAP")
  if (!is.null(coordinates)) {
    coordinate_table <- data.frame(cell = rownames(coordinates), coordinates, row.names = NULL, check.names = FALSE)
    artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
      coordinate_table, output, "pseudotime", "trajectory_umap", "Trajectory UMAP geometry",
      "One row per trajectory cell. Coordinates describe geometry only and do not imply temporal direction.",
      "cell"
    )
  }
  principal_graph_fun <- .fa_pkg_fun("monocle3", "principal_graph")
  graph <- if (!is.null(principal_graph_fun)) tryCatch(principal_graph_fun(cds)[["UMAP"]], error = function(e) NULL) else NULL
  if (!is.null(graph) && .fa_pkg_available("igraph")) {
    graph_to_table <- .fa_pkg_fun("igraph", "as_data_frame")
    edges <- tryCatch(graph_to_table(graph, what = "edges"), error = function(e) NULL)
    vertices <- tryCatch(graph_to_table(graph, what = "vertices"), error = function(e) NULL)
    if (!is.null(edges)) artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
      edges, output, "pseudotime", "principal_graph_edges", "Principal graph edges",
      "One row per undirected principal-graph edge. Vertex IDs describe learned geometry and carry no root direction.",
      "graph edge", list(from = "Source vertex identifier for storage only; the graph is undirected before rooting.",
                         to = "Target vertex identifier for storage only; the graph is undirected before rooting.")
    )
    if (!is.null(vertices)) artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
      vertices, output, "pseudotime", "principal_graph_vertices", "Principal graph vertices",
      "One row per learned principal-graph vertex.", "graph vertex"
    )
  }
  meta <- .seurat_metadata(object)[built$cells, , drop = FALSE]
  root_column <- cfg$root_column %||% .fa_annotation_column(object, cfg) %||% .fa_cluster_column(object, cfg)
  candidate_table <- if (!is.null(root_column) && root_column %in% names(meta)) {
    labels <- as.character(meta[[root_column]])
    rows <- split(seq_along(labels), labels)
    do.call(rbind, lapply(names(rows), function(label) {
      index <- rows[[label]]
      data.frame(
        root_column = root_column,
        candidate_label = label,
        n_cells = length(index),
        example_cell = rownames(meta)[index[[1L]]],
        median_umap_1 = if (!is.null(coordinates)) stats::median(coordinates[rownames(meta)[index], 1L], na.rm = TRUE) else NA_real_,
        status = "candidate_only_requires_biological_selection",
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(
      root_column = NA_character_, candidate_label = NA_character_, n_cells = nrow(meta),
      example_cell = rownames(meta)[1L], median_umap_1 = NA_real_,
      status = "supply_root_cells_or_root_cluster", stringsAsFactors = FALSE
    )
  }
  rownames(candidate_table) <- NULL
  artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
    candidate_table, output, "pseudotime", "root_candidates", "Trajectory root candidates",
    "One row per existing annotation/cluster group that could be considered as a root. These are candidates only; their order and UMAP position are not biological evidence of an origin.",
    "root candidate", list(
      candidate_label = "Existing annotation or cluster label; the package does not select it automatically.",
      median_umap_1 = "Descriptive coordinate only; it is not used as an automatic temporal root.",
      status = "Explicit reminder that biological root selection is still required."
    )
  )
  if (!is.null(graph_test) && isTRUE(cfg$run_graph_test %||% TRUE)) {
    graph_result <- tryCatch(
      graph_test(cds, neighbor_graph = "principal_graph", cores = as.integer(cfg$cores %||% 1L), verbose = FALSE),
      error = function(e) NULL
    )
    if (!is.null(graph_result)) {
      graph_table <- data.frame(feature = rownames(graph_result), as.data.frame(graph_result), row.names = NULL, check.names = FALSE)
      artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
        graph_table, output, "pseudotime", "graph_test", "Trajectory graph-associated features",
        "One row per tested feature. Moran's I tests association with the learned undirected principal graph and does not require or imply a pseudotime root.",
        "feature", list(morans_I = "Moran's I spatial autocorrelation statistic on the principal graph.",
                        q_value = "Multiple-testing adjusted graph-test P value.")
      )
    }
  }
  root <- .fa_resolve_trajectory_root(object, built, meta, cfg, graph, root_column)
  if (!isTRUE(root$valid)) {
    root_issue <- root$message %||% "No biological trajectory root was supplied."
    artifacts[[length(artifacts) + 1L]] <- .fa_write_rds_artifact(
      cds, output, "pseudotime", "trajectory_geometry_unrooted", "Unrooted monocle3 trajectory",
      paste0("monocle3 cell_data_set with reduced coordinates and learned principal graph. It is intentionally unordered. ", root_issue)
    )
    if (!is.null(plot_cells)) {
      geometry_plot <- tryCatch(plot_cells(
        cds, color_cells_by = if (!is.null(root_column) && root_column %in% names(meta)) root_column else "cluster",
        show_trajectory_graph = TRUE, label_groups_by_cluster = FALSE,
        label_leaves = FALSE, label_branch_points = FALSE
      ), error = function(e) NULL)
      if (!is.null(geometry_plot)) artifacts <- c(artifacts, .fa_write_plot_artifacts(
        geometry_plot, output, "pseudotime", "trajectory_geometry_unrooted", "Unrooted trajectory geometry",
        "Learned monocle3 geometry without pseudotime direction. A biological root must be supplied before ordering cells."
      ))
    }
    return(.fa_result(
      object, "needs_input", if (isTRUE(root$supplied)) "trajectory_root_invalid" else "trajectory_root_missing",
      if (isTRUE(root$supplied)) {
        paste0("Trajectory geometry was exported, but pseudotime was not assigned because the supplied root is invalid: ", root_issue)
      } else {
        "Trajectory geometry and root candidates were exported, but pseudotime was not assigned because no biological root was supplied."
      },
      "monocle3", artifacts,
      list(cells = length(built$cells), features = length(built$features), root_column = root_column,
           root_type = root$type, root_issue = root_issue)
    ))
  }
  if (is.null(order_cells) || is.null(pseudotime_fun)) {
    return(.fa_result(object, "partial", "ordering_api_missing",
                      "Trajectory geometry was generated, but the installed monocle3 cannot order cells.",
                      "monocle3", artifacts))
  }
  ordered <- tryCatch(
    if (length(root$root_pr_nodes)) {
      order_cells(cds, reduction_method = "UMAP", root_pr_nodes = root$root_pr_nodes)
    } else {
      order_cells(cds, reduction_method = "UMAP", root_cells = root$root_cells)
    },
    error = function(e) e
  )
  if (inherits(ordered, "error")) {
    return(.fa_result(
      object, "needs_input", "trajectory_root_invalid",
      paste0("The supplied trajectory root could not be applied by monocle3: ", conditionMessage(ordered)),
      "monocle3", artifacts,
      list(root_type = root$type, root_message = root$message)
    ))
  }
  cds <- ordered
  pseudotime <- tryCatch(pseudotime_fun(cds), error = function(e) NULL)
  if (is.null(pseudotime)) {
    return(.fa_result(object, "partial", "pseudotime_unavailable",
                      "A root was supplied, but monocle3 did not return finite pseudotime values.",
                      "monocle3", artifacts))
  }
  pseudotime_table <- data.frame(
    cell = names(pseudotime), pseudotime = as.numeric(pseudotime),
    is_root_cell = names(pseudotime) %in% root$root_cells,
    root_definition = root$type,
    stringsAsFactors = FALSE
  )
  artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
    pseudotime_table, output, "pseudotime", "cell_pseudotime", "Cell pseudotime",
    "One row per trajectory cell. Direction is defined only by the explicitly supplied root cells, marker rule, metadata group, or principal-graph node.",
    "cell", list(pseudotime = "monocle3 pseudotime distance from the explicitly supplied root.",
                 is_root_cell = "TRUE for cells passed directly to monocle3 as roots; FALSE when a principal node defines the root.",
                 root_definition = "Explicit root definition type used for this ordering."),
    units = list(pseudotime = "monocle3 pseudotime units")
  )
  artifacts[[length(artifacts) + 1L]] <- .fa_write_rds_artifact(
    cds, output, "pseudotime", "trajectory_ordered", "Ordered monocle3 trajectory",
    "monocle3 cell_data_set including the explicitly rooted pseudotime ordering."
  )
  if (!is.null(plot_cells)) {
    ordered_plot <- tryCatch(plot_cells(
      cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE,
      label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE
    ), error = function(e) NULL)
    if (!is.null(ordered_plot)) artifacts <- c(artifacts, .fa_write_plot_artifacts(
      ordered_plot, output, "pseudotime", "trajectory_pseudotime", "Rooted pseudotime",
      "monocle3 trajectory colored by pseudotime defined from the explicitly supplied root."
    ))
  }
  .fa_result(
    object, "completed", "explicitly_rooted_pseudotime_completed",
    paste0("Pseudotime was assigned using an explicit ", root$type, " root. ", root$message),
    "monocle3", artifacts,
    list(cells = length(built$cells), features = length(built$features),
         root_column = root_column, root_type = root$type,
         n_root_cells = length(root$root_cells), root_pr_nodes = root$root_pr_nodes)
  )
}

# Cell-cell communication ----------------------------------------------------

.fa_matrix_to_long <- function(matrix, value_name) {
  if (is.null(matrix) || !length(matrix)) return(data.frame())
  table <- as.data.frame(as.table(matrix), stringsAsFactors = FALSE)
  names(table) <- c("source", "target", value_name)
  table
}

.fa_module_communication <- function(object, output, cfg, seed, verbose, species) {
  if (!.fa_pkg_available("CellChat")) {
    return(.fa_result(object, "skipped", "dependency_missing",
                      "CellChat is not installed; communication analysis was not run."))
  }
  annotation_column <- .fa_annotation_column(object, cfg)
  if (is.null(annotation_column)) {
    return(.fa_result(
      object, "needs_input", "annotation_missing",
      "CellChat requires an existing cell annotation column. Cluster IDs were not silently relabeled as cell types."
    ))
  }
  resources <- .fa_species_resources(species, cfg)
  database <- cfg$database
  database_name <- cfg$database_name %||% resources$cellchat
  if (is.null(database) && !is.null(database_name)) database <- .fa_pkg_object("CellChat", database_name)
  if (is.null(database)) {
    return(.fa_result(
      object, "needs_input", "species_database_missing",
      paste0("No CellChat database was supplied for species '", species,
             "'. The package will not substitute the human or mouse database for another species.")
    ))
  }
  assay <- cfg$assay %||% .fa_default_assay(object)
  data <- .fa_matrix(object, "data", assay)
  if (is.null(data)) {
    normalize <- .fa_pkg_fun("Seurat", "NormalizeData")
    if (is.null(normalize)) {
      return(.fa_result(object, "skipped", "normalized_data_missing",
                        "CellChat requires normalized expression data."))
    }
    object <- normalize(object, assay = assay, verbose = FALSE)
    data <- .fa_matrix(object, "data", assay)
  }
  meta <- .seurat_metadata(object)[colnames(data), , drop = FALSE]
  labels <- as.character(meta[[annotation_column]])
  valid <- !is.na(labels) & nzchar(labels)
  minimum_cells <- as.integer(cfg$min_cells %||% 10L)
  valid_labels <- names(which(table(labels[valid]) >= minimum_cells))
  valid <- valid & labels %in% valid_labels
  if (length(valid_labels) < 2L) {
    return(.fa_result(
      object, "skipped", "fewer_than_two_valid_groups",
      sprintf("CellChat needs at least two annotation groups with %s or more cells.", minimum_cells)
    ))
  }
  max_per_group <- as.integer(cfg$max_cells_per_group %||% 1000L)
  set.seed(as.integer(seed))
  selected <- unlist(lapply(valid_labels, function(label) {
    candidates <- which(valid & labels == label)
    if (length(candidates) > max_per_group) sample(candidates, max_per_group) else candidates
  }), use.names = FALSE)
  selected <- sort(unique(selected))
  data <- data[, selected, drop = FALSE]
  meta <- meta[selected, , drop = FALSE]
  symbols <- .fa_feature_symbols(object, assay)
  orgdb <- .fa_orgdb(resources, cfg)
  mapping <- .fa_feature_mapping(
    rownames(data), symbols, orgdb,
    preferred_keytype = cfg$feature_keytype %||% resources$feature_keytype,
    symbol_column = cfg$symbol_column %||% resources$symbol_column %||% "SYMBOL"
  )
  mapped <- as.character(mapping$SYMBOL[match(rownames(data), mapping$feature)])
  feature_keep <- !is.na(mapped) & nzchar(mapped) & !duplicated(toupper(mapped))
  data <- data[feature_keep, , drop = FALSE]
  rownames(data) <- mapped[feature_keep]
  create_cellchat <- .fa_pkg_fun("CellChat", "createCellChat")
  functions <- lapply(
    c("subsetData", "identifyOverExpressedGenes", "identifyOverExpressedInteractions",
      "computeCommunProb", "filterCommunication", "computeCommunProbPathway",
      "aggregateNet", "subsetCommunication"),
    function(name) .fa_pkg_fun("CellChat", name)
  )
  names(functions) <- c("subsetData", "identifyOverExpressedGenes", "identifyOverExpressedInteractions",
                        "computeCommunProb", "filterCommunication", "computeCommunProbPathway",
                        "aggregateNet", "subsetCommunication")
  if (is.null(create_cellchat) || any(vapply(functions, is.null, logical(1)))) {
    return(.fa_result(object, "failed", "cellchat_api_missing",
                      "The installed CellChat does not expose the required analysis functions."))
  }
  cellchat_meta <- data.frame(
    annotation = as.character(meta[[annotation_column]]),
    row.names = rownames(meta), stringsAsFactors = FALSE
  )
  cellchat <- create_cellchat(object = data, meta = cellchat_meta, group.by = "annotation")
  methods::slot(cellchat, "DB") <- database
  pipeline_error <- NULL
  cellchat <- tryCatch({
    value <- functions$subsetData(cellchat)
    value <- functions$identifyOverExpressedGenes(value)
    value <- functions$identifyOverExpressedInteractions(value)
    value <- functions$computeCommunProb(
      value, type = cfg$probability_method %||% "triMean",
      raw.use = isTRUE(cfg$raw_use), population.size = isTRUE(cfg$population_size)
    )
    value <- functions$filterCommunication(value, min.cells = minimum_cells)
    value <- functions$computeCommunProbPathway(value)
    functions$aggregateNet(value)
  }, error = function(e) {
    pipeline_error <<- conditionMessage(e)
    cellchat
  })
  if (!is.null(pipeline_error)) {
    failed_groups <- as.data.frame(table(annotation = cellchat_meta$annotation), stringsAsFactors = FALSE)
    names(failed_groups)[2L] <- "n_cells_analyzed"
    failure_artifact <- .fa_write_table_artifact(
      failed_groups, output, "communication", "cellchat_groups", "CellChat analyzed groups",
      "One row per annotation group supplied to CellChat after deterministic downsampling. No interaction table was produced.",
      "annotation group"
    )
    no_signal <- grepl("no rows|subscript out of bounds|dimension|0 interactions", pipeline_error, ignore.case = TRUE)
    return(.fa_result(
      object, if (no_signal) "skipped" else "failed",
      if (no_signal) "no_detectable_communication" else "cellchat_failed",
      if (no_signal) {
        paste0("CellChat found no analyzable ligand-receptor signal after filtering: ", pipeline_error)
      } else {
        paste0("CellChat failed: ", pipeline_error)
      },
      "CellChat", list(failure_artifact),
      list(species = species, database = database_name %||% "user supplied",
           annotation_column = annotation_column, cells = ncol(data), groups = valid_labels)
    ))
  }
  interactions <- tryCatch(functions$subsetCommunication(cellchat), error = function(e) data.frame())
  pathways <- tryCatch(functions$subsetCommunication(cellchat, slot.name = "netP"), error = function(e) data.frame())
  net <- .slot_or_null(cellchat, "net") %||% list()
  count_table <- .fa_matrix_to_long(net$count, "interaction_count")
  weight_table <- .fa_matrix_to_long(net$weight, "interaction_weight")
  artifacts <- list()
  if (nrow(interactions)) artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
    interactions, output, "communication", "ligand_receptor_interactions", "CellChat ligand-receptor interactions",
    "One row per inferred source-target ligand-receptor record returned by CellChat.",
    "communication record", list(source = "Sending annotation group.", target = "Receiving annotation group.",
                                 ligand = "Ligand gene.", receptor = "Receptor gene.",
                                 prob = "CellChat communication probability.", pval = "CellChat permutation P value.")
  )
  if (nrow(pathways)) artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
    pathways, output, "communication", "pathway_communication", "CellChat pathway communication",
    "One row per source-target signaling-pathway record returned by CellChat.",
    "pathway communication record"
  )
  if (nrow(count_table)) artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
    count_table, output, "communication", "network_counts", "CellChat network counts",
    "One row per sender-receiver pair; interaction_count is the number of retained inferred interactions.",
    "sender-receiver pair"
  )
  if (nrow(weight_table)) artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
    weight_table, output, "communication", "network_weights", "CellChat network weights",
    "One row per sender-receiver pair; interaction_weight is the aggregated CellChat probability weight.",
    "sender-receiver pair"
  )
  group_table <- as.data.frame(table(annotation = cellchat_meta$annotation), stringsAsFactors = FALSE)
  names(group_table)[2L] <- "n_cells_analyzed"
  group_table$original_n_cells <- as.integer(table(factor(labels[valid], levels = group_table$annotation)))
  artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
    group_table, output, "communication", "cellchat_groups", "CellChat analyzed groups",
    "One row per cell annotation group. n_cells_analyzed reflects deterministic per-group downsampling; original_n_cells is the eligible count before downsampling.",
    "annotation group"
  )
  artifacts[[length(artifacts) + 1L]] <- .fa_write_rds_artifact(
    cellchat, output, "communication", "cellchat_object", "CellChat analysis object",
    "CellChat object containing the database, inferred interactions, pathways, and aggregated network."
  )
  if (nrow(weight_table)) {
    weight_table$source <- factor(weight_table$source)
    weight_table$target <- factor(weight_table$target)
    heatmap <- ggplot2::ggplot(weight_table, ggplot2::aes(x = target, y = source, fill = interaction_weight)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#0072B2") +
      ggplot2::labs(x = "Receiver", y = "Sender", fill = "Weight",
                    title = "Aggregated cell-cell communication") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    artifacts <- c(artifacts, .fa_write_plot_artifacts(
      heatmap, output, "communication", "network_weights", "Communication network weights",
      "CellChat aggregated interaction weights for every sender-receiver annotation pair."
    ))
  }
  .fa_result(
    object, "completed", "cellchat_completed",
    sprintf("CellChat analyzed %s cells from %s annotation groups after deterministic per-group downsampling.", ncol(data), length(valid_labels)),
    "CellChat", artifacts,
    list(species = species, database = database_name %||% "user supplied",
         annotation_column = annotation_column, cells = ncol(data), groups = valid_labels)
  )
}

# Cell-cycle scoring ----------------------------------------------------------

.fa_cell_cycle_genes <- function(species, cfg) {
  if (length(cfg$s_genes) && length(cfg$g2m_genes)) {
    mapping <- rbind(
      data.frame(phase_set = "S", source_human_symbol = NA_character_,
                 target_symbol = as.character(cfg$s_genes), stringsAsFactors = FALSE),
      data.frame(phase_set = "G2M", source_human_symbol = NA_character_,
                 target_symbol = as.character(cfg$g2m_genes), stringsAsFactors = FALSE)
    )
    return(list(s.genes = as.character(cfg$s_genes), g2m.genes = as.character(cfg$g2m_genes),
                source = "user", mapping = mapping))
  }
  resources <- .fa_species_resources(species, cfg)
  strategy <- resources$cell_cycle_strategy %||% "user_supplied"
  if (identical(strategy, "user_supplied")) return(NULL)
  genes <- .fa_pkg_object("Seurat", "cc.genes.updated.2019")
  if (is.null(genes) || is.null(genes$s.genes) || is.null(genes$g2m.genes)) return(NULL)
  if (identical(species, "human")) {
    mapping <- rbind(
      data.frame(phase_set = "S", source_human_symbol = genes$s.genes,
                 target_symbol = genes$s.genes, stringsAsFactors = FALSE),
      data.frame(phase_set = "G2M", source_human_symbol = genes$g2m.genes,
                 target_symbol = genes$g2m.genes, stringsAsFactors = FALSE)
    )
    return(list(s.genes = genes$s.genes, g2m.genes = genes$g2m.genes,
                source = "Seurat cc.genes.updated.2019; native human symbols", mapping = mapping))
  }
  if (.fa_pkg_available("babelgene") && !is.null(resources$scientific_name)) {
    orthologs <- .fa_pkg_fun("babelgene", "orthologs")
    requested <- unique(c(genes$s.genes, genes$g2m.genes))
    mapped <- tryCatch(
      orthologs(
        requested, species = resources$scientific_name, human = TRUE,
        min_support = as.integer(cfg$ortholog_min_support %||% 3L), top = TRUE
      ),
      error = function(e) NULL
    )
    if (!is.null(mapped) && nrow(mapped) && all(c("human_symbol", "symbol") %in% names(mapped))) {
      s_target <- as.character(mapped$symbol[match(genes$s.genes, mapped$human_symbol)])
      g_target <- as.character(mapped$symbol[match(genes$g2m.genes, mapped$human_symbol)])
      mapping <- rbind(
        data.frame(phase_set = "S", source_human_symbol = genes$s.genes,
                   target_symbol = s_target, stringsAsFactors = FALSE),
        data.frame(phase_set = "G2M", source_human_symbol = genes$g2m.genes,
                   target_symbol = g_target, stringsAsFactors = FALSE)
      )
      mapping <- mapping[!is.na(mapping$target_symbol) & nzchar(mapping$target_symbol), , drop = FALSE]
      return(list(
        s.genes = unique(s_target[!is.na(s_target) & nzchar(s_target)]),
        g2m.genes = unique(g_target[!is.na(g_target) & nzchar(g_target)]),
        source = paste0(
          "Seurat cc.genes.updated.2019 mapped by babelgene to ",
          resources$scientific_name, " (minimum ortholog support ",
          as.integer(cfg$ortholog_min_support %||% 3L), ")"
        ),
        mapping = mapping
      ))
    }
  }
  if (identical(species, "mouse")) {
    mapping <- rbind(
      data.frame(phase_set = "S", source_human_symbol = genes$s.genes,
                 target_symbol = genes$s.genes, stringsAsFactors = FALSE),
      data.frame(phase_set = "G2M", source_human_symbol = genes$g2m.genes,
                 target_symbol = genes$g2m.genes, stringsAsFactors = FALSE)
    )
    return(list(
      s.genes = genes$s.genes, g2m.genes = genes$g2m.genes,
      source = "Seurat cc.genes.updated.2019; mouse case-insensitive symbol fallback (babelgene unavailable)",
      mapping = mapping
    ))
  }
  NULL
}

.fa_module_cell_cycle <- function(object, output, cfg, seed, verbose, species) {
  genes <- .fa_cell_cycle_genes(species, cfg)
  if (is.null(genes)) {
    strategy <- .fa_species_resources(species, cfg)$cell_cycle_strategy %||% "user_supplied"
    needs_user <- identical(strategy, "user_supplied")
    return(.fa_result(
      object, if (needs_user) "needs_input" else "skipped",
      if (needs_user) "cell_cycle_genes_missing" else "cell_cycle_ortholog_resource_unavailable",
      if (needs_user) {
        paste0("No validated cell-cycle gene set is configured for species '", species,
               "'. Supply s_genes and g2m_genes; another species will not be substituted.")
      } else {
        paste0("The registered cell-cycle strategy for '", species,
               "' could not load its Seurat/babelgene resource. Install the optional dependency or supply s_genes and g2m_genes.")
      }
    ))
  }
  assay <- cfg$assay %||% .fa_default_assay(object)
  symbols_from_object <- .fa_feature_symbols(object, assay)
  resources <- .fa_species_resources(species, cfg)
  orgdb <- .fa_orgdb(resources, cfg)
  mapping <- .fa_feature_mapping(
    names(symbols_from_object), symbols_from_object, orgdb,
    preferred_keytype = cfg$feature_keytype %||% resources$feature_keytype,
    symbol_column = cfg$symbol_column %||% resources$symbol_column %||% "SYMBOL"
  )
  symbols <- stats::setNames(as.character(mapping$SYMBOL), mapping$feature)
  s_features <- .fa_match_features(genes$s.genes, symbols)
  g2m_features <- .fa_match_features(genes$g2m.genes, symbols)
  if (length(s_features) < as.integer(cfg$min_genes_per_set %||% 5L) ||
      length(g2m_features) < as.integer(cfg$min_genes_per_set %||% 5L)) {
    return(.fa_result(
      object, "needs_input", "insufficient_cell_cycle_gene_overlap",
      sprintf("Only %s S-phase and %s G2M genes matched the assay; supply species-appropriate genes or feature symbols.",
              length(s_features), length(g2m_features))
    ))
  }
  meta <- .seurat_metadata(object)
  score_columns <- c("S.Score", "G2M.Score", "Phase")
  used_existing <- all(score_columns %in% names(meta)) && !isTRUE(cfg$force)
  if (!used_existing) {
    scoring <- .fa_pkg_fun("Seurat", "CellCycleScoring")
    if (is.null(scoring)) {
      return(.fa_result(object, "skipped", "seurat_scoring_missing",
                        "Seurat::CellCycleScoring is not available."))
    }
    if (is.null(.fa_matrix(object, "data", assay))) {
      normalize <- .fa_pkg_fun("Seurat", "NormalizeData")
      if (is.null(normalize)) {
        return(.fa_result(object, "skipped", "normalized_data_missing",
                          "Cell-cycle scoring requires normalized expression data."))
      }
      object <- normalize(object, assay = assay, verbose = FALSE)
    }
    control_features <- as.integer(cfg$ctrl %||% min(100L, max(1L, floor(nrow(object) / 24L) - 1L)))
    object <- scoring(
      object = object, s.features = s_features, g2m.features = g2m_features,
      ctrl = control_features, set.ident = FALSE
    )
    meta <- .seurat_metadata(object)
  }
  sample_column <- .fa_sample_column(object, cfg)
  annotation_column <- .fa_annotation_column(object, cfg)
  scores <- data.frame(
    cell = rownames(meta), S.Score = as.numeric(meta$S.Score),
    G2M.Score = as.numeric(meta$G2M.Score), Phase = as.character(meta$Phase),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!is.null(sample_column)) scores$sample <- as.character(meta[[sample_column]])
  if (!is.null(annotation_column)) scores$annotation <- as.character(meta[[annotation_column]])
  matched_genes <- rbind(
    data.frame(phase_set = "S", feature = s_features, symbol = unname(symbols[s_features]), stringsAsFactors = FALSE),
    data.frame(phase_set = "G2M", feature = g2m_features, symbol = unname(symbols[g2m_features]), stringsAsFactors = FALSE)
  )
  artifacts <- list(
    .fa_write_table_artifact(
      scores, output, "cell_cycle", "cell_cycle_scores", "Cell-cycle scores",
      "One row per cell. S.Score, G2M.Score, and Phase are reused from the RDS when already present, otherwise calculated by Seurat from a species-appropriate matched gene set.",
      "cell", list(S.Score = "Seurat S-phase module score.", G2M.Score = "Seurat G2M-phase module score.",
                   Phase = "Phase assigned by comparing S and G2M scores."),
      units = list(S.Score = "Seurat module score", G2M.Score = "Seurat module score")
    ),
    .fa_write_table_artifact(
      matched_genes, output, "cell_cycle", "matched_genes", "Matched cell-cycle genes",
      "One row per feature used for S or G2M scoring. Feature is the exact assay row name and symbol is the matched biological symbol.",
      "matched feature"
    )
  )
  if (!is.null(genes$mapping) && nrow(genes$mapping)) {
    artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
      genes$mapping, output, "cell_cycle", "cell_cycle_gene_resource",
      "Cell-cycle gene resource mapping",
      "One row per source cell-cycle gene and phase set. Non-human target symbols are species-matched orthologs when babelgene is available; this table records the exact mapping used.",
      "cell-cycle source gene",
      list(source_human_symbol = "Human source symbol from Seurat cc.genes.updated.2019; missing for user-supplied target genes.",
           target_symbol = "Species-matched symbol requested for feature matching.")
    )
  }
  if (!is.null(sample_column)) {
    valid <- !is.na(scores$sample) & nzchar(scores$sample) & !is.na(scores$Phase) & nzchar(scores$Phase)
    composition <- as.data.frame(table(sample = scores$sample[valid], phase = scores$Phase[valid]), stringsAsFactors = FALSE)
    composition <- composition[composition$Freq > 0, , drop = FALSE]
    names(composition)[3L] <- "n_cells"
    totals <- tapply(composition$n_cells, composition$sample, sum)
    composition$fraction <- composition$n_cells / totals[composition$sample]
    artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
      composition, output, "cell_cycle", "sample_phase_composition", "Sample cell-cycle composition",
      "One row per observed sample-by-phase combination; fractions sum to one within each sample.",
      "sample-phase", list(fraction = "Phase cell count divided by all phase-assigned cells in the sample."),
      units = list(fraction = "fraction")
    )
  }
  reductions <- .seurat_reductions(object)
  reductions <- unique(c(reductions[grepl("umap", reductions, ignore.case = TRUE)], reductions))
  if (length(reductions)) {
    embedding <- .embedding_table(object, reductions[[1L]])
    if (!is.null(embedding) && ncol(embedding) >= 3L) {
      plotted <- merge(embedding[c("cell", names(embedding)[2:3])], scores, by = "cell", sort = FALSE)
      names(plotted)[2:3] <- c("dimension_1", "dimension_2")
      max_plot <- as.integer(cfg$max_plot_cells %||% 100000L)
      if (nrow(plotted) > max_plot) {
        set.seed(as.integer(seed))
        plotted <- plotted[sort(sample(seq_len(nrow(plotted)), max_plot)), , drop = FALSE]
      }
      plotted$Phase <- factor(plotted$Phase)
      phase_plot <- ggplot2::ggplot(plotted, ggplot2::aes(x = dimension_1, y = dimension_2, color = Phase)) +
        ggplot2::geom_point(size = 0.42, alpha = 0.78) +
        ggplot2::scale_color_manual(values = .fa_palette(nlevels(plotted$Phase))) +
        ggplot2::coord_equal() +
        ggplot2::labs(x = names(embedding)[2L], y = names(embedding)[3L], title = "Cell-cycle phase") +
        ggplot2::theme_bw(base_size = 12) + ggplot2::theme(panel.grid = ggplot2::element_blank())
      artifacts <- c(artifacts, .fa_write_plot_artifacts(
        phase_plot, output, "cell_cycle", "phase_embedding", "Cell-cycle phase embedding",
        "Stored reduction coordinates colored by Seurat cell-cycle phase."
      ))
      for (score in c("S.Score", "G2M.Score")) {
        plotted$score_value <- plotted[[score]]
        score_plot <- ggplot2::ggplot(plotted, ggplot2::aes(x = dimension_1, y = dimension_2, color = score_value)) +
          ggplot2::geom_point(size = 0.42, alpha = 0.78) +
          ggplot2::scale_color_gradient2(low = "#0072B2", mid = "#F7F7F7", high = "#D55E00", midpoint = 0) +
          ggplot2::coord_equal() +
          ggplot2::labs(x = names(embedding)[2L], y = names(embedding)[3L], color = score, title = score) +
          ggplot2::theme_bw(base_size = 12) + ggplot2::theme(panel.grid = ggplot2::element_blank())
        artifacts <- c(artifacts, .fa_write_plot_artifacts(
          score_plot, output, "cell_cycle", paste0(score, "_embedding"), paste(score, "embedding"),
          paste("Stored reduction coordinates colored by", score, ".")
        ))
      }
    }
  }
  .fa_result(
    object, "completed", if (used_existing) "existing_cell_cycle_preserved" else "cell_cycle_scored",
    if (used_existing) "Existing S.Score, G2M.Score, and Phase columns were exported without recalculation." else
      "Cell-cycle scores were calculated from species-appropriate, case-insensitively matched features.",
    "Seurat::CellCycleScoring", artifacts,
    list(species = species, gene_source = genes$source, s_genes_matched = length(s_features),
         g2m_genes_matched = length(g2m_features), reused_existing = used_existing)
  )
}

# Transcription-factor expression --------------------------------------------

.fa_default_tf_genes <- function(species) {
  if (!species %in% c("human", "mouse", "rat", "pig", "cattle", "dog", "macaque")) return(character())
  c(
    "AR", "ARID1A", "ATF3", "BATF", "BCL6", "CEBPA", "CEBPB", "CREB1", "CTCF",
    "E2F1", "E2F2", "E2F3", "EBF1", "ELF1", "ELK1", "ERG", "ESR1", "ETS1",
    "FOS", "FOSB", "FOXA1", "FOXP1", "FOXP3", "GATA1", "GATA2", "GATA3",
    "HIF1A", "IKZF1", "IRF1", "IRF3", "IRF4", "IRF7", "IRF8", "JUN", "JUNB",
    "JUND", "KLF2", "KLF4", "KLF6", "MAF", "MAFB", "MAX", "MEF2C", "MITF",
    "MYC", "MYCN", "NANOG", "NFE2L2", "NFATC1", "NFATC2", "NFKB1", "NFKB2",
    "NR3C1", "PAX5", "POU2F1", "POU2F2", "PPARG", "PRDM1", "RELA", "RELB",
    "REST", "RUNX1", "RUNX2", "RUNX3", "SMAD2", "SMAD3", "SMAD4", "SOX2",
    "SPI1", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
    "TBX21", "TCF7", "TFEB", "TGIF1", "TP53", "XBP1", "YY1", "ZEB1"
  )
}

.fa_orgdb_tf_genes <- function(species, cfg) {
  resources <- .fa_species_resources(species, cfg)
  orgdb <- .fa_orgdb(resources, cfg)
  if (is.null(orgdb) || !.fa_pkg_available("AnnotationDbi")) return(character())
  select_fun <- .fa_pkg_fun("AnnotationDbi", "select")
  keytypes_fun <- .fa_pkg_fun("AnnotationDbi", "keytypes")
  columns_fun <- .fa_pkg_fun("AnnotationDbi", "columns")
  keytypes <- tryCatch(keytypes_fun(orgdb), error = function(e) character())
  columns <- tryCatch(columns_fun(orgdb), error = function(e) character())
  if (!"GOALL" %in% keytypes || !"SYMBOL" %in% columns) return(character())
  requested_columns <- intersect(c("SYMBOL", "ENTREZID", "EVIDENCEALL"), columns)
  table <- tryCatch(
    suppressMessages(select_fun(
      orgdb, keys = "GO:0003700", keytype = "GOALL", columns = requested_columns
    )),
    error = function(e) NULL
  )
  if (is.null(table) || !nrow(table) || !"SYMBOL" %in% names(table)) return(character())
  unique(as.character(table$SYMBOL[!is.na(table$SYMBOL) & nzchar(table$SYMBOL)]))
}

.fa_tf_catalog <- function(species, cfg) {
  genes <- cfg$tf_genes
  source <- NULL
  if (is.data.frame(genes)) {
    column <- intersect(c("symbol", "gene", "tf", "TF"), names(genes))
    genes <- if (length(column)) genes[[column[[1L]]]] else character()
    source <- "user data frame"
  }
  if (is.character(genes) && length(genes) == 1L && file.exists(genes)) {
    table <- tryCatch(utils::read.delim(genes, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(table)) {
      column <- intersect(c("symbol", "gene", "tf", "TF"), names(table))
      genes <- if (length(column)) table[[column[[1L]]]] else table[[1L]]
      source <- paste0("user file: ", basename(cfg$tf_genes))
    }
  }
  genes <- unique(as.character(genes %||% character()))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes)) {
    attr(genes, "source") <- source %||% "user vector"
    return(genes)
  }
  genes <- .fa_orgdb_tf_genes(species, cfg)
  if (length(genes)) {
    attr(genes, "source") <- paste0(
      "species OrgDb GOALL:GO:0003700 (", .fa_species_resources(species, cfg)$orgdb, ")"
    )
    return(genes)
  }
  genes <- .fa_default_tf_genes(species)
  if (length(genes)) attr(genes, "source") <- "built-in curated TF seed (fallback, not a complete catalog)"
  genes
}

.fa_module_tf <- function(object, output, cfg, seed, verbose, species) {
  catalog <- .fa_tf_catalog(species, cfg)
  catalog_source <- attr(catalog, "source") %||% "unknown"
  if (!length(catalog)) {
    strategy <- .fa_species_resources(species, cfg)$tf_catalog_strategy %||% "user_supplied"
    needs_user <- identical(strategy, "user_supplied")
    return(.fa_result(
      object, if (needs_user) "needs_input" else "skipped",
      if (needs_user) "tf_catalog_missing" else "species_orgdb_tf_catalog_unavailable",
      if (needs_user) {
        paste0("No validated TF catalog is configured for species '", species,
               "'. Supply tf_genes; another species catalog will not be substituted.")
      } else {
        paste0("The registered OrgDb TF strategy for '", species,
               "' is unavailable in this R library. Install the matching OrgDb or supply tf_genes.")
      }
    ))
  }
  assay <- cfg$assay %||% .fa_default_assay(object)
  counts <- .fa_matrix(object, "counts", assay)
  if (is.null(counts)) {
    return(.fa_result(object, "skipped", "counts_missing", "TF expression summaries require a counts layer."))
  }
  symbols_from_object <- .fa_feature_symbols(object, assay)
  resources <- .fa_species_resources(species, cfg)
  orgdb <- .fa_orgdb(resources, cfg)
  mapping <- .fa_feature_mapping(
    names(symbols_from_object), symbols_from_object, orgdb,
    preferred_keytype = cfg$feature_keytype %||% resources$feature_keytype,
    symbol_column = cfg$symbol_column %||% resources$symbol_column %||% "SYMBOL"
  )
  symbols <- stats::setNames(as.character(mapping$SYMBOL), mapping$feature)
  features <- .fa_match_features(catalog, symbols)
  if (length(features) < as.integer(cfg$min_tf_genes %||% 5L)) {
    return(.fa_result(
      object, "needs_input", "insufficient_tf_overlap",
      sprintf("Only %s TF genes matched assay features; supply a species- and identifier-appropriate catalog.", length(features))
    ))
  }
  matched <- data.frame(
    requested_symbol = catalog,
    matched_feature = names(symbols)[match(toupper(catalog), toupper(symbols))],
    matched_symbol = unname(symbols[match(toupper(catalog), toupper(symbols))]),
    stringsAsFactors = FALSE
  )
  matched$matched <- !is.na(matched$matched_feature)
  artifacts <- list(.fa_write_table_artifact(
    matched, output, "tf", "tf_catalog_mapping", "TF catalog mapping",
    paste0("One row per requested TF symbol. matched_feature is the exact assay row used for expression summaries. Catalog source: ", catalog_source, "."),
    "TF catalog entry", list(matched = "TRUE when the TF symbol matched an assay feature case-insensitively.")
  ))
  tf_counts <- counts[features, , drop = FALSE]
  meta <- .seurat_metadata(object)[colnames(tf_counts), , drop = FALSE]
  sample_column <- .fa_sample_column(object, cfg)
  group_column <- .fa_group_column(object, cfg)
  annotation_column <- .fa_annotation_column(object, cfg)
  schemes <- list()
  if (!is.null(sample_column)) schemes$sample <- as.character(meta[[sample_column]])
  if (!is.null(group_column)) schemes$group <- as.character(meta[[group_column]])
  if (!is.null(annotation_column)) schemes$annotation <- as.character(meta[[annotation_column]])
  if (!is.null(sample_column) && !is.null(annotation_column)) {
    schemes$sample_by_annotation <- interaction(
      as.character(meta[[sample_column]]), as.character(meta[[annotation_column]]),
      sep = " | ", drop = TRUE
    )
  }
  if (!length(schemes)) schemes$all_cells <- rep("all_cells", ncol(tf_counts))
  averages <- lapply(names(schemes), function(scheme) {
    values <- .fa_long_average(tf_counts, schemes[[scheme]], feature_symbols = symbols[features])
    values$grouping_scheme <- scheme
    values
  })
  averages <- do.call(rbind, averages)
  averages <- averages[c("grouping_scheme", "group", "feature", "symbol", "average_log_normalized_expression")]
  artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
    averages, output, "tf", "tf_average_expression", "TF average expression",
    "One row per TF, grouping scheme, and group. Values are log1p counts-per-10,000 from sparse group aggregation; they describe TF expression, not regulatory activity.",
    "TF-group", list(
      grouping_scheme = "Metadata level used to aggregate cells: sample, group, annotation, or sample_by_annotation.",
      group = "Observed value or value combination within the grouping scheme.",
      average_log_normalized_expression = "log1p of group-aggregated counts normalized to 10,000 total counts."
    ), units = list(average_log_normalized_expression = "log1p counts per 10,000")
  )
  heat_scheme <- if ("sample" %in% averages$grouping_scheme) "sample" else names(schemes)[[1L]]
  heat <- averages[averages$grouping_scheme == heat_scheme, , drop = FALSE]
  variability <- tapply(heat$average_log_normalized_expression, heat$symbol, stats::sd, na.rm = TRUE)
  top_symbols <- names(sort(variability, decreasing = TRUE, na.last = NA))
  top_symbols <- utils::head(top_symbols, as.integer(cfg$heatmap_top_n %||% 40L))
  heat <- heat[heat$symbol %in% top_symbols, , drop = FALSE]
  if (nrow(heat)) {
    z_values <- stats::ave(heat$average_log_normalized_expression, heat$symbol, FUN = function(x) {
      value <- as.numeric(scale(x))
      value[!is.finite(value)] <- 0
      value
    })
    heat$z_expression <- z_values
    heat$symbol <- factor(heat$symbol, levels = rev(top_symbols))
    tf_plot <- ggplot2::ggplot(heat, ggplot2::aes(x = group, y = symbol, fill = z_expression)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "#0072B2", mid = "#F7F7F7", high = "#D55E00", midpoint = 0) +
      ggplot2::labs(x = tools::toTitleCase(heat_scheme), y = "Transcription factor",
                    fill = "Row z-score", title = "Variable TF expression") +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    artifacts <- c(artifacts, .fa_write_plot_artifacts(
      tf_plot, output, "tf", "tf_expression_heatmap", "TF expression heatmap",
      "Row-standardized average TF expression for the most variable TFs. This is expression, not inferred TF activity.",
      width = max(8, 0.35 * length(unique(heat$group)) + 4), height = max(6, 0.18 * length(top_symbols) + 2)
    ))
  }
  fallback_catalog <- grepl("fallback|curated TF seed", catalog_source, ignore.case = TRUE)
  .fa_result(
    object, if (fallback_catalog) "partial" else "completed",
    if (fallback_catalog) "curated_tf_seed_fallback" else "tf_expression_completed",
    sprintf(
      "Expression summaries were generated for %s matched TF genes from %s; no TF activity was inferred.",
      length(features), catalog_source
    ),
    "sparse pseudobulk expression", artifacts,
    list(species = species, requested_tfs = length(catalog), matched_tfs = length(features),
         catalog_source = catalog_source, grouping_schemes = names(schemes))
  )
}

# inferCNV --------------------------------------------------------------------

.fa_standardize_gene_order <- function(x) {
  if (is.null(x)) return(NULL)
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  aliases <- list(
    gene = c("gene", "gene_id", "gene_name", "symbol", "V1"),
    chromosome = c("chromosome", "chr", "seqnames", "seqname", "V2"),
    start = c("start", "gene_start", "V3"),
    stop = c("stop", "end", "gene_end", "V4")
  )
  selected <- vapply(aliases, function(candidates) {
    hit <- candidates[candidates %in% names(x)]
    if (length(hit)) hit[[1L]] else NA_character_
  }, character(1))
  if (anyNA(selected) && ncol(x) >= 4L) selected <- names(x)[seq_len(4L)]
  if (anyNA(selected)) return(NULL)
  output <- data.frame(
    gene = as.character(x[[selected[["gene"]]]]),
    chromosome = as.character(x[[selected[["chromosome"]]]]),
    start = suppressWarnings(as.numeric(x[[selected[["start"]]]])),
    stop = suppressWarnings(as.numeric(x[[selected[["stop"]]]])),
    stringsAsFactors = FALSE
  )
  output <- output[
    !is.na(output$gene) & nzchar(output$gene) & !is.na(output$chromosome) & nzchar(output$chromosome) &
      is.finite(output$start) & is.finite(output$stop), , drop = FALSE
  ]
  output <- output[!duplicated(output$gene), , drop = FALSE]
  output[order(output$chromosome, output$start, output$stop), , drop = FALSE]
}

.fa_txdb_object <- function(value) {
  if (is.null(value)) return(NULL)
  if (!is.character(value)) return(value)
  value <- trimws(as.character(value[[1L]]))
  if (!nzchar(value)) return(NULL)
  if (grepl("::", value, fixed = TRUE)) {
    fields <- strsplit(value, "::", fixed = TRUE)[[1L]]
    if (length(fields) == 2L) return(.fa_pkg_object(fields[[1L]], fields[[2L]]))
  }
  if (!.fa_pkg_available(value)) return(NULL)
  object <- .fa_pkg_object(value, value)
  if (!is.null(object)) return(object)
  exports <- tryCatch(getNamespaceExports(value), error = function(e) character())
  candidates <- exports[grepl("^TxDb", exports)]
  if (length(candidates)) .fa_pkg_object(value, candidates[[1L]]) else NULL
}

.fa_first_character <- function(x) {
  vapply(seq_along(x), function(i) {
    value <- tryCatch(as.character(x[[i]]), error = function(e) character())
    value <- value[!is.na(value) & nzchar(value)]
    if (length(value)) value[[1L]] else NA_character_
  }, character(1))
}

.fa_txdb_gene_order <- function(txdb, orgdb = NULL) {
  txdb <- .fa_txdb_object(txdb)
  genes <- .fa_pkg_fun("GenomicFeatures", "genes")
  if (is.null(txdb) || is.null(genes)) return(NULL)
  ranges <- tryCatch(
    suppressMessages(genes(txdb, columns = "gene_id", single.strand.genes.only = TRUE)),
    error = function(e) tryCatch(
      suppressMessages(genes(txdb, single.strand.genes.only = TRUE)),
      error = function(e2) NULL
    )
  )
  if (is.null(ranges) || !length(ranges)) return(NULL)
  table <- tryCatch(as.data.frame(ranges, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(table) || !nrow(table)) return(NULL)
  ids <- if ("gene_id" %in% names(table)) {
    .fa_first_character(table$gene_id)
  } else {
    as.character(names(ranges))
  }
  if (length(ids) != nrow(table) || all(is.na(ids) | !nzchar(ids))) {
    ids <- as.character(names(ranges))
  }
  if (length(ids) != nrow(table)) ids <- rep(NA_character_, nrow(table))
  symbols <- ids
  if (!is.null(orgdb) && .fa_pkg_available("AnnotationDbi")) {
    select <- .fa_pkg_fun("AnnotationDbi", "select")
    keytypes <- .fa_pkg_fun("AnnotationDbi", "keytypes")
    available <- if (!is.null(keytypes)) tryCatch(keytypes(orgdb), error = function(e) character()) else character()
    likely_entrez <- mean(grepl("^[0-9]+$", ids[!is.na(ids)])) >= 0.8
    keytype <- if (likely_entrez && "ENTREZID" %in% available) {
      "ENTREZID"
    } else if ("ENSEMBL" %in% available) {
      "ENSEMBL"
    } else if ("ENTREZID" %in% available) {
      "ENTREZID"
    } else {
      NULL
    }
    if (!is.null(select) && !is.null(keytype)) {
      mapping <- tryCatch(
        suppressMessages(select(orgdb, keys = unique(ids[!is.na(ids) & nzchar(ids)]),
                                keytype = keytype, columns = "SYMBOL")),
        error = function(e) NULL
      )
      if (!is.null(mapping) && all(c(keytype, "SYMBOL") %in% names(mapping))) {
        mapping <- mapping[!is.na(mapping$SYMBOL) & nzchar(mapping$SYMBOL), c(keytype, "SYMBOL"), drop = FALSE]
        mapping <- mapping[!duplicated(mapping[[keytype]]), , drop = FALSE]
        matched <- match(ids, as.character(mapping[[keytype]]))
        use <- !is.na(matched)
        symbols[use] <- as.character(mapping$SYMBOL[matched[use]])
      }
    }
  }
  chromosome <- if ("seqnames" %in% names(table)) table$seqnames else if ("seqname" %in% names(table)) table$seqname else NA
  stop_column <- if ("end" %in% names(table)) table$end else if ("stop" %in% names(table)) table$stop else NA
  .fa_standardize_gene_order(data.frame(
    gene = symbols,
    chromosome = as.character(chromosome),
    start = table$start,
    stop = stop_column,
    stringsAsFactors = FALSE
  ))
}

.fa_gene_order <- function(cfg) {
  order <- cfg$gene_order
  if (is.character(order) && length(order) == 1L && file.exists(order)) {
    table <- tryCatch(utils::read.delim(order, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
                      error = function(e) NULL)
    standardized <- .fa_standardize_gene_order(table)
    if (is.null(standardized)) {
      table <- tryCatch(utils::read.delim(order, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE),
                        error = function(e) NULL)
      standardized <- .fa_standardize_gene_order(table)
    }
    return(standardized)
  }
  standardized <- .fa_standardize_gene_order(order)
  if (!is.null(standardized)) return(standardized)
  gtf <- cfg$gtf
  if (is.character(gtf) && length(gtf) == 1L && file.exists(gtf) && .fa_pkg_available("rtracklayer")) {
    import <- .fa_pkg_fun("rtracklayer", "import")
    ranges <- tryCatch(import(gtf), error = function(e) NULL)
    if (!is.null(ranges)) {
      table <- as.data.frame(ranges)
      if ("type" %in% names(table)) table <- table[table$type == "gene", , drop = FALSE]
      gene_column <- intersect(c("gene_name", "gene_id", "Name", "ID"), names(table))
      if (length(gene_column)) {
        return(.fa_standardize_gene_order(data.frame(
          gene = table[[gene_column[[1L]]]], chromosome = table$seqnames,
          start = table$start, stop = table$end, stringsAsFactors = FALSE
        )))
      }
    }
  }
  if (!is.null(cfg$txdb)) {
    resources <- list(orgdb = cfg$orgdb)
    orgdb <- .fa_orgdb(resources, cfg)
    standardized <- .fa_txdb_gene_order(cfg$txdb, orgdb)
    if (!is.null(standardized) && nrow(standardized)) return(standardized)
  }
  NULL
}

.fa_write_infercnv_inputs <- function(annotation, gene_order, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  annotation_path <- file.path(directory, "infercnv_annotations.tsv")
  order_path <- file.path(directory, "infercnv_gene_order.tsv")
  utils::write.table(annotation, annotation_path, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE, na = "")
  utils::write.table(gene_order[c("gene", "chromosome", "start", "stop")], order_path,
                     sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")
  list(annotation = annotation_path, gene_order = order_path)
}

.fa_module_cnv <- function(object, output, cfg, seed, verbose, species) {
  reference_groups <- unique(as.character(cfg$reference_groups %||% character()))
  reference_groups <- reference_groups[!is.na(reference_groups) & nzchar(reference_groups)]
  if (!length(reference_groups)) {
    return(.fa_result(
      object, "needs_input", "cnv_reference_missing",
      "inferCNV was not run because no biological reference_groups were supplied. The package never guesses a normal reference population."
    ))
  }
  if (!.fa_pkg_available("infercnv")) {
    return(.fa_result(object, "skipped", "dependency_missing", "infercnv is not installed."))
  }
  annotation_column <- cfg$annotation_column %||% .fa_annotation_column(object, cfg)
  if (is.null(annotation_column)) {
    return(.fa_result(
      object, "needs_input", "cnv_annotation_missing",
      "inferCNV requires an existing annotation column containing the explicitly selected reference groups."
    ))
  }
  resources <- .fa_species_resources(species, cfg)
  for (field in c("gene_order", "gtf", "txdb", "orgdb")) {
    if (is.null(cfg[[field]]) && !is.null(resources[[field]])) cfg[[field]] <- resources[[field]]
  }
  gene_order <- .fa_gene_order(cfg)
  if (is.null(gene_order) || nrow(gene_order) < as.integer(cfg$min_ordered_genes %||% 100L)) {
    return(.fa_result(
      object, "needs_input", "gene_order_missing",
      paste0("inferCNV requires a species/build-matched gene_order data frame/file, GTF, or installed TxDb resource. ",
             "No usable gene order could be resolved for species '", species, "'.")
    ))
  }
  assay <- cfg$assay %||% .fa_default_assay(object)
  counts <- .fa_matrix(object, "counts", assay)
  if (is.null(counts)) {
    return(.fa_result(object, "skipped", "counts_missing", "inferCNV requires raw counts."))
  }
  symbols <- .fa_feature_symbols(object, assay)
  feature_to_gene <- rownames(counts)
  direct_overlap <- mean(feature_to_gene %in% gene_order$gene)
  if (direct_overlap < 0.2) feature_to_gene <- unname(symbols[rownames(counts)])
  gene_index <- match(feature_to_gene, gene_order$gene)
  keep_features <- !is.na(gene_index) & !duplicated(feature_to_gene)
  if (sum(keep_features) < as.integer(cfg$min_ordered_genes %||% 100L)) {
    return(.fa_result(
      object, "needs_input", "insufficient_gene_order_overlap",
      sprintf("Only %s assay features matched the supplied gene order; check species, genome build, and identifier type.", sum(keep_features))
    ))
  }
  counts <- counts[keep_features, , drop = FALSE]
  rownames(counts) <- feature_to_gene[keep_features]
  gene_order <- gene_order[match(rownames(counts), gene_order$gene), , drop = FALSE]
  ordered <- order(gene_order$chromosome, gene_order$start, gene_order$stop)
  counts <- counts[ordered, , drop = FALSE]
  gene_order <- gene_order[ordered, , drop = FALSE]
  meta <- .seurat_metadata(object)[colnames(counts), , drop = FALSE]
  labels <- as.character(meta[[annotation_column]])
  if (!all(reference_groups %in% unique(labels))) {
    missing <- setdiff(reference_groups, unique(labels))
    return(.fa_result(
      object, "needs_input", "cnv_reference_not_found",
      paste0("Reference groups are absent from '", annotation_column, "': ", paste(missing, collapse = ", "), ".")
    ))
  }
  include_groups <- unique(c(reference_groups, as.character(cfg$include_groups %||% unique(labels))))
  valid <- !is.na(labels) & nzchar(labels) & labels %in% include_groups
  max_per_group <- as.integer(cfg$max_cells_per_group %||% 300L)
  set.seed(as.integer(seed))
  selected <- unlist(lapply(sort(unique(labels[valid])), function(label) {
    candidates <- which(valid & labels == label)
    if (length(candidates) > max_per_group) sample(candidates, max_per_group) else candidates
  }), use.names = FALSE)
  selected <- sort(unique(selected))
  counts <- counts[, selected, drop = FALSE]
  meta <- meta[selected, , drop = FALSE]
  labels <- as.character(meta[[annotation_column]])
  annotation <- data.frame(cell = colnames(counts), group = labels, stringsAsFactors = FALSE)
  directory <- .fa_module_directory(output, "cnv")
  inputs <- .fa_write_infercnv_inputs(annotation, gene_order, directory)
  artifacts <- list(
    .fa_write_table_artifact(
      annotation, output, "cnv", "cell_annotations", "inferCNV cell annotations",
      "One row per analyzed cell. group is copied from the selected existing annotation column; reference status comes only from explicit reference_groups.",
      "cell", list(group = paste0("Value from metadata column '", annotation_column, "'."))
    ),
    .fa_write_table_artifact(
      gene_order, output, "cnv", "gene_order", "inferCNV gene order",
      "One row per matched expression feature in genomic order. Coordinates come only from the explicitly supplied gene-order/GTF resource.",
      "gene", list(chromosome = "Chromosome or contig in the supplied genome build.",
                   start = "Genomic start coordinate.", stop = "Genomic end coordinate."),
      units = list(start = "base pairs", stop = "base pairs")
    ),
    .fa_artifact_record(
      "cnv", "infercnv_input", inputs$annotation, output, "inferCNV annotation input",
      "Headerless tab-separated input passed directly to infercnv::CreateInfercnvObject.",
      rows = nrow(annotation), columns = 2L, row_unit = "cell", column_unit = "field",
      column_dictionary = list(cell = "Cell name.", group = "Existing annotation group."), complete = TRUE
    ),
    .fa_artifact_record(
      "cnv", "infercnv_input", inputs$gene_order, output, "inferCNV gene-order input",
      "Headerless tab-separated gene-order input passed directly to infercnv::CreateInfercnvObject.",
      rows = nrow(gene_order), columns = 4L, row_unit = "gene", column_unit = "field",
      column_dictionary = list(gene = "Matched gene name.", chromosome = "Genome-build contig.",
                               start = "Start coordinate.", stop = "End coordinate."), complete = TRUE
    )
  )
  create_infercnv <- .fa_pkg_fun("infercnv", "CreateInfercnvObject")
  run_infercnv <- .fa_pkg_fun("infercnv", "run")
  if (is.null(create_infercnv) || is.null(run_infercnv)) {
    return(.fa_result(object, "failed", "infercnv_api_missing",
                      "The installed infercnv does not expose CreateInfercnvObject and run.",
                      artifacts = artifacts))
  }
  infer_object <- create_infercnv(
    raw_counts_matrix = counts,
    annotations_file = inputs$annotation,
    delim = "\t",
    gene_order_file = inputs$gene_order,
    ref_group_names = reference_groups,
    min_max_counts_per_cell = cfg$min_max_counts_per_cell %||% c(100, Inf)
  )
  run_directory <- file.path(directory, "infercnv_run")
  dir.create(run_directory, recursive = TRUE, showWarnings = FALSE)
  infer_object <- run_infercnv(
    infer_object,
    cutoff = as.numeric(cfg$cutoff %||% 0.1),
    out_dir = run_directory,
    cluster_by_groups = isTRUE(cfg$cluster_by_groups %||% TRUE),
    denoise = isTRUE(cfg$denoise),
    HMM = isTRUE(cfg$HMM),
    num_threads = as.integer(cfg$num_threads %||% 1L),
    no_plot = isTRUE(cfg$no_plot %||% FALSE)
  )
  artifacts[[length(artifacts) + 1L]] <- .fa_write_rds_artifact(
    infer_object, output, "cnv", "infercnv_object", "inferCNV analysis object",
    "infercnv object generated from a bounded per-group cell subset and an explicitly supplied biological reference and gene order."
  )
  expression <- .slot_or_null(infer_object, "expr.data")
  if (!is.null(expression) && nrow(expression) && ncol(expression)) {
    baseline <- as.numeric(cfg$signal_baseline %||% 1)
    signal <- vapply(seq_len(ncol(expression)), function(i) {
      mean(abs(as.numeric(expression[, i]) - baseline), na.rm = TRUE)
    }, numeric(1))
    cell_signal <- data.frame(
      cell = colnames(expression),
      annotation = labels[match(colnames(expression), colnames(counts))],
      mean_absolute_cnv_signal = signal,
      reference = labels[match(colnames(expression), colnames(counts))] %in% reference_groups,
      stringsAsFactors = FALSE
    )
    artifacts[[length(artifacts) + 1L]] <- .fa_write_table_artifact(
      cell_signal, output, "cnv", "cell_cnv_signal", "Per-cell inferCNV signal",
      "One row per inferCNV cell. Signal is the mean absolute deviation of the final inferCNV expression profile from the configured neutral baseline; it is a descriptive summary, not a mutation call.",
      "cell", list(mean_absolute_cnv_signal = "Mean absolute deviation from signal_baseline across ordered genes.",
                   reference = "TRUE only for explicitly supplied reference groups."),
      units = list(mean_absolute_cnv_signal = "inferCNV relative-expression units")
    )
    signal_plot <- ggplot2::ggplot(cell_signal, ggplot2::aes(x = annotation, y = mean_absolute_cnv_signal, fill = annotation)) +
      ggplot2::geom_boxplot(outlier.size = 0.35, linewidth = 0.35) +
      ggplot2::scale_fill_manual(values = .fa_palette(length(unique(cell_signal$annotation)))) +
      ggplot2::labs(x = "Annotation", y = "Mean absolute inferCNV signal", title = "Per-cell inferCNV signal") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    artifacts <- c(artifacts, .fa_write_plot_artifacts(
      signal_plot, output, "cnv", "cell_cnv_signal", "Per-cell inferCNV signal",
      "Distribution of descriptive inferCNV signal by existing annotation group; this is not a mutation or malignancy classifier."
    ))
  }
  .fa_result(
    object, "completed", "infercnv_completed",
    sprintf("inferCNV analyzed %s cells and %s ordered genes using explicit reference groups: %s.",
            ncol(counts), nrow(counts), paste(reference_groups, collapse = ", ")),
    "infercnv", artifacts,
    list(species = species, annotation_column = annotation_column,
         reference_groups = reference_groups, cells = ncol(counts), genes = nrow(counts),
         downsampling_max_per_group = max_per_group)
  )
}

# Download/artifact registry and public internal orchestrator -----------------

.fa_artifact_table <- function(artifacts) {
  if (!length(artifacts)) return(data.frame())
  rows <- lapply(artifacts, function(artifact) {
    dictionary <- tryCatch(
      jsonlite::toJSON(artifact$column_dictionary %||% list(), auto_unbox = TRUE, null = "null"),
      error = function(e) "{}"
    )
    units <- tryCatch(
      jsonlite::toJSON(artifact$units %||% list(), auto_unbox = TRUE, null = "null"),
      error = function(e) "{}"
    )
    data.frame(
      artifact_id = artifact$artifact_id %||% NA_character_,
      module = artifact$module %||% NA_character_,
      type = artifact$type %||% NA_character_,
      format = artifact$format %||% NA_character_,
      path = artifact$path %||% NA_character_,
      path_is_relative = isTRUE(artifact$path_is_relative),
      label = artifact$label %||% NA_character_,
      description = artifact$description %||% NA_character_,
      rows = artifact$rows %||% NA_integer_,
      columns = artifact$columns %||% NA_integer_,
      row_unit = artifact$row_unit %||% NA_character_,
      column_unit = artifact$column_unit %||% NA_character_,
      column_dictionary_json = as.character(dictionary),
      units_json = as.character(units),
      bytes = artifact$bytes %||% NA_real_,
      sha256 = artifact$sha256 %||% NA_character_,
      complete = isTRUE(artifact$complete),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

.fa_module_table <- function(modules) {
  if (!length(modules)) return(data.frame())
  rows <- lapply(modules, function(module) {
    data.frame(
      module = module$id %||% NA_character_,
      requested = isTRUE(module$requested),
      eligible = isTRUE(module$eligible),
      status = module$status %||% NA_character_,
      reason_code = module$reason_code %||% NA_character_,
      message = module$message %||% NA_character_,
      engine = module$engine %||% NA_character_,
      engine_version = module$engine_version %||% NA_character_,
      elapsed_seconds = module$elapsed_seconds %||% NA_real_,
      cells_used = module$cells_used %||% NA_integer_,
      features_used = module$features_used %||% NA_integer_,
      artifact_count = length(module$artifact_ids %||% character()),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.fa_module_downloads <- function(object, output, cfg, seed, verbose, modules, artifacts) {
  artifact_table <- .fa_artifact_table(artifacts)
  module_table <- .fa_module_table(modules)
  output_artifacts <- list()
  output_artifacts[[length(output_artifacts) + 1L]] <- .fa_write_table_artifact(
    artifact_table, output, "downloads", "artifact_index", "Full-analysis artifact index",
    "One row per artifact created before this index. Every table artifact links to a JSON-encoded column dictionary and records size, checksum, completeness, and row/column units.",
    "artifact", list(
      artifact_id = "Stable identifier unique within this analysis run.",
      module = "Canonical report module that created the artifact.",
      path = "Path relative to the user-selected output directory when path_is_relative is TRUE.",
      column_dictionary_json = "JSON object explaining every column of a table artifact.",
      sha256 = "SHA-256 checksum when the optional digest package is available.",
      complete = "TRUE when the artifact writer completed successfully."
    )
  )
  output_artifacts[[length(output_artifacts) + 1L]] <- .fa_write_table_artifact(
    module_table, output, "downloads", "module_status", "Full-analysis module status",
    "One row per scientific module executed before the downloads registry. skipped and needs_input rows retain their exact reason rather than disappearing from the report.",
    "module", list(
      requested = "Whether the module was enabled in config$modules.",
      eligible = "TRUE when the module completed fully or partially with scientifically valid output.",
      status = "completed, partial, skipped, needs_input, or failed.",
      reason_code = "Stable machine-readable reason for the status.",
      message = "Human-readable explanation suitable for a report chapter."
    )
  )
  .fa_result(
    object, "completed", "artifact_registry_completed",
    sprintf("Registered %s analysis artifacts from %s preceding modules.", nrow(artifact_table), nrow(module_table)),
    "scRDSreport artifact registry", output_artifacts,
    list(artifact_count_before_registry = nrow(artifact_table), module_count_before_registry = nrow(module_table))
  )
}

#' Run the optional full single-cell analysis layer
#'
#' This internal API is intentionally separated from `running()` so the core
#' RDS export path remains lightweight and failure-tolerant.  It always returns
#' all twelve canonical module records, including explicit skipped or
#' needs_input records when a scientifically necessary input is unavailable.
#'
#' @param object A Seurat object with joined Seurat v5 layers.
#' @param output User-selected output directory.
#' @param design Sample design returned by `infer_sample_design()` or its
#'   sample-level design data frame.
#' @param species_info Species provenance list or selected species name.
#' @param config Nested configuration. `config$modules` may be a named logical
#'   vector; per-module options may be placed at `config[[module]]` or
#'   `config$module_options[[module]]`.
#' @param seed One reproducible integer seed.
#' @param verbose Whether to print progress messages.
#' @return A list with `object`, named `modules`, standard `artifacts`, and
#'   captured `warnings`.
.run_full_analysis <- function(object, output, design = NULL, species_info = NULL,
                               config = list(), seed = 11L, verbose = TRUE) {
  if (!.is_seurat(object)) .sc_stop("Full analysis requires a Seurat object.")
  if (!is.list(config)) .sc_stop("Full-analysis config must be a list.")
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
    .sc_stop("Full-analysis seed must be one finite number.")
  }
  output <- normalizePath(output, mustWork = FALSE)
  dir.create(file.path(output, "analysis"), recursive = TRUE, showWarnings = FALSE)
  if (is.list(design) && !is.null(design$design) && !is.null(design$cell_sample)) {
    meta <- .seurat_metadata(object)
    if (!".scRDSreport_sample" %in% names(meta) && length(design$cell_sample) == ncol(object)) {
      object <- .attach_design(object, design)
    }
  }
  species <- .fa_species(species_info, config)
  modules <- list()
  artifacts <- list()
  warnings <- character()
  execute <- function(id, fun) {
    result <- .fa_run_module(id, object, output, config, seed, verbose, fun)
    object <<- result$object
    modules[[id]] <<- result$module
    artifacts <<- c(artifacts, result$artifacts)
    warnings <<- c(warnings, result$warnings)
    invisible(NULL)
  }
  execute("qc", .fa_module_qc)
  execute("reduction", .fa_module_reduction)
  execute("cluster", .fa_module_cluster)
  execute("celltype", function(object, output, cfg, seed, verbose) {
    .fa_module_celltype(object, output, cfg, seed, verbose, species = species)
  })
  execute("differential", .fa_module_differential)
  execute("enrichment", function(object, output, cfg, seed, verbose) {
    .fa_module_enrichment(object, output, cfg, seed, verbose, species = species)
  })
  execute("pseudotime", .fa_module_pseudotime)
  execute("communication", function(object, output, cfg, seed, verbose) {
    .fa_module_communication(object, output, cfg, seed, verbose, species = species)
  })
  execute("cell_cycle", function(object, output, cfg, seed, verbose) {
    .fa_module_cell_cycle(object, output, cfg, seed, verbose, species = species)
  })
  execute("tf", function(object, output, cfg, seed, verbose) {
    .fa_module_tf(object, output, cfg, seed, verbose, species = species)
  })
  execute("cnv", function(object, output, cfg, seed, verbose) {
    .fa_module_cnv(object, output, cfg, seed, verbose, species = species)
  })
  downloads_result <- .fa_run_module(
    "downloads", object, output, config, seed, verbose,
    function(object, output, cfg, seed, verbose) {
      .fa_module_downloads(object, output, cfg, seed, verbose, modules = modules, artifacts = artifacts)
    }
  )
  object <- downloads_result$object
  modules$downloads <- downloads_result$module
  artifacts <- c(artifacts, downloads_result$artifacts)
  warnings <- unique(c(warnings, downloads_result$warnings))
  list(
    object = object,
    modules = modules,
    artifacts = artifacts,
    warnings = warnings,
    species = species,
    schema_version = "2.0"
  )
}
