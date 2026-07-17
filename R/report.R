.render_quarto_report <- function(manifest_path, output, verbose = TRUE) {
  .require_optional("DT", "render tables with CSV and Excel download buttons")
  .require_optional("knitr", "render the Quarto report")
  template <- system.file("quarto", "report.qmd", package = "scRDSreport")
  if (!nzchar(template)) .sc_stop("The installed report template is missing.")
  work_template <- file.path(output, ".report", "report.qmd")
  file.copy(template, work_template, overwrite = TRUE)
  scss <- system.file("quarto", "report.scss", package = "scRDSreport")
  if (nzchar(scss)) file.copy(scss, file.path(output, ".report", "report.scss"), overwrite = TRUE)

  old_locale <- Sys.getlocale("LC_CTYPE")
  old_locale_env <- Sys.getenv(c("LC_ALL", "LANG"), unset = NA_character_)
  locale_candidates <- c("en_US.UTF-8", "en_US.utf8", "zh_CN.UTF-8", "zh_CN.utf8")
  for (candidate in locale_candidates) {
    selected <- suppressWarnings(tryCatch(Sys.setlocale("LC_CTYPE", candidate), error = function(e) ""))
    if (nzchar(selected)) {
      Sys.setenv(LC_ALL = candidate, LANG = candidate)
      break
    }
  }
  on.exit(suppressWarnings(try(Sys.setlocale("LC_CTYPE", old_locale), silent = TRUE)), add = TRUE)
  on.exit({
    for (name in names(old_locale_env)) {
      if (is.na(old_locale_env[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, stats::setNames(list(old_locale_env[[name]]), name))
    }
  }, add = TRUE)

  if (requireNamespace("quarto", quietly = TRUE)) {
    .sc_message(verbose, "Rendering Quarto HTML report...")
    quarto::quarto_render(
      input = work_template,
      output_file = "report.html",
      execute_params = list(manifest = normalizePath(manifest_path)),
      quiet = !isTRUE(verbose),
      quarto_args = c("--output-dir", output)
    )
  } else {
    quarto_bin <- Sys.which("quarto")
    if (!nzchar(quarto_bin)) {
      .sc_stop("Quarto was not found. Data exports are complete; install Quarto and rerun with overwrite = TRUE, or call running(..., render = FALSE).")
    }
    # Conda's Quarto launcher can be relocated without its bundled tools tree.
    # Prefer sibling `deno` and ../share/quarto when they exist and the user has
    # not already selected explicit paths.
    quarto_prefix <- dirname(dirname(quarto_bin))
    sibling_deno <- file.path(dirname(quarto_bin), "deno")
    sibling_sass <- file.path(dirname(quarto_bin), "sass")
    sibling_esbuild <- file.path(dirname(quarto_bin), "esbuild")
    inferred_share <- file.path(quarto_prefix, "share", "quarto")
    restore_env <- character()
    if (!nzchar(Sys.getenv("DENO_DIR"))) {
      deno_cache <- file.path(output, ".report", "deno-cache")
      dir.create(deno_cache, recursive = TRUE, showWarnings = FALSE)
      Sys.setenv(DENO_DIR = deno_cache)
      restore_env <- c(restore_env, "DENO_DIR")
    }
    if (!nzchar(Sys.getenv("XDG_CACHE_HOME"))) {
      quarto_cache <- file.path(output, ".report", "cache")
      dir.create(quarto_cache, recursive = TRUE, showWarnings = FALSE)
      Sys.setenv(XDG_CACHE_HOME = quarto_cache)
      restore_env <- c(restore_env, "XDG_CACHE_HOME")
    }
    if (!nzchar(Sys.getenv("QUARTO_DENO")) && file.exists(sibling_deno)) {
      Sys.setenv(QUARTO_DENO = sibling_deno)
      restore_env <- c(restore_env, "QUARTO_DENO")
    }
    if (!nzchar(Sys.getenv("QUARTO_SHARE_PATH")) && dir.exists(inferred_share)) {
      Sys.setenv(QUARTO_SHARE_PATH = inferred_share)
      restore_env <- c(restore_env, "QUARTO_SHARE_PATH")
    }
    if (!nzchar(Sys.getenv("QUARTO_DART_SASS")) && file.exists(sibling_sass)) {
      Sys.setenv(QUARTO_DART_SASS = sibling_sass)
      restore_env <- c(restore_env, "QUARTO_DART_SASS")
    }
    if (!nzchar(Sys.getenv("QUARTO_ESBUILD")) && file.exists(sibling_esbuild)) {
      Sys.setenv(QUARTO_ESBUILD = sibling_esbuild)
      restore_env <- c(restore_env, "QUARTO_ESBUILD")
    }
    expected_pandoc <- file.path(dirname(quarto_bin), "tools", "x86_64", "pandoc")
    sibling_pandoc <- file.path(dirname(quarto_bin), "pandoc")
    sibling_js <- file.path(dirname(quarto_bin), "quarto.js")
    sibling_vendor <- file.path(dirname(quarto_bin), "vendor")
    if (!file.exists(expected_pandoc) && file.exists(sibling_pandoc) &&
        file.exists(sibling_js) && dir.exists(sibling_vendor)) {
      compat_bin <- file.path(output, ".report", "quarto-cli", "bin")
      compat_tools <- file.path(compat_bin, "tools", "x86_64")
      dir.create(compat_tools, recursive = TRUE, showWarnings = FALSE)
      file.copy(quarto_bin, file.path(compat_bin, "quarto"), overwrite = TRUE, copy.mode = TRUE)
      file.copy(sibling_js, file.path(compat_bin, "quarto.js"), overwrite = TRUE, copy.mode = TRUE)
      compat_vendor <- file.path(compat_bin, "vendor")
      if (!file.exists(compat_vendor) && !dir.exists(compat_vendor)) {
        file.symlink(sibling_vendor, compat_vendor)
      }
      compat_pandoc <- file.path(compat_tools, "pandoc")
      if (!file.exists(compat_pandoc)) file.symlink(sibling_pandoc, compat_pandoc)
      quarto_bin <- file.path(compat_bin, "quarto")
    }
    if (length(restore_env)) on.exit(Sys.unsetenv(restore_env), add = TRUE)
    params_path <- file.path(output, ".report", "params.yml")
    writeLines(c("manifest:", paste0("  ", normalizePath(manifest_path))), params_path)
    old_working_directory <- setwd(dirname(work_template))
    on.exit(setwd(old_working_directory), add = TRUE)
    args <- c(
      "render", basename(work_template), "--output", "report.html",
      "--output-dir", "..", "--execute-params", basename(params_path)
    )
    status <- system2(quarto_bin, args = args)
    if (!identical(status, 0L)) .sc_stop("Quarto rendering failed with status %s.", status)
  }
  file.path(output, "report.html")
}

#' Build a downloadable single-cell report from an RDS file
#'
#' @param input Path to a single-cell RDS file.
#' @param output Output directory for the report and data files.
#' @param sample_col Optional metadata column containing sample IDs.
#' @param sample_map Optional named sample-to-group vector or data frame with
#'   `sample_id`, `group`, and optionally `replicate`.
#' @param filter_raw_barcodes One of `"auto"`, `"always"`, or `"never"`.
#'   Auto detects large unfiltered droplet matrices before SCP analysis.
#' @param min_features,min_counts,max_features,max_counts Cell-level thresholds
#'   used when raw barcode filtering is applied.
#' @param filter_low_expression_features One of `"auto"`, `"always"`, or
#'   `"never"`. Auto removes features detected in too few cells from large raw
#'   objects before SCP analysis; existing analyzed objects are never changed.
#' @param min_cells_per_feature Minimum cells in which a retained feature must
#'   have a non-zero count.
#' @param analyze One of `"auto"`, `"always"`, or `"never"`. Auto runs SCP
#'   unless both a dimensional reduction and a cluster/annotation column exist.
#' @param integration_method SCP integration method, or `"none"` (default).
#'   Automatic integration is deliberately disabled because sample and batch
#'   are not necessarily equivalent.
#' @param scp_args Named list passed to `SCP::Standard_SCP()` or
#'   `SCP::Integration_SCP()`.
#' @param run_markers Whether to run `SCP::RunDEtest()` for cluster markers
#'   after SCP completes a raw or partial object.
#' @param marker_args Named list passed to `SCP::RunDEtest()`.
#' @param annotation_col Optional existing cell-type/annotation metadata column.
#'   When omitted, common annotation column names are detected automatically.
#'   The package never invents or overwrites cell-type annotations.
#' @param species `"auto"` or one non-empty species name, such as `"human"`,
#'   `"mouse"`, `"rat"`, or `"zebrafish"`. This records provenance only; it
#'   does not force a species-specific annotation workflow. Auto uses stable
#'   feature-ID prefixes before gene symbols.
#' @param analysis_max_cells Maximum cells passed to SCP. `Inf` (default) keeps
#'   all retained cells; a finite value creates a deterministic analysis subset
#'   while the original RDS remains available for download.
#' @param analysis_max_features Maximum features passed to SCP. `Inf` (default)
#'   keeps all retained features; finite values keep the most detected features.
#' @param matrix_layers Expression layers to export in Matrix Market format.
#' @param embed_downloads One of `"auto"`, `"always"`, or `"never"`.
#'   Auto embeds downloadable files into the standalone HTML up to the total
#'   `embed_max_mb` budget, prioritizing annotation and result tables. Always
#'   forces every exported file into the HTML and can create a very large file.
#' @param embed_max_mb Total megabytes available to automatic HTML download
#'   embedding. This limit is ignored when `embed_downloads = "always"`.
#' @param render Whether to render the Quarto HTML report.
#' @param overwrite Whether managed files in a non-empty output directory may
#'   be replaced. Unrelated files are never deleted.
#' @param seed Random seed used by SCP.
#' @param verbose Print progress messages.
#' @return Invisibly, a list containing paths, analysis status, design, and manifest.
#' @export
running <- function(input, output, sample_col = NULL, sample_map = NULL,
                    filter_raw_barcodes = c("auto", "always", "never"),
                    min_features = 200L, min_counts = 0L,
                    max_features = Inf, max_counts = Inf,
                    filter_low_expression_features = c("auto", "always", "never"),
                    min_cells_per_feature = 3L,
                    analyze = c("auto", "always", "never"),
                    integration_method = "none", scp_args = list(),
                    run_markers = TRUE, marker_args = list(),
                    annotation_col = NULL, species = "auto",
                    analysis_max_cells = Inf,
                    analysis_max_features = Inf,
                    matrix_layers = c("counts", "data"),
                    embed_downloads = c("auto", "always", "never"),
                    embed_max_mb = 50, render = TRUE,
                    overwrite = FALSE, seed = 11L, verbose = TRUE) {
  analyze <- match.arg(analyze)
  filter_raw_barcodes <- match.arg(filter_raw_barcodes)
  filter_low_expression_features <- match.arg(filter_low_expression_features)
  embed_downloads <- match.arg(embed_downloads)
  if (!is.list(scp_args) || !is.list(marker_args)) {
    .sc_stop("scp_args and marker_args must be named lists.")
  }
  valid_named_list <- function(x) !length(x) || (!is.null(names(x)) && all(nzchar(names(x))))
  if (!valid_named_list(scp_args) || !valid_named_list(marker_args)) {
    .sc_stop("scp_args and marker_args must use non-empty argument names.")
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
    .sc_stop("seed must be one finite number.")
  }
  if (!is.character(matrix_layers) || !length(matrix_layers) || anyNA(matrix_layers) ||
      any(!nzchar(trimws(matrix_layers)))) {
    .sc_stop("matrix_layers must be a non-empty character vector such as c('counts', 'data').")
  }
  if (!is.numeric(embed_max_mb) || length(embed_max_mb) != 1L || is.na(embed_max_mb) ||
      !is.finite(embed_max_mb) || embed_max_mb <= 0) {
    .sc_stop("embed_max_mb must be one finite positive number.")
  }
  thresholds <- list(min_features = min_features, min_counts = min_counts,
                     max_features = max_features, max_counts = max_counts)
  valid_threshold <- vapply(thresholds, function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && x >= 0
  }, logical(1))
  if (!all(valid_threshold)) {
    .sc_stop("Cell filtering thresholds must be non-negative numbers.")
  }
  if (min_features > max_features || min_counts > max_counts) {
    .sc_stop("Minimum cell-filtering thresholds cannot exceed their maximum thresholds.")
  }
  for (limit in list(analysis_max_cells = analysis_max_cells,
                     analysis_max_features = analysis_max_features)) {
    if (!is.numeric(limit) || length(limit) != 1L || is.na(limit) || limit < 10L) {
      .sc_stop("analysis_max_cells and analysis_max_features must be at least 10 or Inf.")
    }
  }
  if (!is.character(integration_method) || length(integration_method) != 1L ||
      is.na(integration_method) || !nzchar(integration_method)) {
    .sc_stop("integration_method must be one non-empty character value.")
  }
  if (length(species) != 1L || is.na(species) || !nzchar(trimws(species))) {
    .sc_stop("species must be 'auto' or one non-empty species name (for example human, mouse, rat, or zebrafish).")
  }
  species <- tolower(trimws(species))
  input <- normalizePath(input, mustWork = TRUE)
  output <- .ensure_output(output, overwrite = overwrite)
  .sc_message(verbose, "Reading %s...", input)
  object <- .read_single_cell_rds(input)
  object <- .as_seurat(object)
  object <- .join_split_layers(object, verbose = verbose)
  raw_matrix_object <- object
  raw_meta <- .seurat_metadata(raw_matrix_object)
  primary_annotation <- .primary_annotation_column(raw_meta, annotation_col = annotation_col)
  original_annotation_columns <- unique(c(primary_annotation, .annotation_columns(raw_meta)))
  original_annotation_columns <- original_annotation_columns[
    !is.na(original_annotation_columns) & nzchar(original_annotation_columns)
  ]
  annotation_status <- list(
    source = if (length(original_annotation_columns)) "rds_metadata" else "none",
    preserved = TRUE,
    primary_column = primary_annotation,
    columns = original_annotation_columns,
    original_columns = original_annotation_columns,
    analysis_columns = character(),
    message = if (length(original_annotation_columns)) {
      "Existing annotation columns were captured from the original RDS before filtering and preserved without modification."
    } else {
      "No existing cell-type annotation column was detected in the original RDS; no annotation was invented."
    }
  )
  detected_species <- .detect_species(rownames(raw_matrix_object))

  before <- .object_inventory(object)
  filtering <- .prefilter_raw_barcodes(
    object,
    mode = filter_raw_barcodes,
    min_features = min_features,
    min_counts = min_counts,
    max_features = max_features,
    max_counts = max_counts,
    verbose = verbose
  )
  object <- filtering$object
  inferred <- infer_sample_design(object, sample_col = sample_col, sample_map = sample_map)
  object <- .attach_design(object, inferred)
  status_before <- .analysis_status(before)
  should_analyze <- identical(analyze, "always") || (identical(analyze, "auto") && status_before != "analyzed")
  feature_filtering <- .prefilter_low_expression_features(
    object,
    mode = filter_low_expression_features,
    min_cells = min_cells_per_feature,
    should_analyze = should_analyze,
    verbose = verbose
  )
  object <- feature_filtering$object
  analysis_subset <- list(applied = FALSE, max_cells = analysis_max_cells,
                          cells_before = ncol(object), cells_after = ncol(object))
  if (is.finite(analysis_max_cells) && analysis_max_cells >= 10L &&
      ncol(object) > as.integer(analysis_max_cells) && should_analyze) {
    set.seed(seed)
    keep <- sample(colnames(object), as.integer(analysis_max_cells))
    object <- object[, keep]
    analysis_subset$applied <- TRUE
    analysis_subset$cells_after <- ncol(object)
    .sc_message(verbose, "Analysis subset: SCP will use %s of %s retained cells.",
                analysis_subset$cells_after, analysis_subset$cells_before)
  }
  feature_subset <- list(applied = FALSE, max_features = analysis_max_features,
                         features_before = nrow(object), features_after = nrow(object))
  if (is.finite(analysis_max_features) && analysis_max_features >= 10L &&
      nrow(object) > as.integer(analysis_max_features) && should_analyze) {
    counts <- .layer_data(object, SeuratObject::DefaultAssay(object), "counts")
    detected <- Matrix::rowSums(counts > 0)
    keep <- order(detected, decreasing = TRUE)[seq_len(as.integer(analysis_max_features))]
    object <- object[keep, ]
    feature_subset$applied <- TRUE
    feature_subset$features_after <- nrow(object)
    .sc_message(verbose, "Feature subset: SCP will use %s of %s retained features.",
                feature_subset$features_after, feature_subset$features_before)
  }
  marker_status <- list(
    requested = isTRUE(run_markers) && should_analyze, ran = FALSE,
    group_by = NULL, engine = NULL, error = NULL, scp_error = NULL
  )
  if (should_analyze) {
    object <- .run_scp(object, integration_method = integration_method, scp_args = scp_args, seed = seed, verbose = verbose)
    if (isTRUE(run_markers)) {
      cluster_prefix <- if (identical(tolower(integration_method), "none")) "scRDSreport" else integration_method
      marker_result <- .run_cluster_markers(
        object, marker_args = marker_args, cluster_prefix = cluster_prefix, verbose = verbose
      )
      object <- marker_result$object
      marker_status <- marker_result$status
    }
  } else {
    .sc_message(verbose, "Preserving existing analysis (status: %s).", status_before)
  }
  meta_after <- .seurat_metadata(object)
  analysis_annotation_columns <- unique(c(
    original_annotation_columns[original_annotation_columns %in% names(meta_after)],
    .annotation_columns(meta_after)
  ))
  analysis_annotation_columns <- analysis_annotation_columns[
    !is.na(analysis_annotation_columns) & nzchar(analysis_annotation_columns)
  ]
  annotation_status$analysis_columns <- analysis_annotation_columns
  selected_species <- if (identical(species, "auto")) detected_species$species else species
  species_status <- list(
    requested = species,
    detected = detected_species$species,
    selected = selected_species,
    confidence = detected_species$confidence,
    basis = detected_species$basis
  )
  after <- .object_inventory(object)

  manifest_files <- .export_object(
    object = object, input = input, output = output,
    sample_design = inferred$design, raw_matrix_object = raw_matrix_object,
    matrix_object = object,
    original_annotation_columns = original_annotation_columns,
    analysis_annotation_columns = analysis_annotation_columns,
    matrix_layers = matrix_layers, verbose = verbose
  )
  manifest <- list(
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    input = input,
    output = output,
    sample_source = inferred$source,
    status_before = status_before,
    status_after = .analysis_status(after),
    scp_was_run = should_analyze,
    marker_analysis = marker_status,
    annotation = annotation_status,
    download_embedding = list(
      mode = embed_downloads,
      max_bytes = as.numeric(embed_max_mb) * 1024^2
    ),
    species = species_status,
    integration_method = integration_method,
    barcode_filtering = filtering$summary,
    feature_filtering = feature_filtering$summary,
    analysis_subset = analysis_subset,
    analysis_feature_subset = feature_subset,
    inventory_before = before,
    inventory_after = after,
    sample_design = inferred$design,
    files = manifest_files,
    session_info = utils::capture.output(utils::sessionInfo())
  )
  manifest_path <- file.path(output, ".report", "manifest.rds")
  saveRDS(manifest, manifest_path)
  jsonlite::write_json(manifest, file.path(output, ".report", "manifest.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")

  report_path <- NULL
  if (isTRUE(render)) report_path <- .render_quarto_report(manifest_path, output, verbose = verbose)
  .sc_message(verbose, "Done: %s", output)
  invisible(list(
    output = output,
    report = report_path,
    manifest = manifest_path,
    status_before = status_before,
    status_after = .analysis_status(after),
    sample_design = inferred$design
  ))
}
