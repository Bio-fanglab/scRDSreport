.embedded_download_indices <- function(manifest, output) {
  files <- as.data.frame(manifest$files %||% data.frame(), stringsAsFactors = FALSE)
  if (!nrow(files) || !"path" %in% names(files)) return(integer())
  mode <- as.character(manifest$download_embedding$mode %||% "never")[[1L]]
  if (identical(mode, "never")) return(integer())
  paths <- file.path(output, as.character(files$path))
  sizes <- if ("bytes" %in% names(files)) suppressWarnings(as.numeric(files$bytes)) else rep(NA_real_, nrow(files))
  missing <- is.na(sizes) & file.exists(paths)
  sizes[missing] <- as.numeric(file.info(paths[missing])$size)
  existing <- which(file.exists(paths) & !is.na(sizes))
  if (!length(existing)) return(integer())
  if (identical(mode, "always")) return(existing)
  budget <- suppressWarnings(as.numeric(manifest$download_embedding$max_bytes %||% 0))[[1L]]
  if (!is.finite(budget) || budget <= 0) return(integer())
  section <- if ("section" %in% names(files)) as.character(files$section) else rep("", nrow(files))
  section_priority <- c(
    paste0("module_", c(
      "qc", "reduction", "cluster", "celltype", "differential", "enrichment",
      "pseudotime", "communication", "cell_cycle", "tf", "cnv", "downloads"
    )),
    "annotation_original", "metadata_original", "annotation_analysis", "metadata_analysis",
    "design", "matrix_preview", "analysis_result", "embedding", "inventory", "features",
    "provenance", "expression_matrix", "rds"
  )
  priority <- match(section, section_priority, nomatch = length(section_priority) + 1L)
  artifact_priority <- if ("embed_priority" %in% names(files)) {
    suppressWarnings(as.numeric(files$embed_priority))
  } else {
    rep(100, nrow(files))
  }
  artifact_priority[is.na(artifact_priority)] <- 100
  candidates <- existing[order(priority[existing], artifact_priority[existing], sizes[existing])]
  selected <- integer()
  used <- 0
  for (index in candidates) {
    if (sizes[[index]] <= budget - used) {
      selected <- c(selected, index)
      used <- used + sizes[[index]]
    }
  }
  selected
}

.write_embedded_download_payload <- function(connection, files, paths, indices,
                                             chunk_bytes = 3L * 1024L * 1024L) {
  for (index in indices) {
    input <- file(paths[[index]], open = "rb")
    tryCatch({
      part <- 0L
      repeat {
        chunk <- readBin(input, what = "raw", n = chunk_bytes)
        if (!length(chunk)) break
        part <- part + 1L
        filename <- htmltools::htmlEscape(basename(as.character(files$path[[index]])), attribute = TRUE)
        tag <- paste0(
          "<script type=\"application/octet-stream\" data-scrds-file=\"scrds-embedded-", index,
          "\" data-part=\"", part, "\" data-filename=\"", filename, "\">",
          jsonlite::base64_enc(chunk), "</script>\n"
        )
        writeBin(charToRaw(enc2utf8(tag)), connection)
      }
    }, finally = close(input))
  }
  invisible(NULL)
}

.raw_fixed_match <- function(haystack, needle) {
  if (!length(needle) || length(haystack) < length(needle)) return(NA_integer_)
  candidates <- which(haystack == needle[[1L]])
  candidates <- candidates[candidates <= length(haystack) - length(needle) + 1L]
  for (candidate in candidates) {
    extent <- candidate + seq_along(needle) - 1L
    if (identical(haystack[extent], needle)) return(candidate)
  }
  NA_integer_
}

.inject_embedded_downloads <- function(report_path, manifest_path, output, verbose = TRUE) {
  manifest <- readRDS(manifest_path)
  indices <- .embedded_download_indices(manifest, output)
  if (!length(indices)) return(report_path)
  files <- as.data.frame(manifest$files, stringsAsFactors = FALSE)
  paths <- file.path(output, as.character(files$path))
  temporary <- tempfile(pattern = "report-with-downloads-", tmpdir = dirname(report_path), fileext = ".html")
  on.exit(unlink(temporary), add = TRUE)
  source <- file(report_path, open = "rb")
  target <- file(temporary, open = "wb")
  on.exit(try(close(source), silent = TRUE), add = TRUE)
  on.exit(try(close(target), silent = TRUE), add = TRUE)
  needle <- charToRaw("</body>")
  carry <- raw()
  injected <- FALSE
  repeat {
    chunk <- readBin(source, what = "raw", n = 1024L * 1024L)
    if (!length(chunk)) break
    buffer <- c(carry, chunk)
    if (injected) {
      writeBin(buffer, target)
      carry <- raw()
      next
    }
    hit <- .raw_fixed_match(buffer, needle)
    if (!is.na(hit)) {
      if (hit > 1L) writeBin(buffer[seq_len(hit - 1L)], target)
      .write_embedded_download_payload(target, files, paths, indices)
      writeBin(buffer[seq.int(hit, length(buffer))], target)
      carry <- raw()
      injected <- TRUE
      next
    }
    keep <- min(length(needle) - 1L, length(buffer))
    flush <- length(buffer) - keep
    if (flush > 0L) writeBin(buffer[seq_len(flush)], target)
    carry <- if (keep > 0L) buffer[seq.int(flush + 1L, length(buffer))] else raw()
  }
  if (!injected) {
    if (length(carry)) writeBin(carry, target)
    .write_embedded_download_payload(target, files, paths, indices)
  }
  close(source)
  close(target)
  if (!file.copy(temporary, report_path, overwrite = TRUE, copy.mode = TRUE)) {
    .sc_stop("Could not inject embedded downloads into %s.", report_path)
  }
  .sc_message(
    verbose, "Embedded %s downloadable files after Quarto rendering (%s).",
    length(indices), format(
      structure(sum(file.info(paths[indices])$size), class = "object_size"), units = "auto"
    )
  )
  report_path
}

.render_quarto_report <- function(manifest_path, output, verbose = TRUE) {
  .require_optional("DT", "render tables with CSV and Excel download buttons")
  .require_optional("knitr", "render the Quarto report")
  template <- system.file("quarto", "report.qmd", package = "scRDSreport")
  if (!nzchar(template)) .sc_stop("The installed report template is missing.")
  work_template <- file.path(output, ".report", "report.qmd")
  file.copy(template, work_template, overwrite = TRUE)
  scss <- system.file("quarto", "report.scss", package = "scRDSreport")
  if (nzchar(scss)) file.copy(scss, file.path(output, ".report", "report.scss"), overwrite = TRUE)
  render_directory <- dirname(work_template)
  rendered_local <- file.path(render_directory, "report.html")
  final_report <- file.path(output, "report.html")
  old_working_directory <- setwd(render_directory)
  on.exit(setwd(old_working_directory), add = TRUE)

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
      input = basename(work_template),
      output_file = "report.html",
      execute_params = list(
        manifest = normalizePath(manifest_path),
        defer_embedded_assets = TRUE
      ),
      quiet = !isTRUE(verbose)
    )
  } else {
    quarto_bin <- Sys.which("quarto")
    if (!nzchar(quarto_bin)) {
      conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
      r_prefix <- normalizePath(file.path(R.home(), "..", ".."), mustWork = FALSE)
      candidates <- unique(c(
        if (nzchar(conda_prefix)) file.path(conda_prefix, "bin", "quarto") else character(),
        file.path(r_prefix, "bin", "quarto")
      ))
      candidates <- candidates[file.exists(candidates)]
      if (length(candidates)) quarto_bin <- candidates[[1L]]
    }
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
    writeLines(c(
      "manifest:", paste0("  ", normalizePath(manifest_path)),
      "defer_embedded_assets: true"
    ), params_path)
    args <- c(
      "render", basename(work_template), "--output", "report.html",
      "--execute-params", basename(params_path)
    )
    status <- system2(quarto_bin, args = args)
    if (!identical(status, 0L)) .sc_stop("Quarto rendering failed with status %s.", status)
  }
  if (!file.exists(rendered_local)) {
    .sc_stop("Quarto reported success but did not create %s.", rendered_local)
  }
  .inject_embedded_downloads(rendered_local, manifest_path, output, verbose = verbose)
  if (!file.copy(rendered_local, final_report, overwrite = TRUE, copy.mode = TRUE)) {
    .sc_stop("Could not copy the rendered standalone report to %s.", final_report)
  }
  final_report
}

.subset_top_detected_features <- function(object, max_features) {
  max_features <- min(as.integer(max_features), nrow(object))
  if (max_features >= nrow(object)) return(object)
  counts <- .layer_data(object, SeuratObject::DefaultAssay(object), "counts")
  detected <- Matrix::rowSums(counts > 0)
  keep <- order(detected, decreasing = TRUE)[seq_len(max_features)]
  # SeuratObject 5 currently restores the assay's original feature order but
  # can retain a reordered request for meta.features. Sort the selected row
  # positions before subsetting so both data and metadata receive one order.
  keep <- sort(keep)
  features <- rownames(counts)[keep]
  # Character identities avoid a separate integer-index feature metadata bug.
  object[features, ]
}

#' Build a downloadable single-cell report from an RDS file
#'
#' @param input Path to a single-cell RDS file.
#' @param output Output directory for the report and data files.
#' @param sample_col Optional metadata column containing sample IDs.
#' @param sample_map Optional named sample-to-group vector or data frame with
#'   `sample_id`, `group`, and optionally `replicate`.
#' @param filter_raw_barcodes One of `"auto"`, `"always"`, or `"never"`.
#'   Auto applies the supplied QC thresholds before SCP for raw/partial
#'   objects and preserves already analyzed objects.
#' @param min_features,min_counts,max_features,max_counts,max_percent_mt
#'   Cell-level thresholds used by raw-barcode prefiltering and the full QC
#'   module. Mitochondrial percentage is calculated after the inexpensive
#'   count/feature prefilter so multi-million-droplet objects stay tractable.
#' @param filter_low_expression_features One of `"auto"`, `"always"`, or
#'   `"never"`. Auto removes features detected in too few cells from large raw
#'   objects before SCP analysis; existing analyzed objects are never changed.
#' @param min_cells_per_feature Minimum cells in which a retained feature must
#'   have a non-zero count.
#' @param analyze One of `"auto"`, `"always"`, or `"never"`. Auto runs SCP
#'   unless both a dimensional reduction and a cluster/annotation column exist;
#'   without an explicit `config`, analyzed inputs switch to report-only export.
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
#'   Existing values are never overwritten. With the default
#'   `annotation_mode = "auto_if_missing"`, a new species-matched SingleR column
#'   is attempted only when the input contains no annotation.
#' @param species `"auto"` or one non-empty species name, such as `"human"`,
#'   `"mouse"`, `"rat"`, `"zebrafish"`, `"pig"`, `"cattle"`, `"chicken"`,
#'   `"dog"`, or `"macaque"`. Auto uses stable feature-ID prefixes before gene
#'   symbols. Species-specific modules run only with a matching registered or
#'   explicitly supplied resource; another species is never substituted.
#' @param config A configuration created by \code{\link{report_config}()}. When omitted,
#'   raw/partial objects use the full profile and analyzed objects use the
#'   report-only profile. Explicit full profiles can add advanced modules to an
#'   existing analysis. Existing annotations always take precedence under the
#'   default `auto_if_missing` policy.
#' @param analysis_max_cells Maximum cells passed to SCP. `Inf` (default) keeps
#'   all retained cells; a finite value creates a deterministic analysis subset
#'   while the original RDS remains available for download.
#' @param analysis_max_features Maximum features passed to SCP. `Inf` (default)
#'   keeps all retained features; finite values keep the most detected features.
#' @param matrix_layers Expression layers to export in Matrix Market format.
#' @param embed_downloads One of `"auto"`, `"always"`, or `"never"`.
#'   Auto embeds downloadable files into the standalone HTML up to the total
#'   `embed_max_mb` budget, prioritizing module results, annotation and compact
#'   tables. Always
#'   forces every exported file into the HTML and can create a very large file.
#' @param embed_max_mb Total megabytes available to automatic HTML download
#'   embedding. This limit is ignored when `embed_downloads = "always"`.
#' @param render Whether to render the Quarto HTML report.
#' @param overwrite Whether managed files in a recognized scRDSreport output
#'   directory may be replaced. A non-empty unrecognized directory is refused
#'   even when `TRUE`; unrelated files are never deleted. The input RDS must not
#'   be stored under a managed subdirectory of the same output.
#' @param seed Random seed used by SCP.
#' @param verbose Print progress messages.
#' @return Invisibly, a list containing paths, analysis status, design, and manifest.
#' @export
running <- function(input, output, sample_col = NULL, sample_map = NULL,
                    filter_raw_barcodes = c("auto", "always", "never"),
                    min_features = 200L, min_counts = 500L,
                    max_features = 7500L, max_counts = Inf,
                    max_percent_mt = 20,
                    filter_low_expression_features = c("auto", "always", "never"),
                    min_cells_per_feature = 3L,
                    analyze = c("auto", "always", "never"),
                    integration_method = "none", scp_args = list(),
                    run_markers = TRUE, marker_args = list(),
                    annotation_col = NULL, species = "auto", config = NULL,
                    analysis_max_cells = Inf,
                    analysis_max_features = Inf,
                    matrix_layers = c("counts", "data"),
                    embed_downloads = c("auto", "always", "never"),
                    embed_max_mb = 50, render = TRUE,
                    overwrite = FALSE, seed = 11L, verbose = TRUE) {
  config_was_missing <- is.null(config)
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
  thresholds <- list(
    min_features = min_features, min_counts = min_counts,
    max_features = max_features, max_counts = max_counts,
    max_percent_mt = max_percent_mt
  )
  valid_threshold <- vapply(thresholds, function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && x >= 0
  }, logical(1))
  if (!all(valid_threshold)) {
    .sc_stop("Cell filtering thresholds must be non-negative numbers.")
  }
  if (min_features > max_features || min_counts > max_counts) {
    .sc_stop("Minimum cell-filtering thresholds cannot exceed their maximum thresholds.")
  }
  if (is.finite(max_percent_mt) && max_percent_mt > 100) {
    .sc_stop("max_percent_mt must be between 0 and 100.")
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
  species <- .normalize_species_name(species)
  if (is.null(config)) {
    config <- report_config(
      profile = if (identical(analyze, "never")) "report_only" else "full"
    )
  }
  if (!is.list(config)) {
    .sc_stop("config must be created by report_config().")
  }
  if (exists(".validate_report_config", mode = "function", inherits = TRUE)) {
    config <- .validate_report_config(config)
  }
  if (is.infinite(analysis_max_cells) &&
      is.finite(config$limits$analysis_max_cells %||% Inf)) {
    analysis_max_cells <- config$limits$analysis_max_cells
  }
  if (is.infinite(analysis_max_features) &&
      is.finite(config$limits$analysis_max_features %||% Inf)) {
    analysis_max_features <- config$limits$analysis_max_features
  }
  if (identical(embed_max_mb, 50) &&
      is.numeric(config$limits$embed_max_mb %||% NULL)) {
    embed_max_mb <- config$limits$embed_max_mb
  }
  input <- normalizePath(input, mustWork = TRUE, winslash = "/")
  output <- .normalize_output_path(output)
  .assert_input_outside_managed_output(input, output)
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
  configured_species <- config$resource_overrides$species %||% NULL
  species_status <- .species_selection_status(
    requested = species,
    detected = detected_species,
    configured = configured_species
  )
  selected_species <- species_status$selected
  resource_status <- .resolve_species_resources(config, selected_species)

  before <- .object_inventory(object)
  status_before <- .analysis_status(before)
  if (config_was_missing && identical(analyze, "auto") && identical(status_before, "analyzed")) {
    config <- report_config(profile = "report_only")
    .sc_message(verbose, "Input already contains reductions and annotations/clusters; exporting existing analysis without recomputation.")
  }
  effective_barcode_filter <- if (
    (identical(config$profile, "report_only") || identical(analyze, "never")) &&
      identical(filter_raw_barcodes, "auto")
  ) {
    "never"
  } else if (identical(filter_raw_barcodes, "auto") && !identical(status_before, "analyzed")) {
    # Raw and partial objects must be quality filtered before SCP so reductions,
    # markers and the final object all use the same cell population.
    "always"
  } else {
    filter_raw_barcodes
  }
  filtering <- .prefilter_raw_barcodes(
    object,
    mode = effective_barcode_filter,
    min_features = min_features,
    min_counts = min_counts,
    max_features = max_features,
    max_counts = max_counts,
    max_percent_mt = max_percent_mt,
    mitochondrial_pattern = resource_status$mitochondrial_pattern %||% "^MT-",
    pattern_ignore_case = resource_status$pattern_ignore_case %||% TRUE,
    orgdb = resource_status$orgdb %||% NULL,
    feature_keytype = resource_status$feature_keytype %||% NULL,
    symbol_column = resource_status$symbol_column %||% "SYMBOL",
    verbose = verbose
  )
  object <- filtering$object
  inferred <- infer_sample_design(object, sample_col = sample_col, sample_map = sample_map)
  object <- .attach_design(object, inferred)
  should_analyze <- !identical(config$profile, "report_only") &&
    (identical(analyze, "always") ||
       (identical(analyze, "auto") && status_before != "analyzed"))
  config$qc <- utils::modifyList(
    list(
      filter = should_analyze && !isTRUE(filtering$summary$applied),
      min_features = min_features,
      min_counts = min_counts,
      max_features = max_features,
      max_counts = max_counts,
      max_percent_mt = max_percent_mt,
      max_plot_cells = config$limits$plot_max_cells %||% 50000L
    ),
    config$qc %||% list(),
    keep.null = TRUE
  )
  config$celltype <- utils::modifyList(
    list(
      mode = config$annotation$mode %||% "preserve",
      annotation_column = primary_annotation,
      manual_markers = config$resource_overrides$manual_markers %||% NULL,
      max_plot_cells = config$limits$plot_max_cells %||% 50000L
    ),
    config$celltype %||% list(),
    keep.null = TRUE
  )
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
    object <- .subset_top_detected_features(object, analysis_max_features)
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
  config$qc <- utils::modifyList(
    list(
      orgdb = resource_status$orgdb %||% NULL,
      feature_keytype = resource_status$feature_keytype %||% NULL,
      symbol_column = resource_status$symbol_column %||% "SYMBOL",
      mitochondrial_pattern = resource_status$mitochondrial_pattern %||% "^MT-",
      ribosomal_pattern = resource_status$ribosomal_pattern %||% "^RP[SL]",
      hemoglobin_pattern = resource_status$hemoglobin_pattern %||% "^HB[ABDEGMQZ]",
      pattern_ignore_case = resource_status$pattern_ignore_case %||% TRUE
    ),
    config$qc %||% list(),
    keep.null = TRUE
  )
  config$celltype <- utils::modifyList(
    list(
      auto_annotation_reference = resource_status$auto_annotation_reference %||% NULL,
      allow_reference_download = TRUE
    ),
    config$celltype %||% list(),
    keep.null = TRUE
  )
  full_analysis <- tryCatch(
    .run_full_analysis(
      object = object,
      output = output,
      design = inferred,
      species_info = species_status,
      config = config,
      seed = as.integer(seed),
      verbose = verbose
    ),
    error = function(e) {
      warning(
        sprintf("The full analysis layer could not be initialized: %s", conditionMessage(e)),
        call. = FALSE
      )
      plan <- .build_analysis_plan(
        config,
        context = list(
          species = selected_species,
          has_annotation = length(original_annotation_columns) > 0L,
          has_clusters = length(.cluster_columns(.seurat_metadata(object))) > 0L,
          n_samples = nrow(inferred$design),
          n_cells = ncol(object)
        ),
        resources = resource_status
      )
      for (id in names(plan$modules)) {
        if (isTRUE(plan$modules[[id]]$requested) && !identical(id, "downloads")) {
          plan$modules[[id]]$status <- "failed"
          plan$modules[[id]]$reason <- "full_analysis_initialization_error"
          plan$modules[[id]]$message <- conditionMessage(e)
          plan$modules[[id]]$error <- conditionMessage(e)
        }
      }
      list(
        object = object,
        modules = plan$modules,
        artifacts = list(),
        warnings = conditionMessage(e),
        species = selected_species,
        schema_version = "2.0"
      )
    }
  )
  object <- full_analysis$object
  artifact_id_map <- stats::setNames(
    vapply(full_analysis$artifacts, function(x) as.character(x$path %||% ""), character(1)),
    vapply(full_analysis$artifacts, function(x) as.character(x$artifact_id %||% x$id %||% ""), character(1))
  )
  for (id in names(full_analysis$modules)) {
    module <- full_analysis$modules[[id]]
    if (is.null(module$reason)) module$reason <- module$reason_code %||% ""
    if (is.null(module$summary)) module$summary <- module$details %||% list()
    artifact_ids <- module$artifact_ids %||% character()
    artifact_paths <- unname(artifact_id_map[artifact_ids])
    module$artifacts <- unique(c(artifact_ids, artifact_paths[!is.na(artifact_paths) & nzchar(artifact_paths)]))
    full_analysis$modules[[id]] <- module
  }
  meta_after <- .seurat_metadata(object)
  sample_design_report <- inferred$design
  if ("n_cells" %in% names(sample_design_report)) {
    names(sample_design_report)[names(sample_design_report) == "n_cells"] <- "n_cells_post_qc"
  }
  analysis_sample_counts <- if (".scRDSreport_sample" %in% names(meta_after)) {
    table(as.character(meta_after$.scRDSreport_sample))
  } else {
    integer()
  }
  sample_design_report$n_cells_analysis <- as.integer(
    analysis_sample_counts[as.character(sample_design_report$sample_id)]
  )
  sample_design_report$n_cells_analysis[is.na(sample_design_report$n_cells_analysis)] <- 0L
  celltype_details <- full_analysis$modules$celltype$details %||% list()
  generated_annotation_column <- celltype_details$annotation_column %||% NULL
  generated_annotation_mode <- celltype_details$mode %||% NULL
  is_generated_annotation <- !is.null(generated_annotation_column) &&
    length(generated_annotation_column) == 1L &&
    generated_annotation_column %in% names(meta_after) &&
    !identical(generated_annotation_mode, "preserve")
  analysis_annotation_columns <- unique(c(
    original_annotation_columns[original_annotation_columns %in% names(meta_after)],
    .annotation_columns(meta_after),
    if (is_generated_annotation) generated_annotation_column else character()
  ))
  analysis_annotation_columns <- analysis_annotation_columns[
    !is.na(analysis_annotation_columns) & nzchar(analysis_annotation_columns)
  ]
  annotation_status$analysis_columns <- analysis_annotation_columns
  annotation_status$generated_columns <- if (is_generated_annotation) {
    generated_annotation_column
  } else {
    character()
  }
  if (is_generated_annotation) {
    annotation_status$source <- if (identical(generated_annotation_mode, "manual")) {
      "manual_mapping"
    } else {
      "species_reference_SingleR"
    }
    annotation_status$primary_column <- generated_annotation_column
  } else if (length(original_annotation_columns)) {
    annotation_status$source <- "rds_metadata"
  } else {
    annotation_status$source <- "none"
  }
  celltype_message <- full_analysis$modules$celltype$message %||% NULL
  if (!is.null(celltype_message) && nzchar(celltype_message)) {
    annotation_status$message <- celltype_message
  }
  after <- .object_inventory(object)

  manifest_files <- .export_object(
    object = object, input = input, output = output,
    sample_design = sample_design_report, raw_matrix_object = raw_matrix_object,
    matrix_object = object,
    original_annotation_columns = original_annotation_columns,
    analysis_annotation_columns = analysis_annotation_columns,
    matrix_layers = matrix_layers, verbose = verbose
  )
  analysis_files <- .artifacts_to_manifest_rows(full_analysis$artifacts, output)
  manifest_files <- .bind_rows_fill(manifest_files, analysis_files)
  resource_dependency_names <- unlist(lapply(
    list(resource_status$orgdb, resource_status$txdb),
    function(resource) {
      if (!is.character(resource)) return(character())
      resource[!is.na(resource) & nzchar(resource)]
    }
  ), use.names = FALSE)
  dependency_names <- c(
    "scRDSreport", "SCP", "Seurat", "SeuratObject", "Matrix", "DT",
    "AnnotationDbi", "SingleR", "celldex", "edgeR", "clusterProfiler", "GSVA",
    "msigdbr", "babelgene",
    "monocle3", "SeuratWrappers", "CellChat", "infercnv",
    "ComplexHeatmap", "plotly", "quarto", resource_dependency_names
  )
  dependency_names <- unique(dependency_names[!is.na(dependency_names) & nzchar(dependency_names)])
  dependency_table <- data.frame(
    package = dependency_names,
    installed = vapply(dependency_names, requireNamespace, logical(1), quietly = TRUE),
    version = vapply(dependency_names, function(package) {
      if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
      as.character(utils::packageVersion(package))
    }, character(1)),
    stringsAsFactors = FALSE
  )
  package_version <- dependency_table$version[dependency_table$package == "scRDSreport"]
  if (is.na(package_version) || !length(package_version)) package_version <- "0.3.0"
  created_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  input_sha256 <- .sha256_file(input)
  hash_token <- if (length(input_sha256) && !is.na(input_sha256) && nzchar(input_sha256)) {
    substr(input_sha256, 1L, 12L)
  } else {
    "nohash"
  }
  manifest <- list(
    manifest_schema_version = "2.0",
    run_id = paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "-", hash_token),
    package_version = package_version,
    pipeline_version = "2.0",
    input_sha256 = input_sha256,
    author = list(name = "Anbuengsi", email = "an.bunengsi@qq.com"),
    config = .fa_sanitize_config(config),
    species_resources = .fa_sanitize_config(resource_status),
    dependencies = dependency_table,
    modules = full_analysis$modules,
    artifacts = full_analysis$artifacts,
    warnings = unique(full_analysis$warnings %||% character()),
    created_at = created_at,
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
    sample_design = sample_design_report,
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
    sample_design = sample_design_report,
    modules = full_analysis$modules,
    artifacts = full_analysis$artifacts,
    config = config
  ))
}
