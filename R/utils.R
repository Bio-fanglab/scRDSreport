# The default Quarto template evaluates these packages at report-render time.
# Importing one stable symbol from each package makes the runtime dependencies
# explicit to R package tooling as well as to users installing scRDSreport.
#' @importFrom DT datatable
#' @importFrom ggplot2 ggplot
#' @importFrom ggsci pal_aaas
#' @importFrom htmltools tags
#' @importFrom knitr opts_chunk
NULL

.sc_stop <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

.sc_message <- function(verbose, ...) {
  if (isTRUE(verbose)) message(sprintf(...))
}

.require_optional <- function(package, purpose) {
  if (!requireNamespace(package, quietly = TRUE)) {
    .sc_stop("Package '%s' is required to %s. Please install it first.", package, purpose)
  }
}

.safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^[_ .-]+|[_ .-]+$", "", x)
  ifelse(nzchar(x), x, "unnamed")
}

.is_seurat <- function(x) inherits(x, "Seurat")

.is_sce <- function(x) inherits(x, "SingleCellExperiment")

.nonempty_dir <- function(path) {
  dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)) > 0L
}

.output_owner_marker <- ".scRDSreport-output"

.managed_output_names <- c(
  "analysis", "downloads", "tables", "matrices", "figures",
  ".report", "report.html"
)

.normalize_output_path <- function(output) {
  if (!is.character(output) || length(output) != 1L || is.na(output) ||
      !nzchar(trimws(output))) {
    .sc_stop("output must be one non-empty directory path.")
  }
  raw_output <- gsub("\\\\", "/", trimws(output))
  components <- strsplit(raw_output, "/", fixed = TRUE)[[1L]]
  if (".." %in% components) {
    .sc_stop(
      "output must not contain a '..' parent-directory component. Use a direct dedicated result path."
    )
  }
  output <- normalizePath(output, mustWork = FALSE, winslash = "/")
  forbidden <- unique(c(
    normalizePath(path.expand("~"), mustWork = FALSE, winslash = "/"),
    if (.Platform$OS.type == "windows") character() else "/"
  ))
  is_filesystem_root <- identical(output, "/") ||
    grepl("^[A-Za-z]:/?$", output)
  if (is_filesystem_root || output %in% forbidden) {
    .sc_stop(
      "Refusing broad output directory '%s'. Choose a dedicated scRDSreport result directory.",
      output
    )
  }
  output
}

.is_scRDSreport_output <- function(output) {
  marker <- file.path(output, .output_owner_marker)
  if (file.exists(marker)) {
    first_line <- tryCatch(
      readLines(marker, n = 1L, warn = FALSE, encoding = "UTF-8"),
      error = function(e) character()
    )
    return(length(first_line) == 1L && identical(first_line, "scRDSreport-output-v1"))
  }

  # Version 0.1.x did not write an owner marker. Recognize a legacy output only
  # when its JSON manifest declares the same absolute output path. JSON is used
  # here instead of unserializing an untrusted RDS file.
  legacy_manifest <- file.path(output, ".report", "manifest.json")
  if (!file.exists(legacy_manifest) || !requireNamespace("jsonlite", quietly = TRUE)) {
    return(FALSE)
  }
  manifest <- tryCatch(
    jsonlite::fromJSON(legacy_manifest, simplifyVector = FALSE),
    error = function(e) NULL
  )
  legacy_fields <- c(
    "created_at", "input", "output", "status_before", "status_after", "files"
  )
  recognizable <- is.list(manifest) && !is.null(manifest$output) &&
    (!is.null(manifest$manifest_schema_version) || all(legacy_fields %in% names(manifest)))
  if (!recognizable) {
    return(FALSE)
  }
  recorded_output <- tryCatch(
    normalizePath(as.character(manifest$output)[[1L]], mustWork = FALSE, winslash = "/"),
    error = function(e) ""
  )
  identical(recorded_output, output)
}

.write_output_owner_marker <- function(output) {
  marker <- file.path(output, .output_owner_marker)
  writeLines(
    c(
      "scRDSreport-output-v1",
      "This directory is managed by scRDSreport; unrelated root-level files are preserved."
    ),
    marker,
    useBytes = TRUE
  )
  if (!file.exists(marker)) {
    .sc_stop("Could not create the scRDSreport output owner marker in %s.", output)
  }
  invisible(marker)
}

.assert_input_outside_managed_output <- function(input, output) {
  input <- normalizePath(input, mustWork = TRUE, winslash = "/")
  output <- .normalize_output_path(output)
  managed <- normalizePath(
    file.path(output, .managed_output_names),
    mustWork = FALSE,
    winslash = "/"
  )
  comparable_input <- if (.Platform$OS.type == "windows") tolower(input) else input
  comparable_managed <- if (.Platform$OS.type == "windows") tolower(managed) else managed
  inside <- vapply(comparable_managed, function(path) {
    identical(comparable_input, path) || startsWith(comparable_input, paste0(path, "/"))
  }, logical(1))
  if (any(inside)) {
    .sc_stop(
      paste0(
        "Input RDS is inside a managed scRDSreport output path: %s. ",
        "Move or copy the input outside '%s' before reusing this output directory."
      ),
      input,
      output
    )
  }
  invisible(TRUE)
}

.ensure_output <- function(output, overwrite) {
  output <- .normalize_output_path(output)
  if (.nonempty_dir(output) && !isTRUE(overwrite)) {
    .sc_stop("Output directory is not empty: %s. Use overwrite = TRUE to replace managed files.", output)
  }
  if (.nonempty_dir(output) && isTRUE(overwrite)) {
    if (!.is_scRDSreport_output(output)) {
      .sc_stop(
        paste0(
          "Output directory is non-empty but is not recognized as an scRDSreport output: %s. ",
          "For safety, choose a new or empty directory; unrelated directories are never cleaned."
        ),
        output
      )
    }
    managed <- file.path(output, .managed_output_names)
    managed <- managed[file.exists(managed) | dir.exists(managed)]
    if (length(managed)) unlink(managed, recursive = TRUE, force = TRUE)
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output)) .sc_stop("Could not create output directory: %s.", output)
  .write_output_owner_marker(output)
  dirs <- file.path(output, c(
    "analysis", "downloads", "tables", "matrices", "figures", ".report"
  ))
  for (path in dirs) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  output
}

.write_csv_gz <- function(x, path, row.names = FALSE) {
  con <- gzfile(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  utils::write.csv(x, con, row.names = row.names, na = "")
  invisible(path)
}

.gzip_file <- function(source, target, chunk_size = 1024L * 1024L) {
  input <- file(source, open = "rb")
  output <- gzfile(target, open = "wb")
  on.exit({
    close(input)
    close(output)
  }, add = TRUE)
  repeat {
    bytes <- readBin(input, what = "raw", n = chunk_size)
    if (!length(bytes)) break
    writeBin(bytes, output)
  }
  invisible(target)
}

.relative_path <- function(path, root) {
  path <- normalizePath(path, mustWork = FALSE, winslash = "/")
  root <- paste0(normalizePath(root, mustWork = FALSE, winslash = "/"), "/")
  sub(paste0("^", gsub("([][{}()+*^$|\\?.])", "\\\\\\1", root)), "", path)
}

.manifest_row <- function(section, label, path, root, rows = NA_integer_,
                          columns = NA_integer_, note = "", module = NA_character_,
                          type = "file", row_unit = NA_character_,
                          column_unit = NA_character_, column_dictionary = NA_character_,
                          complete = TRUE) {
  sha256 <- NA_character_
  if (file.exists(path) && requireNamespace("digest", quietly = TRUE)) {
    sha256 <- tryCatch(
      digest::digest(file = path, algo = "sha256", serialize = FALSE),
      error = function(e) NA_character_
    )
  }
  data.frame(
    artifact_id = paste0(
      ifelse(is.na(module) || !nzchar(module), .safe_name(section), .safe_name(module)),
      "__", .safe_name(label)
    ),
    module = module,
    section = section,
    type = type,
    label = label,
    description = note,
    path = .relative_path(path, root),
    format = sub("^.*\\.", "", path),
    rows = as.integer(rows),
    columns = as.integer(columns),
    row_unit = row_unit,
    column_unit = column_unit,
    column_dictionary = column_dictionary,
    bytes = if (file.exists(path)) as.numeric(file.info(path)$size) else NA_real_,
    sha256 = sha256,
    complete = isTRUE(complete),
    note = note,
    stringsAsFactors = FALSE
  )
}

.bind_rows_fill <- function(...) {
  values <- list(...)
  if (length(values) == 1L && is.list(values[[1L]]) &&
      !is.data.frame(values[[1L]])) values <- values[[1L]]
  values <- Filter(function(x) is.data.frame(x) && nrow(x), values)
  if (!length(values)) return(data.frame())
  columns <- unique(unlist(lapply(values, names), use.names = FALSE))
  values <- lapply(values, function(x) {
    missing <- setdiff(columns, names(x))
    for (name in missing) x[[name]] <- NA
    x[columns]
  })
  output <- do.call(rbind, values)
  rownames(output) <- NULL
  output
}

.artifacts_to_manifest_rows <- function(artifacts, output) {
  if (is.null(artifacts) || !length(artifacts)) return(data.frame())
  if (inherits(artifacts, "scRDSreport_artifact_registry")) {
    artifacts <- artifacts$artifacts
  }
  if (inherits(artifacts, "scRDSreport_artifact")) artifacts <- list(artifacts)
  artifacts <- Filter(function(x) is.list(x) && !is.null(x$module) && !is.null(x$path), artifacts)
  if (!length(artifacts)) return(data.frame())

  rows <- lapply(artifacts, function(artifact) {
    path <- as.character(artifact$path %||% "")
    if (!nzchar(path)) return(NULL)
    absolute <- grepl("^(/|[A-Za-z]:[/\\\\])", path)
    full_path <- if (absolute) path else file.path(output, path)
    dictionary <- artifact$column_dictionary %||% list()
    dictionary_text <- if (!length(dictionary)) {
      ""
    } else if (requireNamespace("jsonlite", quietly = TRUE)) {
      jsonlite::toJSON(dictionary, auto_unbox = TRUE, null = "null")
    } else {
      paste(names(dictionary), unlist(dictionary, use.names = FALSE), sep = ": ", collapse = "; ")
    }
    data.frame(
      artifact_id = as.character(artifact$id %||% artifact$artifact_id %||%
                                   paste0(artifact$module, "__", .safe_name(artifact$label %||% basename(path)))),
      module = as.character(artifact$module),
      section = paste0("module_", as.character(artifact$module)),
      type = as.character(artifact$type %||% "file"),
      label = as.character(artifact$label %||% basename(path)),
      description = as.character(artifact$description %||% ""),
      path = .relative_path(full_path, output),
      format = as.character(artifact$format %||% sub("^.*\\.", "", path)),
      rows = as.integer(artifact$rows %||% NA_integer_),
      columns = as.integer(artifact$columns %||% NA_integer_),
      row_unit = as.character(artifact$row_unit %||% ""),
      column_unit = as.character(artifact$column_unit %||% ""),
      column_dictionary = dictionary_text,
      bytes = as.numeric(artifact$bytes %||%
                           if (file.exists(full_path)) file.info(full_path)$size else NA_real_),
      sha256 = as.character(artifact$sha256 %||% NA_character_),
      embed_priority = as.numeric(artifact$embed_priority %||% 100),
      complete = isTRUE(artifact$complete %||% file.exists(full_path)),
      note = as.character(artifact$description %||% ""),
      stringsAsFactors = FALSE
    )
  })
  .bind_rows_fill(rows)
}

.slot_or_null <- function(x, slot) {
  if (isS4(x) && slot %in% methods::slotNames(x)) methods::slot(x, slot) else NULL
}
