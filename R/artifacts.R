.artifact_schema_fields <- function() {
  c(
    "id", "module", "type", "format", "path", "label", "description",
    "rows", "columns", "row_unit", "column_unit", "column_dictionary",
    "bytes", "sha256", "preview", "embed_priority", "complete"
  )
}

.artifact_token <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^[_ .-]+|[_ .-]+$", "", x)
  if (!nzchar(x)) "unnamed" else x
}

.artifact_format <- function(path) {
  name <- basename(path)
  compound <- regmatches(
    tolower(name),
    regexpr("(csv|tsv|mtx|txt|tar)[.]gz$", tolower(name), perl = TRUE)
  )
  if (length(compound) && nzchar(compound)) return(compound)
  extension <- tools::file_ext(name)
  if (nzchar(extension)) tolower(extension) else "file"
}

.sha256_file <- function(path) {
  if (!is.character(path) || length(path) != 1L || !file.exists(path) || dir.exists(path)) {
    return(NA_character_)
  }
  tools_namespace <- asNamespace("tools")
  if (exists("sha256sum", envir = tools_namespace, inherits = FALSE)) {
    value <- tryCatch(
      unname(get("sha256sum", envir = tools_namespace, inherits = FALSE)(path)[[1L]]),
      error = function(e) NA_character_
    )
    if (is.character(value) && length(value) == 1L && grepl("^[0-9a-fA-F]{64}$", value)) {
      return(tolower(value))
    }
  }
  sha256_binary <- Sys.which("sha256sum")
  if (nzchar(sha256_binary)) {
    output <- tryCatch(
      system2(sha256_binary, args = shQuote(path), stdout = TRUE, stderr = FALSE),
      error = function(e) character()
    )
    value <- if (length(output)) sub("[[:space:]].*$", "", output[[1L]]) else NA_character_
    if (is.character(value) && length(value) == 1L && grepl("^[0-9a-fA-F]{64}$", value)) {
      return(tolower(value))
    }
  }
  shasum_binary <- Sys.which("shasum")
  if (nzchar(shasum_binary)) {
    output <- tryCatch(
      system2(shasum_binary, args = c("-a", "256", shQuote(path)), stdout = TRUE, stderr = FALSE),
      error = function(e) character()
    )
    value <- if (length(output)) sub("[[:space:]].*$", "", output[[1L]]) else NA_character_
    if (is.character(value) && length(value) == 1L && grepl("^[0-9a-fA-F]{64}$", value)) {
      return(tolower(value))
    }
  }
  NA_character_
}

.normalize_column_dictionary <- function(x) {
  if (is.null(x)) return(list())
  if (is.data.frame(x)) {
    if (!all(c("name", "description") %in% names(x))) {
      stop("A data-frame column_dictionary needs name and description columns.", call. = FALSE)
    }
    return(stats::setNames(as.list(as.character(x$description)), as.character(x$name)))
  }
  if (is.character(x) && !is.null(names(x))) x <- as.list(x)
  if (!is.list(x) || (length(x) && (is.null(names(x)) || any(!nzchar(names(x)))))) {
    stop("column_dictionary must be a named list, named character vector, or name/description data frame.", call. = FALSE)
  }
  x
}

.artifact_default_id <- function(module, type, path, label) {
  paste(
    .artifact_token(module),
    .artifact_token(type),
    .artifact_token(if (nzchar(path)) basename(path) else label),
    sep = ":"
  )
}

.new_artifact <- function(
    module, type, path, label,
    description = "", rows = NA_real_, columns = NA_real_,
    row_unit = "", column_unit = "", column_dictionary = list(),
    format = NULL, id = NULL, bytes = NULL, sha256 = NULL,
    preview = NULL, embed_priority = 100L, complete = NULL) {
  scalar_text <- list(module = module, type = type, path = path, label = label, description = description)
  valid_text <- vapply(scalar_text, function(value) {
    is.character(value) && length(value) == 1L && !is.na(value)
  }, logical(1))
  if (!all(valid_text) || !nzchar(module) || !nzchar(type) || !nzchar(label)) {
    stop("module, type, path, label, and description must be scalar text; module, type, and label cannot be empty.", call. = FALSE)
  }
  if (!is.numeric(rows) || length(rows) != 1L || (!is.na(rows) && rows < 0) ||
      !is.numeric(columns) || length(columns) != 1L || (!is.na(columns) && columns < 0)) {
    stop("rows and columns must be non-negative scalar numbers or NA.", call. = FALSE)
  }
  if (!is.numeric(embed_priority) || length(embed_priority) != 1L || is.na(embed_priority)) {
    stop("embed_priority must be one number.", call. = FALSE)
  }
  if (is.null(format)) format <- .artifact_format(path)
  if (!is.character(format) || length(format) != 1L || is.na(format) || !nzchar(format)) {
    stop("format must be one non-empty character value.", call. = FALSE)
  }
  if (is.null(id)) id <- .artifact_default_id(module, type, path, label)
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
    stop("id must be one non-empty character value.", call. = FALSE)
  }
  if (is.null(complete)) complete <- file.exists(path) && !dir.exists(path)
  if (!is.logical(complete) || length(complete) != 1L || is.na(complete)) {
    stop("complete must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(bytes)) {
    bytes <- if (file.exists(path) && !dir.exists(path)) as.numeric(file.info(path)$size) else NA_real_
  }
  if (!is.numeric(bytes) || length(bytes) != 1L || (!is.na(bytes) && bytes < 0)) {
    stop("bytes must be a non-negative scalar number or NA.", call. = FALSE)
  }
  if (is.null(sha256)) sha256 <- if (isTRUE(complete)) .sha256_file(path) else NA_character_
  if (!is.character(sha256) || length(sha256) != 1L) {
    stop("sha256 must be one character value or NA.", call. = FALSE)
  }

  artifact <- list(
    id = id,
    module = module,
    type = type,
    format = format,
    path = path,
    label = label,
    description = description,
    rows = as.numeric(rows),
    columns = as.numeric(columns),
    row_unit = as.character(row_unit),
    column_unit = as.character(column_unit),
    column_dictionary = .normalize_column_dictionary(column_dictionary),
    bytes = as.numeric(bytes),
    sha256 = sha256,
    preview = preview,
    embed_priority = as.numeric(embed_priority),
    complete = complete
  )
  structure(artifact, class = c("scRDSreport_artifact", "list"))
}

.validate_artifact <- function(artifact) {
  if (!inherits(artifact, "scRDSreport_artifact")) {
    stop("artifact must be created by .new_artifact().", call. = FALSE)
  }
  if (!identical(names(artifact), .artifact_schema_fields())) {
    stop("artifact does not match the required schema.", call. = FALSE)
  }
  invisible(artifact)
}

.new_artifact_registry <- function() {
  structure(
    list(schema_version = "1.0", artifacts = list()),
    class = c("scRDSreport_artifact_registry", "list")
  )
}

.register_artifact <- function(registry, artifact, replace = FALSE) {
  if (!inherits(registry, "scRDSreport_artifact_registry")) {
    stop("registry must be created by .new_artifact_registry().", call. = FALSE)
  }
  .validate_artifact(artifact)
  ids <- vapply(registry$artifacts, function(x) x$id, character(1))
  existing <- match(artifact$id, ids)
  if (!is.na(existing) && !isTRUE(replace)) {
    stop("Artifact ID is already registered: ", artifact$id, call. = FALSE)
  }
  if (is.na(existing)) {
    registry$artifacts[[length(registry$artifacts) + 1L]] <- artifact
  } else {
    registry$artifacts[[existing]] <- artifact
  }
  registry
}

.artifact_ids <- function(x) {
  if (inherits(x, "scRDSreport_artifact_registry")) x <- x$artifacts
  if (inherits(x, "scRDSreport_artifact")) x <- list(x)
  if (!is.list(x) || !length(x)) return(character())
  unique(vapply(x, function(artifact) {
    if (inherits(artifact, "scRDSreport_artifact")) artifact$id else NA_character_
  }, character(1), USE.NAMES = FALSE)[!vapply(x, function(artifact) {
    !inherits(artifact, "scRDSreport_artifact")
  }, logical(1))])
}

.artifact_manifest <- function(registry) {
  if (!inherits(registry, "scRDSreport_artifact_registry")) {
    stop("registry must be created by .new_artifact_registry().", call. = FALSE)
  }
  stats::setNames(unclass(registry$artifacts), .artifact_ids(registry))
}

.artifact_table <- function(registry) {
  if (!inherits(registry, "scRDSreport_artifact_registry")) {
    stop("registry must be created by .new_artifact_registry().", call. = FALSE)
  }
  if (!length(registry$artifacts)) {
    return(data.frame(
      id = character(), module = character(), type = character(),
      format = character(), path = character(), label = character(),
      description = character(), rows = numeric(), columns = numeric(),
      row_unit = character(), column_unit = character(),
      column_dictionary = I(vector("list", 0L)), bytes = numeric(),
      sha256 = character(), preview = I(vector("list", 0L)),
      embed_priority = numeric(), complete = logical(),
      stringsAsFactors = FALSE
    ))
  }
  output <- do.call(rbind, lapply(registry$artifacts, function(x) {
    data.frame(
      id = x$id, module = x$module, type = x$type, format = x$format,
      path = x$path, label = x$label, description = x$description,
      rows = x$rows, columns = x$columns, row_unit = x$row_unit,
      column_unit = x$column_unit, bytes = x$bytes, sha256 = x$sha256,
      embed_priority = x$embed_priority, complete = x$complete,
      stringsAsFactors = FALSE
    )
  }))
  output$column_dictionary <- I(lapply(registry$artifacts, `[[`, "column_dictionary"))
  output$preview <- I(lapply(registry$artifacts, `[[`, "preview"))
  output <- output[.artifact_schema_fields()]
  rownames(output) <- NULL
  output
}

.module_result <- function(value = NULL, artifacts = list()) {
  if (inherits(artifacts, "scRDSreport_artifact")) artifacts <- list(artifacts)
  if (!inherits(artifacts, "scRDSreport_artifact_registry") && !is.list(artifacts)) {
    stop("artifacts must be an artifact, artifact list, or artifact registry.", call. = FALSE)
  }
  structure(list(value = value, artifacts = artifacts), class = c("scRDSreport_module_result", "list"))
}
