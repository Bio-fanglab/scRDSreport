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

.ensure_output <- function(output, overwrite) {
  output <- normalizePath(output, mustWork = FALSE)
  if (.nonempty_dir(output) && !isTRUE(overwrite)) {
    .sc_stop("Output directory is not empty: %s. Use overwrite = TRUE to replace managed files.", output)
  }
  if (.nonempty_dir(output) && isTRUE(overwrite)) {
    managed <- file.path(output, c("downloads", "tables", "matrices", "figures", ".report", "report.html"))
    managed <- managed[file.exists(managed) | dir.exists(managed)]
    if (length(managed)) unlink(managed, recursive = TRUE, force = TRUE)
  }
  dirs <- file.path(output, c("downloads", "tables", "matrices", "figures", ".report"))
  for (path in c(output, dirs)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
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

.manifest_row <- function(section, label, path, root, rows = NA_integer_, columns = NA_integer_, note = "") {
  data.frame(
    section = section,
    label = label,
    path = .relative_path(path, root),
    format = sub("^.*\\.", "", path),
    rows = as.integer(rows),
    columns = as.integer(columns),
    bytes = if (file.exists(path)) as.numeric(file.info(path)$size) else NA_real_,
    note = note,
    stringsAsFactors = FALSE
  )
}

.slot_or_null <- function(x, slot) {
  if (isS4(x) && slot %in% methods::slotNames(x)) methods::slot(x, slot) else NULL
}
