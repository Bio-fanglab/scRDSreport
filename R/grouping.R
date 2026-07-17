.sample_column_score <- function(name, values, n_cells) {
  values <- as.character(values)
  n_unique <- length(unique(values[!is.na(values) & nzchar(values)]))
  if (!n_unique || n_unique > max(100L, ceiling(n_cells * 0.2))) return(-Inf)
  lname <- tolower(name)
  exact <- c(sample_id = 120, sample = 115, orig.ident = 110, sample_name = 105,
             library_id = 95, library = 90, donor = 85, patient = 85,
             subject = 85, individual = 80, batch = 60)
  score <- if (lname %in% names(exact)) unname(exact[[lname]]) else 0
  if (grepl("sample|library|donor|patient|subject|individual", lname)) score <- score + 50
  if (grepl("cluster|cell.?type|annotation|percent|count|feature", lname)) score <- score - 100
  score + min(n_unique, 20L) - (n_unique == n_cells) * 100
}

.barcode_sample <- function(cell_names) {
  if (is.null(cell_names) || !length(cell_names)) return(NULL)
  prefix <- sub("[_:|].*$", "", cell_names)
  counts <- table(prefix)
  if (length(counts) >= 2L && all(counts >= 10L)) prefix else NULL
}

.choose_sample_ids <- function(object, sample_col = NULL) {
  if (is.character(object) && is.null(dim(object))) {
    return(list(cell_sample = object, source = "provided_names"))
  }
  meta <- .seurat_metadata(object)
  if (!is.null(sample_col)) {
    if (!sample_col %in% names(meta)) .sc_stop("sample_col '%s' is not present in cell metadata.", sample_col)
    return(list(cell_sample = as.character(meta[[sample_col]]), source = sample_col))
  }
  if (ncol(meta)) {
    scores <- vapply(names(meta), function(nm) .sample_column_score(nm, meta[[nm]], nrow(meta)), numeric(1))
    if (any(is.finite(scores) & scores > 0)) {
      best <- names(which.max(scores))
      if (length(unique(as.character(meta[[best]]))) == 1L) {
        inferred <- .barcode_sample(rownames(meta) %||% colnames(object))
        if (!is.null(inferred)) return(list(cell_sample = inferred, source = "cell_barcode_prefix"))
      }
      return(list(cell_sample = as.character(meta[[best]]), source = best))
    }
  }
  inferred <- .barcode_sample(rownames(meta) %||% colnames(object))
  if (!is.null(inferred)) return(list(cell_sample = inferred, source = "cell_barcode_prefix"))
  project <- if (.is_seurat(object)) .slot_or_null(object, "project.name") else NULL
  value <- if (!is.null(project) && nzchar(project)) project else "sample1"
  list(cell_sample = rep(value, max(1L, nrow(meta))), source = "single_project")
}

.parse_sample_names <- function(sample_ids) {
  sample_ids <- unique(as.character(sample_ids))
  out <- data.frame(
    sample_id = sample_ids,
    group = sample_ids,
    replicate = NA_character_,
    confidence = "low",
    needs_review = TRUE,
    grouping_rule = "independent_sample",
    stringsAsFactors = FALSE
  )

  explicit <- regexec("^(.*?)(?:[_. -]+)(?:rep(?:licate)?|r)[_. -]*([0-9]+)$", sample_ids, ignore.case = TRUE)
  matches <- regmatches(sample_ids, explicit)
  for (i in seq_along(matches)) {
    if (length(matches[[i]]) == 3L && nzchar(matches[[i]][2L])) {
      out$group[i] <- matches[[i]][2L]
      out$replicate[i] <- matches[[i]][3L]
      out$confidence[i] <- "high"
      out$needs_review[i] <- FALSE
      out$grouping_rule[i] <- "explicit_replicate_suffix"
    }
  }

  unresolved <- which(is.na(out$replicate))
  subject_pattern <- "(?:^|[_. -])((?:patient|pt|donor|subject|mouse|animal)[_. -]*[A-Za-z0-9]+)(?:[_. -]|$)"
  for (i in unresolved) {
    hit <- regexpr(subject_pattern, sample_ids[i], ignore.case = TRUE, perl = TRUE)
    if (hit[1L] > 0L) {
      subject <- regmatches(sample_ids[i], hit)
      group <- gsub(subject_pattern, "_", sample_ids[i], ignore.case = TRUE, perl = TRUE)
      group <- gsub("^[_. -]+|[_. -]+$", "", group)
      if (nzchar(group)) {
        out$group[i] <- group
        out$replicate[i] <- gsub("^[_. -]+|[_. -]+$", "", subject)
        out$confidence[i] <- "medium"
        out$needs_review[i] <- TRUE
        out$grouping_rule[i] <- "subject_token_candidate"
      }
    }
  }

  unresolved <- which(is.na(out$replicate))
  bare <- regexec("^(.*?)[_. -]+([0-9]+)$", sample_ids)
  bare_matches <- regmatches(sample_ids, bare)
  bases <- vapply(bare_matches, function(x) if (length(x) == 3L) x[2L] else NA_character_, character(1))
  nonrep_context <- grepl("(?:^|[_. -])(day|d|hour|hr|h|time|dose|t)$", bases, ignore.case = TRUE)
  base_counts <- table(bases[!is.na(bases) & !nonrep_context])
  for (i in unresolved) {
    if (!is.na(bases[i]) && !nonrep_context[i] && base_counts[[bases[i]]] >= 2L) {
      out$group[i] <- bases[i]
      out$replicate[i] <- bare_matches[[i]][3L]
      out$confidence[i] <- "medium"
      out$needs_review[i] <- TRUE
      out$grouping_rule[i] <- "numeric_suffix_candidate"
    }
  }
  out
}

.apply_sample_map <- function(design, sample_map) {
  if (is.null(sample_map)) return(design)
  if (is.atomic(sample_map) && !is.null(names(sample_map))) {
    sample_map <- data.frame(sample_id = names(sample_map), group = as.character(sample_map), stringsAsFactors = FALSE)
  }
  if (!is.data.frame(sample_map) || !all(c("sample_id", "group") %in% names(sample_map))) {
    .sc_stop("sample_map must be a named vector or a data frame containing sample_id and group.")
  }
  if (anyDuplicated(sample_map$sample_id)) .sc_stop("sample_map contains duplicate sample_id values.")
  if (anyNA(sample_map$group) || any(!nzchar(as.character(sample_map$group)))) {
    .sc_stop("sample_map contains a missing or empty group.")
  }
  idx <- match(design$sample_id, sample_map$sample_id)
  if (anyNA(idx)) .sc_stop("sample_map is missing samples: %s", paste(design$sample_id[is.na(idx)], collapse = ", "))
  design$group <- as.character(sample_map$group[idx])
  if ("replicate" %in% names(sample_map)) design$replicate <- as.character(sample_map$replicate[idx])
  design$confidence <- "user"
  design$needs_review <- FALSE
  design$grouping_rule <- "user_map"
  design
}

#' Infer sample groups and replicates conservatively
#'
#' @param x A character vector of sample names or a supported single-cell object.
#' @param sample_col Optional metadata column containing sample IDs.
#' @param sample_map Optional named sample-to-group vector or data frame with
#'   `sample_id`, `group`, and optionally `replicate` columns.
#' @return A list with a sample-level design table, per-cell sample IDs, and the source column.
#' @export
infer_sample_design <- function(x, sample_col = NULL, sample_map = NULL) {
  chosen <- .choose_sample_ids(x, sample_col = sample_col)
  chosen$cell_sample[is.na(chosen$cell_sample) | !nzchar(chosen$cell_sample)] <- "unknown_sample"
  design <- .parse_sample_names(chosen$cell_sample)
  design <- .apply_sample_map(design, sample_map)
  design$n_cells <- as.integer(table(factor(chosen$cell_sample, levels = design$sample_id)))
  list(design = design, cell_sample = chosen$cell_sample, source = chosen$source)
}
