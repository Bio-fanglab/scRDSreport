.detect_qc_columns <- function(meta) {
  feature_candidates <- grep("^nFeature([_.]|$)", names(meta), value = TRUE, ignore.case = TRUE)
  count_candidates <- grep("^nCount([_.]|$)", names(meta), value = TRUE, ignore.case = TRUE)
  list(
    features = if (length(feature_candidates)) feature_candidates[[1L]] else NULL,
    counts = if (length(count_candidates)) count_candidates[[1L]] else NULL
  )
}

.qc_keep <- function(meta, min_features = 200L, min_counts = 0L,
                     max_features = Inf, max_counts = Inf) {
  columns <- .detect_qc_columns(meta)
  if (is.null(columns$features)) return(NULL)
  features <- as.numeric(meta[[columns$features]])
  counts <- if (!is.null(columns$counts)) as.numeric(meta[[columns$counts]]) else rep(Inf, length(features))
  keep <- is.finite(features) & features >= min_features & features <= max_features
  if (!is.null(columns$counts)) {
    keep <- keep & is.finite(counts) & counts >= min_counts & counts <= max_counts
  }
  list(keep = keep, columns = columns, features = features, counts = counts)
}

.looks_like_raw_droplets <- function(features, min_features, n_cells) {
  low_fraction <- mean(features < min_features, na.rm = TRUE)
  n_cells >= 10000L &&
    (stats::median(features, na.rm = TRUE) < min_features || low_fraction >= 0.5)
}

.mitochondrial_percent <- function(object) {
  assay <- tryCatch(SeuratObject::DefaultAssay(object), error = function(e) NULL)
  if (is.null(assay) || !nzchar(assay)) return(NULL)
  counts <- .layer_data(object, assay, "counts")
  if (is.null(counts) || !nrow(counts) || !ncol(counts)) return(NULL)
  symbols <- rownames(counts)
  feature_meta <- tryCatch(as.data.frame(object[[assay]][[]]), error = function(e) NULL)
  if (!is.null(feature_meta) && nrow(feature_meta)) {
    candidates <- names(feature_meta)[tolower(names(feature_meta)) %in%
                                        c("gene_symbols", "gene_symbol", "symbol", "gene_name")]
    if (length(candidates)) {
      mapped <- as.character(feature_meta[[candidates[[1L]]]])
      names(mapped) <- rownames(feature_meta)
      replacement <- mapped[match(rownames(counts), names(mapped))]
      use <- !is.na(replacement) & nzchar(replacement)
      symbols[use] <- replacement[use]
    }
  }
  mitochondrial <- grepl("^MT-", symbols, ignore.case = TRUE)
  if (!any(mitochondrial)) return(NULL)
  totals <- as.numeric(Matrix::colSums(counts))
  list(
    percent = 100 * as.numeric(Matrix::colSums(counts[mitochondrial, , drop = FALSE])) /
      pmax(totals, 1),
    features = sum(mitochondrial)
  )
}

.prefilter_raw_barcodes <- function(object, mode = c("auto", "always", "never"),
                                    min_features = 200L, min_counts = 0L,
                                    max_features = Inf, max_counts = Inf,
                                    max_percent_mt = Inf,
                                    verbose = TRUE) {
  mode <- match.arg(mode)
  meta <- .seurat_metadata(object)
  base_summary <- list(
    mode = mode,
    applied = FALSE,
    reason = "Filtering disabled or not required.",
    feature_column = NULL,
    count_column = NULL,
    min_features = min_features,
    min_counts = min_counts,
    max_features = max_features,
    max_counts = max_counts,
    max_percent_mt = max_percent_mt,
    mitochondrial_features = 0L,
    cells_removed_mitochondrial = 0L,
    cells_before = nrow(meta),
    cells_after = nrow(meta),
    cells_removed = 0L
  )
  if (identical(mode, "never")) return(list(object = object, summary = base_summary))

  qc <- .qc_keep(meta, min_features, min_counts, max_features, max_counts)
  if (is.null(qc)) {
    if (identical(mode, "always")) {
      .sc_stop("Barcode filtering was requested, but no nFeature_* metadata column was found.")
    }
    base_summary$reason <- "No nFeature_* metadata column was available for safe automatic filtering."
    return(list(object = object, summary = base_summary))
  }

  raw_like <- .looks_like_raw_droplets(qc$features, min_features, nrow(meta))
  should_filter <- identical(mode, "always") || raw_like
  if (!should_filter) {
    base_summary$feature_column <- qc$columns$features
    base_summary$count_column <- qc$columns$counts
    base_summary$reason <- "Object does not look like an unfiltered droplet matrix."
    return(list(object = object, summary = base_summary))
  }

  keep_n <- sum(qc$keep, na.rm = TRUE)
  if (keep_n < 10L) {
    .sc_stop(
      "Barcode filtering would retain only %s cells. Lower min_features/min_counts or use filter_raw_barcodes = 'never'.",
      keep_n
    )
  }
  .sc_message(
    verbose,
    "Filtering raw barcodes: retaining %s of %s cells (nFeature >= %s, nCount >= %s)...",
    format(keep_n, big.mark = ","), format(nrow(meta), big.mark = ","), min_features, min_counts
  )
  object <- object[, qc$keep]
  mitochondrial <- NULL
  if (is.finite(max_percent_mt)) {
    mitochondrial <- .mitochondrial_percent(object)
    if (!is.null(mitochondrial)) {
      mt_keep <- is.finite(mitochondrial$percent) &
        mitochondrial$percent <= max_percent_mt
      mt_keep_n <- sum(mt_keep)
      if (mt_keep_n < 10L) {
        .sc_stop(
          "Mitochondrial filtering would retain only %s cells. Increase max_percent_mt or use filter_raw_barcodes = 'never'.",
          mt_keep_n
        )
      }
      if (any(!mt_keep)) {
        .sc_message(
          verbose,
          "Filtering mitochondrial outliers: retaining %s of %s cells (percent.mt <= %s)...",
          format(mt_keep_n, big.mark = ","), format(ncol(object), big.mark = ","),
          max_percent_mt
        )
        object <- object[, mt_keep]
      }
      base_summary$mitochondrial_features <- mitochondrial$features
      base_summary$cells_removed_mitochondrial <- keep_n - mt_keep_n
      keep_n <- mt_keep_n
    }
  }
  base_summary$applied <- TRUE
  base_summary$reason <- if (identical(mode, "always")) {
    "User-requested barcode filtering."
  } else {
    "Large raw-droplet signature detected automatically."
  }
  if (is.finite(max_percent_mt)) {
    base_summary$reason <- paste0(
      base_summary$reason,
      if (is.null(mitochondrial)) {
        " No mitochondrial gene symbols matched, so the mitochondrial cutoff was recorded but not applied."
      } else {
        " The mitochondrial cutoff was evaluated after the inexpensive count/feature prefilter."
      }
    )
  }
  base_summary$feature_column <- qc$columns$features
  base_summary$count_column <- qc$columns$counts
  base_summary$cells_after <- ncol(object)
  base_summary$cells_removed <- nrow(meta) - ncol(object)
  list(object = object, summary = base_summary)
}

.prefilter_low_expression_features <- function(
    object, mode = c("auto", "always", "never"), min_cells = 3L,
    should_analyze = TRUE, verbose = TRUE) {
  mode <- match.arg(mode)
  min_cells <- as.integer(min_cells)
  if (length(min_cells) != 1L || is.na(min_cells) || min_cells < 1L) {
    .sc_stop("min_cells_per_feature must be a positive integer.")
  }
  base_summary <- list(
    mode = mode,
    applied = FALSE,
    reason = "Filtering disabled or not required.",
    assay = NULL,
    layer = "counts",
    min_cells = min_cells,
    features_before = nrow(object),
    features_after = nrow(object),
    features_removed = 0L
  )
  if (identical(mode, "never") || !isTRUE(should_analyze)) {
    if (!isTRUE(should_analyze)) {
      base_summary$reason <- "Existing analysis was preserved; features were not changed."
    }
    return(list(object = object, summary = base_summary))
  }

  assay <- tryCatch(SeuratObject::DefaultAssay(object), error = function(e) NULL)
  if (is.null(assay) || !nzchar(assay)) {
    base_summary$reason <- "No default assay was available for feature filtering."
    return(list(object = object, summary = base_summary))
  }
  counts <- .layer_data(object, assay = assay, layer = "counts")
  if (is.null(counts) || !nrow(counts) || !ncol(counts)) {
    base_summary$assay <- assay
    base_summary$reason <- "No counts layer was available for feature filtering."
    return(list(object = object, summary = base_summary))
  }

  detected <- Matrix::rowSums(counts != 0)
  keep <- is.finite(detected) & detected >= min_cells
  should_filter <- identical(mode, "always") ||
    (identical(mode, "auto") && nrow(counts) >= 25000L && any(!keep))
  base_summary$assay <- assay
  if (!should_filter) {
    base_summary$reason <- "Feature matrix is not large enough to require automatic filtering."
    return(list(object = object, summary = base_summary))
  }
  keep_n <- sum(keep)
  if (keep_n < 10L) {
    .sc_stop(
      "Feature filtering would retain only %s features. Lower min_cells_per_feature or use filter_low_expression_features = 'never'.",
      keep_n
    )
  }
  .sc_message(
    verbose,
    "Filtering low-expression features: retaining %s of %s features (detected in >= %s cells)...",
    format(keep_n, big.mark = ","), format(nrow(counts), big.mark = ","), min_cells
  )
  object <- object[rownames(counts)[keep], ]
  base_summary$applied <- TRUE
  base_summary$reason <- if (identical(mode, "always")) {
    "User-requested feature filtering."
  } else {
    "Large feature matrix detected before analysis."
  }
  base_summary$features_after <- keep_n
  base_summary$features_removed <- nrow(counts) - keep_n
  list(object = object, summary = base_summary)
}
