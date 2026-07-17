.attach_design <- function(object, inferred) {
  if (!.is_seurat(object)) return(object)
  map <- inferred$design
  idx <- match(inferred$cell_sample, map$sample_id)
  confidence <- if ("confidence" %in% names(map)) map$confidence else rep("unknown", nrow(map))
  needs_review <- if ("needs_review" %in% names(map)) map$needs_review else rep(TRUE, nrow(map))
  grouping_rule <- if ("grouping_rule" %in% names(map)) map$grouping_rule else rep("unspecified", nrow(map))
  meta <- data.frame(
    .scRDSreport_sample = inferred$cell_sample,
    .scRDSreport_group = map$group[idx],
    .scRDSreport_replicate = map$replicate[idx],
    .scRDSreport_design_confidence = confidence[idx],
    .scRDSreport_design_needs_review = needs_review[idx],
    .scRDSreport_grouping_rule = grouping_rule[idx],
    row.names = colnames(object),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  object <- SeuratObject::AddMetaData(object, metadata = meta)
  object
}

.analysis_status <- function(inventory) {
  inventory$value[inventory$item == "analysis_status"][[1L]]
}

.base_binding_call <- function(function_name, ...) {
  function_object <- get(function_name, envir = baseenv(), inherits = FALSE)
  function_object(...)
}

.replace_scp_import <- function(name, value) {
  imports <- parent.env(asNamespace("SCP"))
  if (!exists(name, envir = imports, inherits = FALSE)) return(NULL)
  locked <- .base_binding_call("bindingIsLocked", name, imports)
  original <- get(name, envir = imports, inherits = FALSE)
  if (locked) .base_binding_call("unlockBinding", name, imports)
  assign(name, value, envir = imports)
  if (locked) .base_binding_call("lockBinding", name, imports)
  list(name = name, imports = imports, value = original, locked = locked)
}

.restore_scp_imports <- function(states) {
  for (state in rev(states)) {
    if (is.null(state)) next
    if (.base_binding_call("bindingIsLocked", state$name, state$imports)) {
      .base_binding_call("unlockBinding", state$name, state$imports)
    }
    assign(state$name, state$value, envir = state$imports)
    if (isTRUE(state$locked)) .base_binding_call("lockBinding", state$name, state$imports)
  }
  invisible(NULL)
}

.apply_scp_seurat5_compat <- function(verbose = TRUE) {
  if (utils::packageVersion("SeuratObject") < "5.0.0") return(list())

  get_assay_data_compat <- function(object, slot = NULL, layer = NULL, ...) {
    selected <- layer %||% slot %||% "data"
    SeuratObject::GetAssayData(object = object, layer = selected, ...)
  }
  set_assay_data_compat <- function(object, slot = NULL, new.data, assay = NULL, layer = NULL, ...) {
    selected <- layer %||% slot %||% "data"
    SeuratObject::SetAssayData(
      object = object, layer = selected, new.data = new.data, assay = assay, ...
    )
  }
  states <- Filter(Negate(is.null), list(
    .replace_scp_import("GetAssayData", get_assay_data_compat),
    .replace_scp_import("SetAssayData", set_assay_data_compat)
  ))
  if (length(states)) {
    .sc_message(
      verbose,
      "Applied the SCP 0.5.x compatibility adapter for SeuratObject >= 5 (slot -> layer)."
    )
  }
  states
}

.repair_cluster_identity <- function(object, prefix, verbose = TRUE) {
  meta <- .seurat_metadata(object)
  candidates <- setdiff(.cluster_columns(meta), "active_ident")
  preferred <- candidates[grepl(
    paste0("^", .safe_name(prefix), ".*_SNN_res[.]"), candidates, ignore.case = TRUE
  )]
  source <- if (length(preferred)) preferred[[1L]] else if (length(candidates)) candidates[[1L]] else NULL
  if (is.null(source)) {
    warning("SCP completed but no valid multi-level cluster column was found.", call. = FALSE)
    return(object)
  }
  identities <- factor(meta[[source]])
  SeuratObject::Idents(object) <- identities
  canonical <- data.frame(
    .scRDSreport_cluster = identities,
    row.names = colnames(object),
    check.names = FALSE
  )
  object <- SeuratObject::AddMetaData(object, canonical)
  .sc_message(verbose, "Using %s as the canonical cluster identity.", source)
  object
}

.run_scp <- function(object, integration_method = "none", scp_args = list(), seed = 11L, verbose = TRUE) {
  .require_optional("SCP", "analyze a minimally processed single-cell object")
  .require_optional("SeuratObject", "work with a Seurat object")
  compatibility_state <- .apply_scp_seurat5_compat(verbose = verbose)
  if (length(compatibility_state)) {
    on.exit(.restore_scp_imports(compatibility_state), add = TRUE)
  }
  if (!.is_seurat(object)) .sc_stop("SCP analysis requires a Seurat object after conversion.")

  cells <- ncol(object)
  features <- nrow(object)
  if (cells < 10L || features < 10L) .sc_stop("SCP analysis needs at least 10 cells and 10 features.")
  large_object <- cells > 20000L
  default_nhvf <- min(if (large_object) 1000L else 2000L, features)
  default_dims <- min(if (large_object) 30L else 50L, cells - 1L, features - 1L)
  if (large_object) {
    .sc_message(
      verbose,
      "Large-object defaults: %s HVF and %s linear-reduction dimensions (override with scp_args).",
      default_nhvf, default_dims
    )
  }
  defaults <- list(
    srt = object,
    prefix = "scRDSreport",
    do_normalization = !.has_normalized_data(object),
    nHVF = default_nhvf,
    linear_reduction_dims = default_dims,
    neighbor_k = min(20L, cells - 1L),
    force_linear_reduction = FALSE,
    force_nonlinear_reduction = FALSE,
    seed = as.integer(seed)
  )

  if (!identical(tolower(integration_method), "none")) {
    n_samples <- length(unique(object$.scRDSreport_sample))
    if (n_samples < 2L) .sc_stop("Integration requested but only one sample was detected.")
    .sc_message(verbose, "Running SCP integration (%s) using inferred sample IDs as batch...", integration_method)
    args <- utils::modifyList(list(
      srtMerge = object,
      batch = ".scRDSreport_sample",
      integration_method = integration_method,
      do_normalization = !.has_normalized_data(object),
      nHVF = default_nhvf,
      linear_reduction_dims = default_dims,
      neighbor_k = min(20L, cells - 1L),
      seed = as.integer(seed)
    ), scp_args)
    result <- do.call(SCP::Integration_SCP, args)
    return(.repair_cluster_identity(result, prefix = integration_method, verbose = verbose))
  }

  .sc_message(verbose, "Running SCP::Standard_SCP() to complete the analysis...")
  result <- do.call(SCP::Standard_SCP, utils::modifyList(defaults, scp_args))
  .repair_cluster_identity(result, prefix = "scRDSreport", verbose = verbose)
}

.run_cluster_markers <- function(object, marker_args = list(), cluster_prefix = "scRDSreport", verbose = TRUE) {
  meta <- .seurat_metadata(object)
  candidates <- setdiff(.cluster_columns(meta), "active_ident")
  preferred <- c(
    intersect(".scRDSreport_cluster", candidates),
    candidates[grepl(paste0("^", .safe_name(cluster_prefix), ".*clusters$"), candidates, ignore.case = TRUE)]
  )
  preferred <- unique(preferred)
  group_by <- if (length(preferred)) preferred[[1L]] else if (length(candidates)) candidates[[1L]] else NULL

  if (is.null(group_by)) {
    active <- .slot_or_null(object, "active.ident")
    if (is.null(active) || length(unique(active)) < 2L) {
    return(list(
      object = object,
      status = list(requested = TRUE, ran = FALSE, group_by = NULL, engine = NULL,
                    error = "No multi-level cluster identity was found.", scp_error = NULL)
      ))
    }
    group_by <- ".scRDSreport_cluster"
    cluster_meta <- data.frame(value = as.character(active), row.names = colnames(object), check.names = FALSE)
    names(cluster_meta) <- group_by
    object <- SeuratObject::AddMetaData(object, cluster_meta)
  }

  marker_cell_limit <- if (ncol(object) > 20000L) 500L else 2000L
  .sc_message(
    verbose,
    "Running SCP::RunDEtest() for cluster markers (%s; up to %s cells per cluster)...",
    group_by, marker_args$max.cells.per.ident %||% marker_cell_limit
  )
  args <- utils::modifyList(list(
    srt = object,
    group_by = group_by,
    markers_type = "all",
    only.pos = TRUE,
    max.cells.per.ident = marker_cell_limit,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  ), marker_args)
  group_by <- args$group_by
  error_message <- NULL
  result <- tryCatch(
    do.call(SCP::RunDEtest, args),
    error = function(e) {
      error_message <<- conditionMessage(e)
      object
    }
  )
  if (!is.null(error_message)) {
    .sc_message(verbose, "SCP::RunDEtest() failed; trying Seurat::FindAllMarkers() compatibility fallback...")
    fallback_error <- NULL
    result <- tryCatch({
      .require_optional("Seurat", "run the Seurat 5 cluster-marker fallback")
      meta <- .seurat_metadata(object)
      if (!group_by %in% names(meta)) .sc_stop("Marker group column is missing: %s", group_by)
      SeuratObject::Idents(object) <- factor(meta[[group_by]])
      fallback_args <- list(
        object = object,
        features = marker_args$features %||% NULL,
        assay = marker_args$assay %||% NULL,
        logfc.threshold = log(marker_args$fc.threshold %||% 1.5, base = marker_args$base %||% 2),
        test.use = marker_args$test.use %||% "wilcox",
        min.pct = marker_args$min.pct %||% 0.1,
        min.diff.pct = marker_args$min.diff.pct %||% -Inf,
        only.pos = marker_args$only.pos %||% TRUE,
        max.cells.per.ident = marker_args$max.cells.per.ident %||% marker_cell_limit,
        random.seed = marker_args$seed %||% 11L,
        latent.vars = marker_args$latent.vars %||% NULL,
        min.cells.feature = marker_args$min.cells.feature %||% 3,
        min.cells.group = marker_args$min.cells.group %||% 3,
        verbose = verbose
      )
      markers <- do.call(Seurat::FindAllMarkers, fallback_args)
      object@misc$scRDSreport_cluster_markers <- markers
      object
    }, error = function(e) {
      fallback_error <<- conditionMessage(e)
      object
    })
    if (!is.null(fallback_error)) {
      warning(sprintf(
        "Cluster marker analysis failed in both SCP and Seurat fallback: %s",
        fallback_error
      ), call. = FALSE)
    }
    return(list(
      object = result,
      status = list(
        requested = TRUE,
        ran = is.null(fallback_error),
        group_by = group_by,
        engine = if (is.null(fallback_error)) "Seurat::FindAllMarkers fallback" else NULL,
        error = fallback_error,
        scp_error = error_message
      )
    ))
  }
  list(
    object = result,
    status = list(requested = TRUE, ran = TRUE, group_by = group_by,
                  engine = "SCP::RunDEtest", error = NULL, scp_error = NULL)
  )
}

.detect_species <- function(features) {
  features <- sub("[.][0-9]+$", "", as.character(features))
  prefixes <- c(
    mouse = "^ENSMUSG", rat = "^ENSRNOG", zebrafish = "^ENSDARG",
    pig = "^ENSSSCG", cattle = "^ENSBTAG", chicken = "^ENSGALG",
    macaque = "^ENSMMUG", dog = "^ENSCAFG"
  )
  fractions <- vapply(prefixes, function(pattern) mean(grepl(pattern, features, ignore.case = TRUE)), numeric(1))
  if (length(fractions) && max(fractions) >= 0.2) {
    detected <- names(which.max(fractions))
    return(list(
      species = detected, confidence = "high",
      basis = paste0(sub("^\\^", "", prefixes[[detected]]), " feature IDs")
    ))
  }
  human_ensembl <- mean(grepl("^ENSG[0-9]", features, ignore.case = TRUE))
  if (human_ensembl >= 0.2) return(list(species = "human", confidence = "high", basis = "ENSG feature IDs"))
  if (mean(grepl("^FBgn[0-9]", features, ignore.case = TRUE)) >= 0.2) {
    return(list(species = "drosophila", confidence = "high", basis = "FBgn feature IDs"))
  }
  if (mean(grepl("^WBGene[0-9]", features, ignore.case = TRUE)) >= 0.2) {
    return(list(species = "c_elegans", confidence = "high", basis = "WBGene feature IDs"))
  }

  symbols <- features[grepl("^[A-Za-z][A-Za-z0-9._-]*$", features)]
  if (!length(symbols)) return(list(species = "unknown", confidence = "low", basis = "unrecognized feature IDs"))
  human_like <- mean(symbols == toupper(symbols) & grepl("[A-Z]", symbols))
  title_case_like <- mean(grepl("^[A-Z][a-z]", symbols))
  if (human_like >= 0.8 && human_like > title_case_like * 2) {
    return(list(species = "human", confidence = "medium", basis = "gene-symbol capitalization"))
  }
  if (title_case_like >= 0.5 && title_case_like > human_like * 2) {
    return(list(species = "unknown", confidence = "low", basis = "title-case symbols cannot distinguish mouse, rat, and other species"))
  }
  list(species = "unknown", confidence = "low", basis = "ambiguous gene symbols")
}
