.module_catalog <- function() {
  ids <- .report_module_ids()
  labels <- c(
    "Quality control",
    "Dimensionality reduction",
    "Clustering",
    "Cell annotation and composition",
    "Differential analysis",
    "Functional enrichment",
    "Trajectory and pseudotime analysis",
    "Cell-cell communication",
    "Cell-cycle analysis",
    "Transcription-factor expression",
    "Copy-number variation (inferCNV)",
    "Data downloads and reproducibility"
  )
  categories <- c(
    "core", "core", "core", "annotation", "differential", "enrichment",
    "trajectory", "communication", "cell_cycle", "regulation", "cnv", "output"
  )
  dependencies <- list(
    character(), "qc", "reduction", "cluster", "celltype", "differential",
    c("reduction", "celltype"), "celltype", "celltype", "celltype", "celltype",
    character()
  )
  steps <- list(
    c("metrics", "thresholds", "filtering", "summaries", "violin", "scatter"),
    c("hvg", "pca_variance", "elbow", "pca", "umap", "tsne", "three_dimensional"),
    c("neighbors", "clustering", "cluster_reduction_plots"),
    c("annotation", "composition", "sample_composition", "annotation_reduction_plots"),
    c("celltype_one_vs_rest", "sample_or_group_contrasts", "volcano", "upset", "expression_heatmap"),
    c("go", "kegg", "gsea", "gsva", "dotplot", "upset", "enrichment_heatmap"),
    c("trajectory", "root_selection", "ordering", "distribution", "dynamic_genes", "functional_enrichment"),
    c("global_network", "interaction_strength", "ligand_receptor", "pathway", "sender_receiver", "sample_specific"),
    c("scores", "phase", "composition", "reduction_plots"),
    c("feature_mapping", "average_expression", "heatmap", "reduction_plots", "summary"),
    c("inputs", "reference", "infercnv", "signal", "chromosome_summary", "heatmap"),
    c("objects", "matrices", "tables", "figures", "manifest", "session_info")
  )
  data.frame(
    id = ids,
    label = labels,
    category = categories,
    depends_on = I(dependencies),
    steps = I(steps),
    stringsAsFactors = FALSE
  )
}

.module_time <- function(x) {
  if (is.null(x)) return(NULL)
  format(x, "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
}

.module_schema_fields <- function() {
  c(
    "id", "label", "requested", "eligible", "status", "reason", "message",
    "engine", "version", "parameters", "seed", "timing", "warnings", "error",
    "artifact_ids"
  )
}

.new_module_record <- function(
    id, label = id, requested = TRUE, eligible = NA,
    status = "planned", reason = "eligibility_pending", message = "",
    engine = NULL, version = NULL, parameters = list(), seed = NULL,
    started_at = NULL, finished_at = NULL, elapsed_seconds = NA_real_,
    warnings = character(), error = NULL, artifact_ids = character()) {
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
    stop("module id must be one non-empty character value.", call. = FALSE)
  }
  if (!is.logical(requested) || length(requested) != 1L || is.na(requested)) {
    stop("requested must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(eligible) || length(eligible) != 1L) {
    stop("eligible must be TRUE, FALSE, or NA.", call. = FALSE)
  }
  structure(
    list(
      id = id,
      label = as.character(label),
      requested = requested,
      eligible = eligible,
      status = as.character(status),
      reason = as.character(reason),
      message = as.character(message),
      engine = engine,
      version = version,
      parameters = parameters,
      seed = seed,
      timing = list(
        started_at = started_at,
        finished_at = finished_at,
        elapsed_seconds = as.numeric(elapsed_seconds)
      ),
      warnings = as.character(warnings),
      error = error,
      artifact_ids = unique(as.character(artifact_ids))
    ),
    class = c("scRDSreport_module_record", "list")
  )
}

.context_value <- function(context, name, default = NA) {
  if (!is.list(context) || is.null(context[[name]])) default else context[[name]]
}

.known_flag <- function(value) {
  is.logical(value) && length(value) == 1L && !is.na(value)
}

.any_prerequisite_flag <- function(...) {
  values <- list(...)
  known <- vapply(values, .known_flag, logical(1))
  known_values <- unlist(values[known], use.names = FALSE)
  if (length(known_values) && any(known_values)) return(TRUE)
  if (all(known)) return(FALSE)
  NA
}

.resource_field_available <- function(resources, field) {
  if (is.null(resources) || is.null(resources[[field]])) return(FALSE)
  length(resources[[field]]) > 0L && !all(is.na(resources[[field]]))
}

.eligibility_result <- function(eligible = NA, reason = "eligibility_pending", message = "") {
  list(eligible = eligible, reason = reason, message = message)
}

.plan_assembly_tokens <- function(value) {
  value <- unique(trimws(as.character(value %||% character())))
  value <- value[!is.na(value) & nzchar(value)]
  if (!length(value)) return(character())
  pieces <- unique(c(value, unlist(strsplit(value, "[/,;|[:space:]]+"), use.names = FALSE)))
  tokens <- tolower(gsub("[^[:alnum:]]", "", pieces))
  unique(tokens[nzchar(tokens)])
}

.celltype_eligibility <- function(config, context, resources) {
  mode <- config$annotation$mode
  has_annotation <- .context_value(context, "has_annotation")
  celltype_config <- utils::modifyList(
    config$celltype %||% list(),
    config$module_options$celltype %||% list(),
    keep.null = TRUE
  )
  explicit_reference <- !is.null(celltype_config$reference) &&
    length(celltype_config$reference) > 0L
  if (identical(mode, "preserve")) {
    if (.known_flag(has_annotation)) {
      if (isTRUE(has_annotation)) return(.eligibility_result(TRUE, "ready", "Existing annotations will be preserved."))
      return(.eligibility_result(FALSE, "missing_annotation", "No existing annotation is available and preserve mode will not invent one."))
    }
    return(.eligibility_result())
  }
  if (identical(mode, "auto")) {
    if (is.null(resources)) return(.eligibility_result())
    available <- explicit_reference ||
      (.resource_field_available(resources, "auto_annotation_reference") &&
         .resource_field_available(resources, "orgdb"))
    if (available) return(.eligibility_result(TRUE, "ready", "Explicit automatic-annotation mode has compatible resources."))
    return(.eligibility_result(FALSE, "annotation_resources_unavailable", "Automatic annotation has no compatible species resource."))
  }
  if (identical(mode, "auto_if_missing")) {
    if (.known_flag(has_annotation) && isTRUE(has_annotation)) {
      return(.eligibility_result(
        TRUE, "ready",
        "An existing RDS annotation is available and will be preserved."
      ))
    }
    if (is.null(resources)) return(.eligibility_result())
    available <- explicit_reference ||
      (.resource_field_available(resources, "auto_annotation_reference") &&
         .resource_field_available(resources, "orgdb"))
    if (available) {
      return(.eligibility_result(
        TRUE, "ready",
        "No existing annotation is available; a species-matched reference annotation may be generated."
      ))
    }
    return(.eligibility_result(
      FALSE, "annotation_resources_unavailable",
      "No existing annotation or compatible species-matched reference is available."
    ))
  }
  if (is.null(resources)) return(.eligibility_result())
  if (.resource_field_available(resources, "manual_markers")) {
    return(.eligibility_result(TRUE, "ready", "Explicit manual marker definitions were supplied."))
  }
  .eligibility_result(FALSE, "manual_markers_required", "Manual annotation requires explicit marker definitions.")
}

.evaluate_module_eligibility <- function(id, config, context, resources) {
  always_ready <- c("qc", "reduction", "cluster", "downloads")
  if (id %in% always_ready) return(.eligibility_result(TRUE, "ready", "No additional biological prerequisite."))

  has_annotation <- .context_value(context, "has_annotation")
  has_clusters <- .context_value(context, "has_clusters")
  has_group <- .context_value(context, "has_group")
  has_cell_cycle_scores <- .context_value(context, "has_cell_cycle_scores")
  annotation_or_cluster <- .any_prerequisite_flag(has_annotation, has_clusters)
  differential_grouping <- .any_prerequisite_flag(has_group, annotation_or_cluster)
  n_samples <- .context_value(context, "n_samples")
  enough_samples <- if (is.numeric(n_samples) && length(n_samples) == 1L && !is.na(n_samples)) n_samples >= 2L else NA

  if (identical(id, "celltype")) return(.celltype_eligibility(config, context, resources))
  if (identical(id, "differential")) {
    if (identical(config$differential$strategy, "none")) {
      return(.eligibility_result(FALSE, "differential_disabled", "Differential analysis was disabled explicitly."))
    }
    if (isTRUE(differential_grouping) && isTRUE(enough_samples)) {
      return(.eligibility_result(TRUE, "ready", "Annotation/cluster and at least two samples are available."))
    }
    if (isTRUE(differential_grouping)) {
      return(.eligibility_result(
        TRUE,
        "descriptive_fallback_ready",
        paste(
          "Annotation or cluster groups are available. If biological replication is",
          "insufficient, the report will show descriptive group effect sizes without P values."
        )
      ))
    }
    if (isFALSE(differential_grouping)) {
      return(.eligibility_result(FALSE, "annotation_or_cluster_required", "Differential analysis needs an annotation or cluster field."))
    }
    if (isFALSE(enough_samples)) return(.eligibility_result(FALSE, "two_samples_required", "At least two samples are required."))
    return(.eligibility_result())
  }
  if (identical(id, "enrichment")) {
    if (is.null(resources)) return(.eligibility_result())
    has_ora <- .resource_field_available(resources, "orgdb") &&
      .resource_field_available(resources, "kegg_code")
    has_gene_sets <- .resource_field_available(resources, "gene_sets") ||
      .resource_field_available(resources, "gene_sets_strategy")
    if (has_ora && has_gene_sets) {
      return(.eligibility_result(
        TRUE, "ready",
        "Species-matched GO/KEGG and gene-set resources are registered; runtime checks installed packages and resource availability."
      ))
    }
    return(.eligibility_result(FALSE, "enrichment_resources_unavailable", "Complete GO/KEGG/GSEA/GSVA resources are unavailable for this species."))
  }
  if (identical(id, "pseudotime")) {
    n_cells <- .context_value(context, "n_cells")
    if (is.numeric(n_cells) && length(n_cells) == 1L && !is.na(n_cells) && n_cells < 100L) {
      return(.eligibility_result(FALSE, "insufficient_cells", "Trajectory inference requires at least 100 cells."))
    }
    if (is.null(config$trajectory$root)) {
      return(.eligibility_result(
        TRUE,
        "root_candidate_only",
        paste(
          "Trajectory geometry and root candidates may be generated, but directed",
          "pseudotime and dynamic-gene claims require a biologically justified root."
        )
      ))
    }
    return(.eligibility_result(TRUE, "ready", "Trajectory inference and an explicit root are configured."))
  }
  if (identical(id, "communication")) {
    if (.known_flag(has_annotation) && !isTRUE(has_annotation)) {
      if (isTRUE(has_clusters)) {
        return(.eligibility_result(
          TRUE, "cluster_context_only",
          "Cluster-size context can be exported, but communication inference still requires a biological annotation."
        ))
      }
      return(.eligibility_result(FALSE, "annotation_required", "CellChat requires an existing or explicitly requested annotation."))
    }
    if (is.null(resources)) return(.eligibility_result())
    if (!.resource_field_available(resources, "cellchat_db")) {
      return(.eligibility_result(
        TRUE, "communication_context_only",
        "No CellChat database is registered for this species; annotation-group context will still be exported."
      ))
    }
    return(.eligibility_result(
      TRUE, "ready",
      "Annotation and a species-matched CellChat database are registered; runtime checks the installed CellChat package."
    ))
  }
  if (identical(id, "cell_cycle")) {
    if (isTRUE(has_cell_cycle_scores)) {
      return(.eligibility_result(
        TRUE, "existing_scores_ready",
        "Existing S.Score, G2M.Score, and Phase values can be exported without species gene resources."
      ))
    }
    if (is.null(resources)) return(.eligibility_result())
    strategy <- resources$cell_cycle_strategy
    if (!is.null(strategy) && length(strategy) && !identical(strategy, "user_supplied")) {
      return(.eligibility_result(TRUE, "ready", "A species-specific cell-cycle strategy is available."))
    }
    return(.eligibility_result(FALSE, "cell_cycle_resources_unavailable", "Cell-cycle genes must be supplied for this species."))
  }
  if (identical(id, "tf")) {
    if (is.null(resources)) return(.eligibility_result())
    strategy <- resources$tf_catalog_strategy
    if (!is.null(strategy) && length(strategy) && !identical(strategy, "user_supplied")) {
      return(.eligibility_result(TRUE, "ready", "A transcription-factor catalog strategy is available."))
    }
    return(.eligibility_result(FALSE, "tf_catalog_required", "A transcription-factor catalog must be supplied for this species."))
  }
  if (identical(id, "cnv")) {
    cnv_config <- utils::modifyList(
      config$cnv %||% list(), config$module_options$cnv %||% list(), keep.null = TRUE
    )
    reference_groups <- as.character(cnv_config$reference_groups %||% character())
    has_reference_groups <- any(!is.na(reference_groups) & nzchar(trimws(reference_groups)))
    if (!has_reference_groups) {
      return(.eligibility_result(
        FALSE, "cnv_reference_required",
        "inferCNV reference groups are never guessed; runtime can still export a CNV-readiness overview."
      ))
    }
    explicit_coordinates <- any(vapply(
      c("gene_order", "gtf", "txdb"),
      function(field) {
        value <- cnv_config[[field]]
        if (is.null(value) || !length(value)) return(FALSE)
        if (is.character(value)) return(any(!is.na(value) & nzchar(trimws(value))))
        TRUE
      }, logical(1)
    ))
    if (is.null(resources) && !explicit_coordinates) return(.eligibility_result())
    registered <- if (is.null(resources)) NULL else resources$genome_assembly
    requested <- cnv_config$object_genome_assembly %||% cnv_config$requested_genome_assembly %||%
      cnv_config$assembly %||% cnv_config$genome_assembly %||% cnv_config$confirmed_genome_assembly
    confirmation <- cnv_config$genome_assembly_confirmed %||% FALSE
    if (is.character(confirmation) && length(confirmation) && is.null(requested)) requested <- confirmation[[1L]]
    registered_tokens <- .plan_assembly_tokens(registered)
    requested_tokens <- .plan_assembly_tokens(requested)
    assembly_matches <- length(registered_tokens) && length(requested_tokens) &&
      length(intersect(registered_tokens, requested_tokens)) > 0L
    mismatch <- length(registered_tokens) && length(requested_tokens) && !assembly_matches
    assembly_confirmed <- !mismatch && (isTRUE(confirmation) || assembly_matches)
    if (!explicit_coordinates && !assembly_confirmed) {
      return(.eligibility_result(
        FALSE, "cnv_genome_assembly_confirmation_required",
        "A registered TxDb is not sufficient by itself; confirm the input object's genome assembly before inferCNV."
      ))
    }
    coordinates_available <- explicit_coordinates ||
      .resource_field_available(resources, "txdb") || .resource_field_available(resources, "gtf")
    if (coordinates_available) {
      return(.eligibility_result(TRUE, "ready", "CNV reference and explicitly accepted genome coordinates are configured."))
    }
    return(.eligibility_result(FALSE, "cnv_genome_resource_required", "inferCNV needs a compatible TxDb or GTF."))
  }
  .eligibility_result()
}

.build_analysis_plan <- function(config = report_config(), context = list(), resources = NULL) {
  .validate_report_config(config)
  if (!is.list(context)) stop("context must be a named list.", call. = FALSE)
  if (is.null(resources)) {
    species <- .context_value(context, "species", "auto")
    if (!identical(species, "auto")) resources <- .resolve_species_resources(config, species)
  }
  catalog <- .module_catalog()
  records <- lapply(seq_len(nrow(catalog)), function(i) {
    id <- catalog$id[[i]]
    requested <- isTRUE(config$modules[[id]])
    if (!requested) {
      return(.new_module_record(
        id = id, label = catalog$label[[i]], requested = FALSE, eligible = FALSE,
        status = "not_requested", reason = "profile_or_module_selection",
        message = "The module remains in the plan but was not requested by this profile.",
        parameters = list(depends_on = catalog$depends_on[[i]], steps = catalog$steps[[i]])
      ))
    }
    decision <- .evaluate_module_eligibility(id, config, context, resources)
    status <- if (isFALSE(decision$eligible)) "skipped" else "planned"
    .new_module_record(
      id = id, label = catalog$label[[i]], requested = TRUE,
      eligible = decision$eligible, status = status, reason = decision$reason,
      message = decision$message,
      parameters = list(depends_on = catalog$depends_on[[i]], steps = catalog$steps[[i]])
    )
  })
  names(records) <- catalog$id
  structure(
    list(
      schema_version = "1.0",
      profile = config$profile,
      created_at = .module_time(Sys.time()),
      species = if (is.null(resources)) NULL else resources$species,
      modules = records
    ),
    class = c("scRDSreport_analysis_plan", "list")
  )
}

.module_by_id <- function(plan, id) {
  if (!inherits(plan, "scRDSreport_analysis_plan")) {
    stop("plan must be created by .build_analysis_plan().", call. = FALSE)
  }
  if (!id %in% names(plan$modules)) stop("Unknown module in plan: ", id, call. = FALSE)
  plan$modules[[id]]
}

.set_module_eligibility <- function(plan, id, eligible, reason, message = "") {
  record <- .module_by_id(plan, id)
  if (!is.logical(eligible) || length(eligible) != 1L || is.na(eligible)) {
    stop("eligible must be TRUE or FALSE.", call. = FALSE)
  }
  record$eligible <- eligible
  record$reason <- as.character(reason)
  record$message <- as.character(message)
  record$status <- if (!record$requested) "not_requested" else if (eligible) "planned" else "skipped"
  plan$modules[[id]] <- record
  plan
}

.module_engine_version <- function(engine) {
  if (is.null(engine) || !is.character(engine) || length(engine) != 1L || !nzchar(engine)) return(NULL)
  package <- sub("::.*$", "", engine)
  if (identical(package, engine) || !requireNamespace(package, quietly = TRUE)) return(NULL)
  as.character(utils::packageVersion(package))
}

.run_module_safely <- function(
    record, fun, ..., engine = NULL, version = NULL,
    parameters = NULL, seed = NULL) {
  if (!inherits(record, "scRDSreport_module_record")) {
    stop("record must be created by .new_module_record().", call. = FALSE)
  }
  if (!record$requested || isFALSE(record$eligible)) {
    record$status <- if (!record$requested) "not_requested" else "skipped"
    return(structure(
      list(record = record, value = NULL, artifacts = list()),
      class = c("scRDSreport_module_run", "list")
    ))
  }

  started <- Sys.time()
  record$status <- "running"
  record$timing$started_at <- .module_time(started)
  record$engine <- engine
  record$version <- if (is.null(version)) .module_engine_version(engine) else as.character(version)
  if (!is.null(parameters)) record$parameters <- parameters
  if (!is.null(seed)) record$seed <- seed

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (!is.null(record$seed)) set.seed(as.integer(record$seed))

  warning_messages <- character()
  error_message <- NULL
  value <- tryCatch(
    withCallingHandlers(
      {
        if (!is.function(fun)) stop("Module runner requires a function.")
        fun(...)
      },
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
  )

  artifacts <- list()
  if (inherits(value, "scRDSreport_module_result")) {
    artifacts <- value$artifacts
    value <- value$value
  }
  artifact_ids <- .artifact_ids(artifacts)
  finished <- Sys.time()
  record$timing$finished_at <- .module_time(finished)
  record$timing$elapsed_seconds <- as.numeric(difftime(finished, started, units = "secs"))
  record$warnings <- unique(warning_messages)
  record$artifact_ids <- unique(c(record$artifact_ids, artifact_ids))

  if (is.null(error_message)) {
    record$status <- "completed"
    record$reason <- "completed"
    record$message <- if (length(warning_messages)) {
      "Module completed with captured warnings."
    } else {
      "Module completed successfully."
    }
    record$error <- NULL
  } else {
    record$status <- "failed"
    record$reason <- "module_error"
    record$message <- "The optional module failed; remaining modules may continue."
    record$error <- error_message
  }

  structure(
    list(record = record, value = value, artifacts = artifacts),
    class = c("scRDSreport_module_run", "list")
  )
}

.run_plan_module_safely <- function(plan, id, fun, ...) {
  record <- .module_by_id(plan, id)
  run <- .run_module_safely(record, fun, ...)
  plan$modules[[id]] <- run$record
  structure(
    list(plan = plan, record = run$record, value = run$value, artifacts = run$artifacts),
    class = c("scRDSreport_plan_run", "list")
  )
}

.plan_as_data_frame <- function(plan) {
  if (!inherits(plan, "scRDSreport_analysis_plan")) {
    stop("plan must be created by .build_analysis_plan().", call. = FALSE)
  }
  do.call(rbind, lapply(plan$modules, function(record) {
    data.frame(
      id = record$id,
      label = record$label,
      requested = record$requested,
      eligible = record$eligible,
      status = record$status,
      reason = record$reason,
      message = record$message,
      engine = if (is.null(record$engine)) NA_character_ else record$engine,
      version = if (is.null(record$version)) NA_character_ else record$version,
      seed = if (is.null(record$seed)) NA_real_ else as.numeric(record$seed),
      started_at = if (is.null(record$timing$started_at)) NA_character_ else record$timing$started_at,
      finished_at = if (is.null(record$timing$finished_at)) NA_character_ else record$timing$finished_at,
      elapsed_seconds = record$timing$elapsed_seconds,
      warning_count = length(record$warnings),
      error = if (is.null(record$error)) NA_character_ else record$error,
      artifact_ids = paste(record$artifact_ids, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
}
