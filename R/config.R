.report_module_ids <- function() {
  c(
    "qc",
    "reduction",
    "cluster",
    "celltype",
    "differential",
    "enrichment",
    "pseudotime",
    "communication",
    "cell_cycle",
    "tf",
    "cnv",
    "downloads"
  )
}

.infrastructure_module_ids <- function() {
  "downloads"
}

.profile_module_ids <- function(profile) {
  switch(
    profile,
    full = .report_module_ids(),
    core = c(
      "qc", "reduction", "cluster", "celltype", .infrastructure_module_ids()
    ),
    report_only = .infrastructure_module_ids(),
    character()
  )
}

.normalize_module_selection <- function(modules, profile) {
  known <- .report_module_ids()
  selected <- .profile_module_ids(profile)

  if (is.null(modules)) {
    return(stats::setNames(known %in% selected, known))
  }

  if (is.character(modules)) {
    modules <- unique(modules)
    if ("all" %in% modules) modules <- known
    unknown <- setdiff(modules, known)
    if (length(unknown)) {
      stop("Unknown report module(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    selected <- modules
  } else if (is.logical(modules) && !is.null(names(modules))) {
    if (anyNA(modules) || any(!nzchar(names(modules)))) {
      stop("Named logical module selections cannot contain NA or empty names.", call. = FALSE)
    }
    unknown <- setdiff(names(modules), known)
    if (length(unknown)) {
      stop("Unknown report module(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    requested <- stats::setNames(known %in% selected, known)
    requested[names(modules)] <- modules
    selected <- names(requested)[requested]
  } else if (is.list(modules) && all(c("include", "exclude") %in% names(modules))) {
    include <- unique(as.character(modules$include))
    exclude <- unique(as.character(modules$exclude))
    unknown <- setdiff(c(include, exclude), known)
    if (length(unknown)) {
      stop("Unknown report module(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    selected <- union(selected, include)
    selected <- setdiff(selected, exclude)
  } else {
    stop(
      "modules must be NULL, a character vector, a named logical vector, or a list with include/exclude.",
      call. = FALSE
    )
  }

  # A running() call always needs the downloads/provenance module.
  selected <- union(selected, .infrastructure_module_ids())
  stats::setNames(known %in% selected, known)
}

.normalize_trajectory_root <- function(root) {
  if (is.null(root)) return(NULL)
  if ((is.character(root) || is.numeric(root)) && length(root) == 1L && !is.na(root)) {
    value <- trimws(as.character(root))
    if (!nzchar(value)) stop("trajectory_root cannot be empty.", call. = FALSE)
    return(list(type = "principal_node", value = value))
  }
  if (!is.list(root) || !length(root) || is.null(names(root)) || any(!nzchar(names(root)))) {
    stop(
      "trajectory_root must be NULL, one principal-node value, or a named list describing cells, markers, or a metadata group.",
      call. = FALSE
    )
  }
  root
}

.normalize_cnv_reference <- function(reference) {
  if (is.null(reference)) return(NULL)
  if (!is.character(reference) || anyNA(reference)) {
    stop("cnv_reference must be NULL or a character vector of reference groups.", call. = FALSE)
  }
  reference <- unique(trimws(reference))
  reference <- reference[nzchar(reference)]
  if (!length(reference)) stop("cnv_reference cannot be empty.", call. = FALSE)
  reference
}

.default_report_limits <- function() {
  list(
    analysis_max_cells = Inf,
    analysis_max_features = Inf,
    plot_max_cells = 100000,
    marker_max_cells_per_ident = 2000,
    min_cells_per_group = 20,
    workers = 1L,
    embed_max_mb = 50
  )
}

.normalize_report_limits <- function(limits) {
  if (is.null(limits)) limits <- list()
  if (!is.list(limits) || (length(limits) && (is.null(names(limits)) || any(!nzchar(names(limits)))))) {
    stop("limits must be a named list.", call. = FALSE)
  }
  defaults <- .default_report_limits()
  unknown <- setdiff(names(limits), names(defaults))
  if (length(unknown)) {
    stop("Unknown limit(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  output <- utils::modifyList(defaults, limits)
  positive_names <- setdiff(names(defaults), "workers")
  for (name in positive_names) {
    value <- output[[name]]
    if (!is.numeric(value) || length(value) != 1L || is.na(value) || value <= 0) {
      stop(name, " must be one positive number or Inf.", call. = FALSE)
    }
  }
  workers <- output$workers
  if (!is.numeric(workers) || length(workers) != 1L || is.na(workers) ||
      !is.finite(workers) || workers < 1 || workers != as.integer(workers)) {
    stop("workers must be one positive integer.", call. = FALSE)
  }
  output$workers <- as.integer(workers)
  output
}

#' Configure the complete single-cell report pipeline
#'
#' `report_config()` creates the optional configuration object consumed by the
#' one-command report pipeline. The default `full` profile requests every
#' analysis represented in the original report. Modules whose biological or
#' species-specific prerequisites are unavailable remain visible in the plan
#' with an explicit skipped or needs-input reason; they are never silently
#' removed.
#'
#' @param profile One of `"full"`, `"core"`, or `"report_only"`.
#' @param modules Optional module selection. Supply module IDs, a named logical
#'   vector, or a list with `include` and `exclude`. The downloads module is
#'   always retained.
#' @param annotation_mode One of `"auto_if_missing"`, `"preserve"`, `"auto"`,
#'   or `"manual"`. The default preserves an annotation already stored in the
#'   RDS and otherwise attempts a species-matched reference annotation. It
#'   never overwrites an existing metadata column and never substitutes a
#'   reference from another species.
#' @param differential Differential-analysis strategy. `"auto"` selects a
#'   replicate-aware pseudobulk method when possible and otherwise exports a
#'   bounded, explicitly descriptive effect-size fallback without P values or
#'   FDR. Other accepted values are
#'   `"pseudobulk"`, `"wilcox"`, and `"none"`.
#' @param trajectory_root Optional Monocle root definition: one principal-node
#'   value or a named list describing root cells, markers, or a metadata group.
#' @param cnv_reference Optional character vector of inferCNV reference groups.
#'   References are never guessed.
#' @param resource_overrides Named list of explicit species-resource overrides.
#' @param module_options Named list keyed by canonical module ID. Each value is
#'   a named list passed only to that module, for example
#'   `list(qc = list(filter = TRUE), pseudotime = list(max_cells = 5000))`.
#'   Before CNV inference with a built-in coordinate resource, declare the
#'   input build with `cnv = list(object_genome_assembly = "GRCh38/hg38")`
#'   or explicitly confirm it after verification.
#' @param limits Named resource-limit list. Supported fields are
#'   `analysis_max_cells`, `analysis_max_features`, `plot_max_cells`,
#'   `marker_max_cells_per_ident`, `min_cells_per_group`, `workers`, and
#'   `embed_max_mb`.
#' @return An object of class `scRDSreport_config`.
#' @export
report_config <- function(
    profile = c("full", "core", "report_only"),
    modules = NULL,
    annotation_mode = c("auto_if_missing", "preserve", "auto", "manual"),
    differential = c("auto", "pseudobulk", "wilcox", "none"),
    trajectory_root = NULL,
    cnv_reference = NULL,
    resource_overrides = list(),
    module_options = list(),
    limits = list()) {
  profile <- match.arg(profile)
  annotation_mode <- match.arg(annotation_mode)
  differential <- match.arg(differential)
  if (!is.list(resource_overrides) ||
      (length(resource_overrides) && (is.null(names(resource_overrides)) || any(!nzchar(names(resource_overrides)))))) {
    stop("resource_overrides must be a named list.", call. = FALSE)
  }
  if (!is.list(module_options) ||
      (length(module_options) &&
         (is.null(names(module_options)) || any(!nzchar(names(module_options))) ||
            any(!vapply(module_options, is.list, logical(1)))))) {
    stop("module_options must be a named list whose values are named module-option lists.", call. = FALSE)
  }
  unknown_options <- setdiff(names(module_options), .report_module_ids())
  if (length(unknown_options)) {
    stop("Unknown module option target(s): ", paste(unknown_options, collapse = ", "), call. = FALSE)
  }

  output <- list(
    schema_version = "1.0",
    profile = profile,
    modules = .normalize_module_selection(modules, profile),
    annotation = list(mode = annotation_mode),
    differential = list(
      strategy = differential,
      fallback_grouping = c("group", "annotation", "cluster"),
      max_contrasts = 6L,
      run_group_markers = TRUE,
      max_marker_groups = 20L
    ),
    enrichment = list(descriptive_rank_summary = TRUE),
    trajectory = list(
      root = .normalize_trajectory_root(trajectory_root),
      export_geometry_without_root = TRUE
    ),
    communication = list(export_group_overview = TRUE),
    cnv = list(
      reference_groups = .normalize_cnv_reference(cnv_reference),
      export_readiness = TRUE
    ),
    resource_overrides = resource_overrides,
    module_options = module_options,
    limits = .normalize_report_limits(limits)
  )
  for (id in names(module_options)) {
    existing <- output[[id]]
    if (!is.list(existing)) existing <- list()
    output[[id]] <- utils::modifyList(existing, module_options[[id]], keep.null = TRUE)
  }
  structure(output, class = c("scRDSreport_config", "list"))
}

.validate_report_config <- function(config) {
  if (!inherits(config, "scRDSreport_config")) {
    stop("config must be created by report_config().", call. = FALSE)
  }
  known <- .report_module_ids()
  if (!identical(names(config$modules), known) ||
      !is.logical(config$modules) || anyNA(config$modules)) {
    stop("config contains an invalid module selection.", call. = FALSE)
  }
  invisible(config)
}
