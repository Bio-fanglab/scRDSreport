.dependency_row <- function(component, kind = "r_package", source = "cran",
                            modules = "package", species = "*",
                            minimum_version = NA_character_, remote = NA_character_,
                            required = TRUE, reason = "") {
  data.frame(
    component = component,
    kind = kind,
    source = source,
    modules = paste(modules, collapse = "|"),
    species = paste(species, collapse = "|"),
    minimum_version = minimum_version,
    install_spec = remote,
    required = isTRUE(required),
    reason = reason,
    stringsAsFactors = FALSE
  )
}

.species_dependency_modules <- c(
  "celltype", "enrichment", "communication", "cell_cycle", "tf", "cnv"
)

.dependency_catalog <- function() {
  rows <- list(
    .dependency_row("R", "system", "system", minimum_version = "4.2.0",
                    reason = "R runtime required by scRDSreport."),
    .dependency_row("BiocParallel", source = "bioconductor",
                    reason = "Parallel backend used by SCP marker analysis."),
    .dependency_row("digest", reason = "SHA-256 provenance for exported files."),
    .dependency_row("DT", reason = "Interactive downloadable report tables."),
    .dependency_row("Matrix", reason = "Sparse expression matrices."),
    .dependency_row("SCP", source = "github", minimum_version = "0.5.6",
                    remote = "zhanghao-njmu/SCP",
                    reason = "Core completion of raw or partial single-cell objects."),
    .dependency_row("Seurat", reason = "Single-cell analysis and conversion."),
    .dependency_row("SeuratObject", reason = "Seurat object and assay access."),
    .dependency_row("ggplot2", minimum_version = "3.4.0",
                    reason = "Static scientific figures."),
    .dependency_row("ggsci", reason = "Scientific colour palettes."),
    .dependency_row("htmltools", reason = "Standalone HTML report components."),
    .dependency_row("jsonlite", reason = "Manifest and embedded-download metadata."),
    .dependency_row("knitr", reason = "Quarto R execution engine."),

    .dependency_row("quarto", modules = "render", required = FALSE,
                    reason = "Recommended R interface to the Quarto CLI; the CLI is checked separately."),
    .dependency_row("Quarto CLI", "system", "system", modules = "render",
                    minimum_version = "1.3.0",
                    reason = "External renderer required when render = TRUE."),

    .dependency_row("AnnotationDbi", source = "bioconductor",
                    modules = c("celltype", "enrichment", "cell_cycle", "tf", "cnv"),
                    reason = "Identifier mapping and organism annotation databases."),
    .dependency_row("celldex", source = "bioconductor", modules = "celltype",
                    species = c("human", "mouse"),
                    reason = "Registered human and mouse SingleR references."),
    .dependency_row("SingleR", source = "bioconductor", modules = "celltype",
                    species = c("human", "mouse"),
                    reason = "Automatic reference annotation when the RDS lacks labels."),
    .dependency_row("SummarizedExperiment", source = "bioconductor",
                    modules = c("celltype", "pseudotime", "cnv"),
                    reason = "Reference, trajectory, CNV, and SCE-compatible containers."),
    .dependency_row("scrapper", source = "bioconductor", modules = "celltype",
                    species = c("human", "mouse"), required = FALSE,
                    reason = "Recommended cluster-level SingleR aggregation backend."),
    .dependency_row("edgeR", source = "bioconductor", modules = "differential",
                    reason = "Replicate-aware pseudobulk inference."),
    .dependency_row("clusterProfiler", source = "bioconductor", modules = "enrichment",
                    reason = "GO, KEGG, ORA, and GSEA analysis."),
    .dependency_row("GSVA", source = "bioconductor", modules = "enrichment",
                    reason = "Sample-level gene-set scores."),
    .dependency_row("msigdbr", modules = "enrichment",
                    reason = "Species-aware MSigDB gene sets."),
    .dependency_row("monocle3", source = "github", modules = "pseudotime",
                    minimum_version = "1.0.0", remote = "cole-trapnell-lab/monocle3",
                    reason = "Trajectory geometry and rooted pseudotime."),
    .dependency_row("SingleCellExperiment", source = "bioconductor",
                    modules = "pseudotime",
                    reason = "Reduced-dimension access for Monocle 3 objects."),
    .dependency_row("igraph", modules = "pseudotime", required = FALSE,
                    reason = "Recommended trajectory graph export."),
    .dependency_row("HDF5 development libraries", "system", "system",
                    modules = "pseudotime", required = TRUE,
                    reason = "Build prerequisite for the BPCells dependency of current Monocle 3."),
    .dependency_row("CellChat", source = "github", modules = "communication",
                    species = c("human", "mouse"), minimum_version = "1.5.0",
                    remote = "jinworks/CellChat",
                    reason = "Cell-cell communication inference with registered human or mouse databases."),
    .dependency_row("babelgene", modules = "cell_cycle", species = "nonhuman",
                    reason = "Human-to-target ortholog mapping for non-human cell-cycle genes."),
    .dependency_row("infercnv", source = "bioconductor", modules = "cnv",
                    reason = "Expression-derived CNV analysis after an explicit normal reference is supplied."),
    .dependency_row("JAGS", "system", "system", modules = "cnv",
                    minimum_version = "4.0.0",
                    reason = "External runtime required by infercnv through rjags."),
    .dependency_row("GenomicFeatures", source = "bioconductor", modules = "cnv",
                    reason = "Gene coordinates from a TxDb."),
    .dependency_row("rtracklayer", source = "bioconductor", modules = "cnv",
                    required = FALSE,
                    reason = "Recommended GTF import path when a TxDb is unavailable."),

    .dependency_row("org.Hs.eg.db", source = "bioconductor",
                    modules = c("celltype", "enrichment", "cell_cycle", "tf", "cnv"),
                    species = "human", reason = "Human gene identifier and GO mapping."),
    .dependency_row("org.Mm.eg.db", source = "bioconductor",
                    modules = c("celltype", "enrichment", "cell_cycle", "tf", "cnv"),
                    species = "mouse", reason = "Mouse gene identifier and GO mapping."),
    .dependency_row("org.Rn.eg.db", source = "bioconductor",
                    modules = c("enrichment", "cell_cycle", "tf", "cnv"),
                    species = "rat", reason = "Rat gene identifier and GO mapping."),
    .dependency_row("org.Dr.eg.db", source = "bioconductor",
                    modules = c("enrichment", "cell_cycle", "tf", "cnv"),
                    species = "zebrafish", reason = "Zebrafish gene identifier and GO mapping."),
    .dependency_row("org.Ss.eg.db", source = "bioconductor",
                    modules = c("enrichment", "cell_cycle", "tf", "cnv"),
                    species = "pig", reason = "Pig gene identifier and GO mapping."),
    .dependency_row("org.Bt.eg.db", source = "bioconductor",
                    modules = c("enrichment", "cell_cycle", "tf", "cnv"),
                    species = "cattle", reason = "Cattle gene identifier and GO mapping."),
    .dependency_row("org.Gg.eg.db", source = "bioconductor",
                    modules = c("enrichment", "cell_cycle", "tf", "cnv"),
                    species = "chicken", reason = "Chicken gene identifier and GO mapping."),
    .dependency_row("org.Cf.eg.db", source = "bioconductor",
                    modules = c("enrichment", "cell_cycle", "tf", "cnv"),
                    species = "dog", reason = "Dog gene identifier and GO mapping."),
    .dependency_row("org.Mmu.eg.db", source = "bioconductor",
                    modules = c("enrichment", "cell_cycle", "tf", "cnv"),
                    species = "macaque", reason = "Macaque gene identifier and GO mapping."),
    .dependency_row("TxDb.Hsapiens.UCSC.hg38.knownGene", source = "bioconductor",
                    modules = "cnv", species = "human",
                    reason = "Registered hg38 gene coordinates; use only after genome-build confirmation."),
    .dependency_row("TxDb.Mmusculus.UCSC.mm10.knownGene", source = "bioconductor",
                    modules = "cnv", species = "mouse",
                    reason = "Registered mm10 gene coordinates; use only after genome-build confirmation.")
  )
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

.dependency_split <- function(x) {
  if (is.na(x) || !nzchar(x)) character() else strsplit(x, "|", fixed = TRUE)[[1L]]
}

.dependency_input_features <- function(input) {
  object <- .read_single_cell_rds(input)
  if (is.list(object) && !.is_seurat(object) && !.is_sce(object)) {
    key <- intersect(names(object), c("counts", "count", "raw_counts", "matrix", "expr", "expression"))
    if (length(key)) object <- object[[key[[1L]]]]
  }
  rownames(object) %||% character()
}

.resolve_dependency_species <- function(species = "auto", input = NULL) {
  if (!is.character(species) || !length(species) || anyNA(species) ||
      any(!nzchar(trimws(species)))) {
    stop("species must contain one or more non-empty species names, 'auto', or 'all'.", call. = FALSE)
  }
  normalized <- unique(vapply(species, .normalize_species_name, character(1)))
  if ("all" %in% normalized) {
    if (length(normalized) != 1L) stop("'all' cannot be combined with other species.", call. = FALSE)
    return(list(
      species = names(.common_species_registry()), requested = "all",
      detected = NULL, needs_species = FALSE,
      message = "All built-in species resource packages were requested."
    ))
  }
  if ("auto" %in% normalized) {
    if (length(normalized) != 1L) stop("'auto' cannot be combined with other species.", call. = FALSE)
    if (is.null(input)) {
      return(list(
        species = character(), requested = "auto", detected = NULL,
        needs_species = TRUE,
        message = "No input RDS was supplied, so species-specific packages were not guessed."
      ))
    }
    features <- .dependency_input_features(input)
    detected <- .detect_species(features)
    selected <- .normalize_species_name(detected$species %||% "unknown")
    if (selected %in% names(.common_species_registry())) {
      return(list(
        species = selected, requested = "auto", detected = detected,
        needs_species = FALSE,
        message = paste0("Species dependency selection used ", detected$basis, ".")
      ))
    }
    return(list(
      species = character(), requested = "auto", detected = detected,
      needs_species = TRUE,
      message = paste0(
        "The input species could not be resolved safely (", detected$basis,
        "). Supply species explicitly before installing organism databases."
      )
    ))
  }
  supported <- names(.common_species_registry())
  unsupported <- setdiff(normalized, supported)
  list(
    species = intersect(normalized, supported), requested = normalized,
    detected = NULL, needs_species = length(unsupported) > 0L,
    message = if (length(unsupported)) {
      paste0(
        "No built-in dependency resources are registered for: ",
        paste(unsupported, collapse = ", "),
        ". Select only generic modules or supply a supported species."
      )
    } else {
      "Explicit species selection controls dependency installation."
    }
  )
}

.dependency_request <- function(profile = c("full", "core", "report_only"),
                                modules = NULL, species = "auto", input = NULL,
                                render = TRUE) {
  profile <- match.arg(profile)
  if (!is.logical(render) || length(render) != 1L || is.na(render)) {
    stop("render must be TRUE or FALSE.", call. = FALSE)
  }
  module_selection <- .normalize_module_selection(modules, profile)
  selected_modules <- names(module_selection)[module_selection]
  species_info <- .resolve_dependency_species(species, input)
  catalog <- .dependency_catalog()
  selected <- logical(nrow(catalog))
  for (index in seq_len(nrow(catalog))) {
    entry_modules <- .dependency_split(catalog$modules[[index]])
    entry_species <- .dependency_split(catalog$species[[index]])
    module_match <- "package" %in% entry_modules ||
      (isTRUE(render) && "render" %in% entry_modules) ||
      any(entry_modules %in% selected_modules)
    species_match <- "*" %in% entry_species ||
      ("nonhuman" %in% entry_species && any(species_info$species != "human")) ||
      any(entry_species %in% species_info$species)
    selected[[index]] <- module_match && species_match
  }
  catalog <- catalog[selected, , drop = FALSE]
  unresolved_species <- isTRUE(species_info$needs_species) &&
    any(selected_modules %in% .species_dependency_modules)
  list(
    profile = profile,
    modules = selected_modules,
    species = species_info,
    catalog = catalog,
    unresolved_species = unresolved_species,
    render = render
  )
}

.dependency_package_state <- function(package, minimum_version = NA_character_) {
  installed <- nzchar(system.file(package = package))
  disk_version <- if (installed) {
    tryCatch(as.character(utils::packageVersion(package)), error = function(e) NA_character_)
  } else {
    NA_character_
  }
  already_loaded <- package %in% loadedNamespaces()
  loaded_version <- if (already_loaded) {
    tryCatch(as.character(getNamespaceVersion(package)), error = function(e) NA_character_)
  } else {
    NA_character_
  }
  version <- if (already_loaded && !is.na(loaded_version)) loaded_version else disk_version
  version_ok <- is.na(minimum_version) || !nzchar(minimum_version) ||
    (!is.na(version) && utils::compareVersion(version, minimum_version) >= 0L)
  load_error <- NULL
  loadable <- if (installed && !version_ok) {
    already_loaded
  } else {
    tryCatch(
      {
        suppressWarnings(loadNamespace(package))
        TRUE
      },
      error = function(e) {
        load_error <<- conditionMessage(e)
        FALSE
      }
    )
  }
  if (loadable && package %in% loadedNamespaces()) {
    version <- tryCatch(
      as.character(getNamespaceVersion(package)),
      error = function(e) version
    )
    version_ok <- is.na(minimum_version) || !nzchar(minimum_version) ||
      (!is.na(version) && utils::compareVersion(version, minimum_version) >= 0L)
  }
  list(
    installed = installed,
    loadable = loadable,
    version = version,
    version_ok = version_ok,
    ready = installed && loadable && version_ok,
    detail = if (!installed) {
      "Package is not installed in the active R library paths."
    } else if (!version_ok) {
      if (already_loaded && !identical(loaded_version, disk_version)) {
        paste0(
          "Namespace version ", loaded_version, " is loaded while installed files are ",
          disk_version, "; restart R before verification. Required version: ",
          minimum_version, "."
        )
      } else {
        paste0("Installed/runtime version ", version, " is older than required ", minimum_version, ".")
      }
    } else if (!loadable) {
      paste0("Package files exist but the namespace cannot load: ", load_error %||% "unknown error")
    } else {
      "Package namespace loaded successfully."
    }
  )
}

.find_quarto_cli <- function() {
  candidates <- character()
  if (requireNamespace("quarto", quietly = TRUE)) {
    candidates <- c(candidates, tryCatch(quarto::quarto_path(), error = function(e) character()))
  }
  configured <- Sys.getenv("QUARTO_PATH", unset = "")
  if (nzchar(configured)) {
    candidates <- c(
      candidates,
      if (dir.exists(configured)) {
        c(
          file.path(configured, c("quarto", "quarto.exe", "quarto.cmd")),
          file.path(configured, "bin", c("quarto", "quarto.exe", "quarto.cmd"))
        )
      } else configured
    )
  }
  candidates <- c(candidates, Sys.which("quarto"))
  conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
  if (nzchar(conda_prefix)) candidates <- c(
    candidates,
    file.path(conda_prefix, "bin", c("quarto", "quarto.exe", "quarto.cmd")),
    file.path(conda_prefix, "Scripts", c("quarto.exe", "quarto.cmd")),
    file.path(conda_prefix, "Library", "bin", c("quarto.exe", "quarto.cmd"))
  )
  r_prefix <- normalizePath(file.path(R.home(), "..", ".."), mustWork = FALSE)
  candidates <- unique(c(
    candidates,
    file.path(r_prefix, "bin", c("quarto", "quarto.exe", "quarto.cmd"))
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates) & file.exists(candidates)]
  if (length(candidates)) normalizePath(candidates[[1L]], winslash = "/") else ""
}

.system_command_version <- function(command, arguments = "--version") {
  if (!nzchar(command) || !file.exists(command)) return(NA_character_)
  output <- tryCatch(
    suppressWarnings(system2(command, arguments, stdout = TRUE, stderr = TRUE)),
    error = function(e) character()
  )
  token <- regmatches(output, regexpr("[0-9]+(?:[.][0-9]+){1,3}", output, perl = TRUE))
  token <- token[nzchar(token)]
  if (length(token)) token[[1L]] else NA_character_
}

.jags_state <- function(minimum_version, infercnv_ready = FALSE) {
  paths <- unname(Sys.which(c("jags", "jags-terminal")))
  paths <- paths[!is.na(paths) & nzchar(paths)]
  path <- if (length(paths)) paths[[1L]] else ""
  version <- .system_command_version(path)
  detected <- isTRUE(infercnv_ready) || nzchar(path)
  version_ok <- !detected || is.na(version) ||
    utils::compareVersion(version, minimum_version) >= 0L
  list(
    installed = detected,
    loadable = isTRUE(infercnv_ready),
    version = version,
    version_ok = version_ok,
    ready = detected && version_ok,
    detail = if (detected && version_ok) {
      if (isTRUE(infercnv_ready)) "infercnv loaded successfully with its rjags/JAGS dependency." else paste0("JAGS executable found at ", path, ".")
    } else if (detected) {
      paste0("Detected JAGS ", version, ", but version ", minimum_version, " or newer is required.")
    } else {
      "JAGS was not detected. Install JAGS 4.x, then reinstall/load rjags before infercnv."
    }
  )
}

.hdf5_state <- function(monocle_ready = FALSE) {
  if (isTRUE(monocle_ready)) {
    return(list(
      installed = TRUE, loadable = TRUE, version = NA_character_, version_ok = TRUE,
      ready = TRUE,
      detail = "monocle3 is already loadable; its BPCells/HDF5 installation prerequisite does not need to be rebuilt."
    ))
  }
  h5cc <- Sys.which("h5cc")
  pkg_config <- Sys.which("pkg-config")
  pkg_ready <- FALSE
  if (nzchar(pkg_config)) {
    status <- tryCatch(system2(pkg_config, c("--exists", "hdf5")), error = function(e) 1L)
    pkg_ready <- identical(as.integer(status), 0L)
  }
  ready <- nzchar(h5cc) || pkg_ready
  list(
    installed = ready, loadable = ready, version = NA_character_, version_ok = TRUE,
    ready = ready,
    detail = if (ready) {
      "HDF5 development files were detected for BPCells compilation."
    } else {
      "HDF5 development files were not detected. Install the HDF5 library and pkg-config before current Monocle 3/BPCells."
    }
  )
}

.dependency_system_state <- function(component, minimum_version, package_states) {
  if (identical(component, "R")) {
    version <- as.character(getRversion())
    version_ok <- utils::compareVersion(version, minimum_version) >= 0L
    return(list(
      installed = TRUE, loadable = TRUE, version = version, version_ok = version_ok,
      ready = version_ok,
      detail = if (version_ok) "R runtime version is compatible." else paste0("R ", minimum_version, " or newer is required.")
    ))
  }
  if (identical(component, "Quarto CLI")) {
    path <- .find_quarto_cli()
    version <- .system_command_version(path)
    version_ok <- nzchar(path) && !is.na(version) &&
      utils::compareVersion(version, minimum_version) >= 0L
    return(list(
      installed = nzchar(path), loadable = nzchar(path), version = version,
      version_ok = version_ok, ready = version_ok,
      detail = if (version_ok) {
        paste0("Quarto CLI ", version, " found at ", path, ".")
      } else if (nzchar(path)) {
        paste0("Quarto CLI was found at ", path, " but its version could not be verified as >= ", minimum_version, ".")
      } else {
        "Quarto CLI was not found. The R package named 'quarto' does not install the external CLI."
      }
    ))
  }
  if (identical(component, "JAGS")) {
    infercnv_ready <- isTRUE(package_states[["infercnv"]]$ready %||% FALSE)
    return(.jags_state(minimum_version, infercnv_ready))
  }
  if (identical(component, "HDF5 development libraries")) {
    monocle_ready <- isTRUE(package_states[["monocle3"]]$ready %||% FALSE)
    return(.hdf5_state(monocle_ready))
  }
  list(
    installed = FALSE, loadable = FALSE, version = NA_character_, version_ok = FALSE,
    ready = FALSE, detail = "Unknown system dependency."
  )
}

.dependency_status_from_request <- function(request) {
  catalog <- request$catalog
  package_states <- list()
  package_rows <- which(catalog$kind == "r_package")
  for (index in package_rows) {
    package_states[[catalog$component[[index]]]] <- .dependency_package_state(
      catalog$component[[index]], catalog$minimum_version[[index]]
    )
  }
  rows <- lapply(seq_len(nrow(catalog)), function(index) {
    entry <- catalog[index, , drop = FALSE]
    state <- if (identical(entry$kind[[1L]], "r_package")) {
      package_states[[entry$component[[1L]]]]
    } else {
      .dependency_system_state(
        entry$component[[1L]], entry$minimum_version[[1L]], package_states
      )
    }
    data.frame(
      component = entry$component,
      kind = entry$kind,
      source = entry$source,
      modules = gsub("|", ", ", entry$modules, fixed = TRUE),
      species = gsub("|", ", ", entry$species, fixed = TRUE),
      required = entry$required,
      installed = isTRUE(state$installed),
      loadable = isTRUE(state$loadable),
      version = state$version %||% NA_character_,
      minimum_version = entry$minimum_version,
      ready = isTRUE(state$ready),
      status = if (isTRUE(state$ready)) {
        "ready"
      } else if (!isTRUE(entry$required)) {
        "recommended_missing"
      } else if (isTRUE(state$installed) && !isTRUE(state$version_ok)) {
        "outdated"
      } else if (isTRUE(state$installed) && !isTRUE(state$loadable)) {
        "broken"
      } else if (identical(entry$kind[[1L]], "system")) {
        "manual_action"
      } else {
        "missing"
      },
      detail = state$detail %||% entry$reason,
      install_spec = entry$install_spec,
      reason = entry$reason,
      stringsAsFactors = FALSE
    )
  })
  status <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (isTRUE(request$unresolved_species)) {
    status <- rbind(status, data.frame(
      component = "Species-specific resources",
      kind = "configuration",
      source = "user_input",
      modules = paste(
        intersect(request$modules, .species_dependency_modules),
        collapse = ", "
      ),
      species = "unresolved",
      required = TRUE,
      installed = FALSE,
      loadable = FALSE,
      version = NA_character_,
      minimum_version = NA_character_,
      ready = FALSE,
      status = "needs_species",
      detail = request$species$message,
      install_spec = NA_character_,
      reason = "Organism databases are selected only after a safe species choice.",
      stringsAsFactors = FALSE
    ))
  }
  rownames(status) <- NULL
  structure(
    status,
    class = c("scRDSreport_dependency_status", "data.frame"),
    request = request[c("profile", "modules", "species", "render")]
  )
}

#' Inspect dependencies for a report profile
#'
#' `dependency_status()` checks both package presence and whether each package
#' namespace can actually load. It also treats the Quarto CLI, JAGS, and the
#' HDF5 build prerequisite as external system dependencies rather than
#' confusing them with similarly named R packages.
#'
#' @param profile One of `"full"`, `"core"`, or `"report_only"`.
#' @param modules Optional module selection in the same formats accepted by
#'   [report_config()].
#' @param species One or more supported species, `"auto"`, or `"all"`.
#'   Auto detection needs `input`; unresolved species are never guessed.
#' @param input Optional RDS path used only to detect the species when
#'   `species = "auto"`.
#' @param render Whether the dependency plan must include HTML rendering.
#' @return A data frame with installation, namespace-load, version, and system
#'   readiness information.
#' @export
dependency_status <- function(profile = c("full", "core", "report_only"),
                              modules = NULL, species = "auto", input = NULL,
                              render = TRUE) {
  request <- .dependency_request(profile, modules, species, input, render)
  .dependency_status_from_request(request)
}

.dependency_failure_message <- function(status, install_errors = character()) {
  missing <- status[status$required & !status$ready, , drop = FALSE]
  lines <- if (nrow(missing)) paste0(
    "- ", missing$component, " [", missing$status, "]: ", missing$detail
  ) else character()
  errors <- if (length(install_errors)) paste0("- installer: ", install_errors) else character()
  paste(
    c(
      "scRDSreport dependencies are not ready:",
      lines,
      errors,
      "Run install_dependencies() after resolving the listed system/configuration requirements."
    ),
    collapse = "\n"
  )
}

#' Validate dependencies for a report profile
#'
#' @inheritParams dependency_status
#' @param strict Stop when a required dependency is not ready. If `FALSE`, a
#'   warning is emitted instead.
#' @param quiet Suppress the success message.
#' @return Invisibly, the dependency status data frame.
#' @export
check_dependencies <- function(profile = c("full", "core", "report_only"),
                               modules = NULL, species = "auto", input = NULL,
                               render = TRUE, strict = TRUE, quiet = FALSE) {
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict) ||
      !is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("strict and quiet must be TRUE or FALSE.", call. = FALSE)
  }
  status <- dependency_status(profile, modules, species, input, render)
  missing <- status$required & !status$ready
  if (any(missing)) {
    message <- .dependency_failure_message(status)
    if (isTRUE(strict)) stop(message, call. = FALSE) else warning(message, call. = FALSE)
  } else if (!isTRUE(quiet)) {
    message("All required scRDSreport dependencies for this request are ready.")
  }
  invisible(status)
}

.dependency_install_plan <- function(status) {
  action <- rep("none", nrow(status))
  action[status$kind == "r_package" & !status$ready] <- paste0("install_", status$source[status$kind == "r_package" & !status$ready])
  action[status$kind == "system" & !status$ready] <- "manual_system_install"
  action[status$kind == "configuration" & !status$ready] <- "specify_species_or_resource"
  output <- status
  output$action <- action
  output
}

.bootstrap_dependency_managers <- function(lib) {
  missing <- c("BiocManager", "remotes")[
    !vapply(c("BiocManager", "remotes"), requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing)) {
    utils::install.packages(missing, repos = "https://cloud.r-project.org", lib = lib)
  }
  unavailable <- c("BiocManager", "remotes")[
    !vapply(c("BiocManager", "remotes"), requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(unavailable)) {
    stop("Could not install dependency managers: ", paste(unavailable, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.default_dependency_library <- function() {
  paths <- .libPaths()
  writable <- paths[dir.exists(paths) & file.access(paths, 2L) == 0L]
  if (length(writable)) return(writable[[1L]])

  configured <- Sys.getenv("R_LIBS_USER", unset = "")
  configured <- strsplit(configured, .Platform$path.sep, fixed = TRUE)[[1L]]
  configured <- configured[nzchar(trimws(configured))]
  if (length(configured)) return(path.expand(configured[[1L]]))

  r_version <- paste0(R.version$major, ".", sub("[.].*$", "", R.version$minor))
  file.path(path.expand("~"), "R", paste0(R.version$platform, "-library"), r_version)
}

.install_monocle3 <- function(lib, repos, upgrade, verbose,
                              minimum_version = "1.0.0") {
  if (.dependency_package_state("monocle3", minimum_version)$ready) {
    return(invisible(TRUE))
  }
  if (getRversion() < "4.4.1") {
    stop(
      "Current Monocle 3 installation requires R >= 4.4.1. Upgrade R, or install a compatible older Monocle 3 manually.",
      call. = FALSE
    )
  }
  bioc_version <- tryCatch(BiocManager::version(), error = function(e) numeric_version("0"))
  if (bioc_version < "3.21") {
    stop("Current Monocle 3 installation requires Bioconductor >= 3.21.", call. = FALSE)
  }
  prerequisites <- c(
    "BiocGenerics", "DelayedArray", "DelayedMatrixStats", "limma", "lme4",
    "S4Vectors", "SingleCellExperiment", "SummarizedExperiment", "batchelor",
    "HDF5Array", "ggrastr"
  )
  BiocManager::install(
    prerequisites, ask = FALSE, update = isTRUE(upgrade), lib = lib, quiet = !isTRUE(verbose)
  )
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    remotes::install_github(
      "bnprks/BPCells", subdir = "r", dependencies = NA,
      upgrade = if (isTRUE(upgrade)) "always" else "never",
      repos = repos, lib = lib, quiet = !isTRUE(verbose)
    )
  }
  remotes::install_github(
    "cole-trapnell-lab/monocle3", dependencies = NA,
    upgrade = if (isTRUE(upgrade)) "always" else "never",
    repos = repos, lib = lib, quiet = !isTRUE(verbose)
  )
  invisible(TRUE)
}

.install_requested_dependencies <- function(status, lib, upgrade, verbose) {
  failures <- character()
  missing <- status[status$kind == "r_package" & !status$ready, , drop = FALSE]
  if (!nrow(missing)) return(failures)
  standard <- unique(missing$component[missing$source %in% c("cran", "bioconductor")])
  github <- missing[missing$source == "github", , drop = FALSE]
  repos <- BiocManager::repositories()
  old_repos <- getOption("repos")
  options(repos = repos)
  on.exit(options(repos = old_repos), add = TRUE)

  if (length(standard)) {
    tryCatch(
      BiocManager::install(
        standard, ask = FALSE, update = isTRUE(upgrade), lib = lib,
        quiet = !isTRUE(verbose)
      ),
      error = function(e) failures <<- c(failures, paste0("CRAN/Bioconductor batch: ", conditionMessage(e)))
    )
  }
  if (nrow(github)) {
    for (index in seq_len(nrow(github))) {
      package <- github$component[[index]]
      remote <- github$install_spec[[index]]
      tryCatch(
        {
          if (identical(package, "monocle3")) {
            .install_monocle3(
              lib, repos, upgrade, verbose,
              minimum_version = github$minimum_version[[index]]
            )
          } else {
            remotes::install_github(
              remote, dependencies = NA,
              upgrade = if (isTRUE(upgrade)) "always" else "never",
              repos = repos, lib = lib, quiet = !isTRUE(verbose)
            )
          }
        },
        error = function(e) failures <<- c(failures, paste0(package, ": ", conditionMessage(e)))
      )
    }
  }
  unique(failures)
}

#' Install dependencies for a report profile
#'
#' The installer uses Bioconductor repositories for CRAN/Bioconductor packages,
#' installs declared GitHub-only analysis engines separately, and then reruns
#' namespace-load and version checks. External programs such as Quarto, JAGS,
#' and HDF5 development libraries are reported with actionable instructions but
#' are not silently installed by an R package.
#'
#' @inheritParams dependency_status
#' @param upgrade Whether to update already installed packages. The default is
#'   `FALSE` to avoid changing unrelated working environments.
#' @param lib R library in which missing packages should be installed.
#'   `NULL` selects the first writable active library, then falls back to the
#'   configured user library.
#' @param dry_run Return the installation plan without changing the library.
#' @param strict Stop after installation if any required package, system
#'   dependency, or species selection is still unresolved.
#' @param verbose Show installer progress.
#' @return Invisibly, the final dependency status; with `dry_run = TRUE`, the
#'   returned data frame has an additional `action` column.
#' @export
install_dependencies <- function(profile = c("full", "core", "report_only"),
                                 modules = NULL, species = "auto", input = NULL,
                                 render = TRUE, upgrade = FALSE,
                                 lib = NULL, dry_run = FALSE,
                                 strict = TRUE, verbose = TRUE) {
  if (!is.logical(upgrade) || length(upgrade) != 1L || is.na(upgrade) ||
      !is.logical(dry_run) || length(dry_run) != 1L || is.na(dry_run) ||
      !is.logical(strict) || length(strict) != 1L || is.na(strict) ||
      !is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("upgrade, dry_run, strict, and verbose must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(lib)) lib <- .default_dependency_library()
  if (!is.character(lib) || length(lib) != 1L || is.na(lib) || !nzchar(trimws(lib))) {
    stop("lib must be one non-empty R library path.", call. = FALSE)
  }
  request <- .dependency_request(profile, modules, species, input, render)
  old_paths <- .libPaths()
  restore_paths <- TRUE
  on.exit(if (isTRUE(restore_paths)) .libPaths(old_paths), add = TRUE)
  if (dir.exists(lib)) {
    .libPaths(unique(c(normalizePath(lib, winslash = "/"), old_paths)))
  }
  before <- .dependency_status_from_request(request)
  plan <- .dependency_install_plan(before)
  if (isTRUE(dry_run)) {
    return(invisible(plan))
  }

  needs_r_install <- any(before$kind == "r_package" & !before$ready)
  failures <- character()
  if (needs_r_install) {
    dir.create(lib, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(lib) || file.access(lib, 2L) != 0L) {
      stop("R library is not writable: ", lib, call. = FALSE)
    }
    .libPaths(unique(c(normalizePath(lib, winslash = "/"), .libPaths())))
    .bootstrap_dependency_managers(lib)
    failures <- .install_requested_dependencies(before, lib, upgrade, verbose)
  }

  after <- .dependency_status_from_request(request)
  missing <- after$required & !after$ready
  recommended_missing <- !after$required & !after$ready
  if (any(missing)) {
    message <- .dependency_failure_message(after, failures)
    if (isTRUE(strict)) stop(message, call. = FALSE) else warning(message, call. = FALSE)
  } else if (any(recommended_missing) || length(failures)) {
    recommended <- after[recommended_missing, , drop = FALSE]
    recommended_lines <- if (nrow(recommended)) paste0(
      "- ", recommended$component, " [", recommended$status, "]: ",
      recommended$detail
    ) else character()
    warning(
      paste(
        c(
          "All required dependencies are ready, but some recommended dependencies remain unavailable:",
          recommended_lines,
          if (length(failures)) paste0("- installer: ", failures) else character()
        ),
        collapse = "\n"
      ),
      call. = FALSE
    )
  } else if (isTRUE(verbose)) {
    message("All requested scRDSreport R dependencies are installed and loadable.")
  }
  restore_paths <- FALSE
  invisible(after)
}
