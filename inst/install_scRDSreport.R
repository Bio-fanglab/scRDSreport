# Standalone bootstrap for a new R environment. Sourcing this file only defines
# the function; installation starts after the user calls install_scRDSreport().
install_scRDSreport <- function(input = NULL, profile = "full",
                                species = "auto", render = TRUE,
                                lib = NULL, upgrade = FALSE) {
  profile <- match.arg(profile, c("full", "core", "report_only"))
  if (!is.logical(render) || length(render) != 1L || is.na(render) ||
      !is.logical(upgrade) || length(upgrade) != 1L || is.na(upgrade)) {
    stop("render and upgrade must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(species) || !length(species) || anyNA(species) ||
      any(!nzchar(trimws(species)))) {
    stop("species must contain one or more non-empty species names, 'auto', or 'all'.", call. = FALSE)
  }
  special_species <- tolower(trimws(species)) %in% c("auto", "all")
  if (any(special_species) && length(species) != 1L) {
    stop("'auto' and 'all' cannot be combined with other species.", call. = FALSE)
  }
  if (!is.null(input)) {
    if (!is.character(input) || length(input) != 1L || is.na(input) ||
        !nzchar(trimws(input)) || !file.exists(input)) {
      stop("input must be NULL or one existing RDS path.", call. = FALSE)
    }
  }
  if (identical(tolower(trimws(species[[1L]])), "auto") &&
      is.null(input) && !identical(profile, "report_only")) {
    stop(
      "species = 'auto' needs input for full/core installation; supply input or an explicit species.",
      call. = FALSE
    )
  }
  if ("scRDSreport" %in% loadedNamespaces()) {
    stop(
      "scRDSreport is already loaded. Restart R and run this bootstrap in a clean session so the newly installed namespace is used.",
      call. = FALSE
    )
  }
  if (is.null(lib)) {
    paths <- .libPaths()
    writable <- paths[dir.exists(paths) & file.access(paths, 2L) == 0L]
    if (length(writable)) {
      lib <- writable[[1L]]
    } else {
      configured <- Sys.getenv("R_LIBS_USER", unset = "")
      configured <- strsplit(configured, .Platform$path.sep, fixed = TRUE)[[1L]]
      configured <- configured[nzchar(trimws(configured))]
      if (length(configured)) {
        lib <- path.expand(configured[[1L]])
      } else {
        r_version <- paste0(
          R.version$major, ".", sub("[.].*$", "", R.version$minor)
        )
        lib <- file.path(
          path.expand("~"), "R",
          paste0(R.version$platform, "-library"), r_version
        )
      }
    }
  }
  if (!is.character(lib) || length(lib) != 1L || is.na(lib) || !nzchar(trimws(lib))) {
    stop("lib must be one non-empty R library path.", call. = FALSE)
  }
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(lib) || file.access(lib, 2L) != 0L) {
    stop("R library is not writable: ", lib, call. = FALSE)
  }
  .libPaths(unique(c(normalizePath(lib, winslash = "/"), .libPaths())))

  managers <- c("BiocManager", "remotes")
  missing_managers <- managers[
    !vapply(managers, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_managers)) {
    utils::install.packages(
      missing_managers, repos = "https://cloud.r-project.org", lib = lib
    )
  }
  unavailable <- managers[
    !vapply(managers, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(unavailable)) {
    stop("Could not install: ", paste(unavailable, collapse = ", "), call. = FALSE)
  }

  repositories <- BiocManager::repositories()
  remotes::install_github(
    "Bio-fanglab/scRDSreport",
    dependencies = NA,
    upgrade = if (isTRUE(upgrade)) "always" else "never",
    repos = repositories,
    lib = lib
  )
  scRDSreport::install_dependencies(
    profile = profile,
    species = species,
    input = input,
    render = render,
    upgrade = upgrade,
    lib = lib,
    strict = TRUE
  )
  scRDSreport::check_dependencies(
    profile = profile,
    species = species,
    input = input,
    render = render,
    strict = TRUE
  )
}
