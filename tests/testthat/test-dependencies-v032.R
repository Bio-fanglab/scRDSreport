ready_dependency_state <- function(version = "99.0.0") {
  list(
    installed = TRUE,
    loadable = TRUE,
    version = version,
    version_ok = TRUE,
    ready = TRUE,
    detail = "mock dependency ready"
  )
}

missing_dependency_state <- function(installed = FALSE, loadable = FALSE,
                                     version_ok = TRUE) {
  list(
    installed = installed,
    loadable = loadable,
    version = if (installed) "0.0.1" else NA_character_,
    version_ok = version_ok,
    ready = FALSE,
    detail = if (installed) "mock namespace load failure" else "mock package missing"
  )
}

test_that("dependency catalogue packages are declared exactly as package dependencies", {
  description <- utils::packageDescription("scRDSreport")
  field_packages <- function(field) {
    value <- description[[field]]
    if (is.null(value) || is.na(value)) return(character())
    tokens <- trimws(unlist(strsplit(value, ",", fixed = TRUE)))
    sub("[[:space:]]*[(].*$", "", tokens)
  }
  declared <- unique(c(field_packages("Imports"), field_packages("Suggests")))
  catalog <- scRDSreport:::.dependency_catalog()
  catalog_packages <- catalog$component[catalog$kind == "r_package"]
  expect_true(all(catalog_packages %in% declared))
  expect_true(all(c("BiocManager", "remotes") %in% declared))
  expect_false(any(c("KEGG.db", "SeuratWrappers") %in% declared))
})

test_that("dependency selection is module- and species-specific", {
  mouse <- scRDSreport:::.dependency_request(
    profile = "full", species = "mouse", render = FALSE
  )$catalog
  expect_true(all(c(
    "org.Mm.eg.db", "TxDb.Mmusculus.UCSC.mm10.knownGene", "CellChat"
  ) %in% mouse$component))
  expect_false(any(c(
    "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene", "Quarto CLI"
  ) %in% mouse$component))

  rat <- scRDSreport:::.dependency_request(
    profile = "full", species = "rat", render = FALSE
  )$catalog
  expect_true(all(c("org.Rn.eg.db", "babelgene") %in% rat$component))
  expect_false(any(c("CellChat", "org.Mm.eg.db", "org.Hs.eg.db") %in% rat$component))

  all_species <- scRDSreport:::.dependency_request(
    profile = "full", species = "all", render = FALSE
  )$catalog
  expect_equal(sum(grepl("^org[.]", all_species$component)), 9L)
  expect_equal(sum(grepl("^TxDb[.]", all_species$component)), 2L)
})

test_that("auto species uses an input RDS and never guesses unresolved symbols", {
  mouse_matrix <- matrix(
    1,
    nrow = 5,
    ncol = 3,
    dimnames = list(paste0("ENSMUSG", sprintf("%011d", 1:5)), paste0("cell", 1:3))
  )
  input <- tempfile(fileext = ".rds")
  saveRDS(mouse_matrix, input)
  request <- scRDSreport:::.dependency_request(
    profile = "full", species = "auto", input = input, render = FALSE
  )
  expect_equal(request$species$species, "mouse")
  expect_false(request$unresolved_species)
  expect_true("org.Mm.eg.db" %in% request$catalog$component)

  unresolved <- scRDSreport:::.dependency_request(
    profile = "full", species = "auto", render = FALSE
  )
  expect_true(unresolved$unresolved_species)
  expect_length(unresolved$species$species, 0L)

  unresolved_core <- scRDSreport:::.dependency_request(
    profile = "core", species = "auto", render = FALSE
  )
  expect_true(unresolved_core$unresolved_species)

  unsupported <- scRDSreport:::.dependency_request(
    profile = "full", species = "unsupported_species", render = FALSE
  )
  expect_true(unsupported$unresolved_species)
  expect_length(unsupported$species$species, 0L)
})

test_that("status distinguishes broken namespaces, old versions, and unresolved species", {
  testthat::local_mocked_bindings(
    .dependency_package_state = function(package, minimum_version = NA_character_) {
      if (identical(package, "SCP")) return(missing_dependency_state(TRUE, FALSE))
      if (identical(package, "ggplot2")) return(missing_dependency_state(TRUE, TRUE, FALSE))
      ready_dependency_state()
    },
    .dependency_system_state = function(...) ready_dependency_state(),
    .package = "scRDSreport"
  )
  status <- dependency_status(
    profile = "full", species = "auto", render = FALSE
  )
  expect_equal(status$status[status$component == "SCP"], "broken")
  expect_equal(status$status[status$component == "ggplot2"], "outdated")
  expect_true("Species-specific resources" %in% status$component)
  expect_equal(
    status$status[status$component == "Species-specific resources"],
    "needs_species"
  )
})

test_that("render FALSE removes the Quarto requirements", {
  no_render <- scRDSreport:::.dependency_request(
    profile = "report_only", species = "auto", render = FALSE
  )$catalog$component
  with_render <- scRDSreport:::.dependency_request(
    profile = "report_only", species = "auto", render = TRUE
  )$catalog$component
  expect_false(any(c("quarto", "Quarto CLI") %in% no_render))
  expect_true(all(c("quarto", "Quarto CLI") %in% with_render))
})

test_that("dry-run reports source-specific actions without installing anything", {
  testthat::local_mocked_bindings(
    .dependency_package_state = function(...) missing_dependency_state(),
    .dependency_system_state = function(...) missing_dependency_state(),
    .bootstrap_dependency_managers = function(...) stop("installer must not run"),
    .install_requested_dependencies = function(...) stop("installer must not run"),
    .package = "scRDSreport"
  )
  plan <- install_dependencies(
    profile = "full", species = "mouse", render = TRUE, dry_run = TRUE
  )
  expect_equal(plan$action[plan$component == "SCP"], "install_github")
  expect_equal(
    plan$action[plan$component == "BiocParallel"],
    "install_bioconductor"
  )
  expect_equal(plan$action[plan$component == "Seurat"], "install_cran")
  expect_equal(
    plan$action[plan$component == "Quarto CLI"],
    "manual_system_install"
  )
})

test_that("strict dependency checks report every required failure", {
  testthat::local_mocked_bindings(
    .dependency_package_state = function(package, ...) {
      if (package %in% c("SCP", "Seurat")) missing_dependency_state() else ready_dependency_state()
    },
    .dependency_system_state = function(...) ready_dependency_state(),
    .package = "scRDSreport"
  )
  error <- tryCatch(
    {
      check_dependencies(
      profile = "report_only", species = "auto", render = FALSE,
      strict = TRUE, quiet = TRUE
      )
      NULL
    },
    error = identity
  )
  expect_s3_class(error, "error")
  expect_match(conditionMessage(error), "SCP", fixed = TRUE)
  expect_match(conditionMessage(error), "Seurat", fixed = TRUE)
})

test_that("standalone bootstrap validates arguments before installing", {
  bootstrap <- new.env(parent = baseenv())
  bootstrap_path <- system.file(
    "install_scRDSreport.R", package = "scRDSreport", mustWork = TRUE
  )
  sys.source(bootstrap_path, envir = bootstrap)
  expect_true(is.function(bootstrap$install_scRDSreport))
  expect_error(
    bootstrap$install_scRDSreport(profile = "invalid"),
    "arg"
  )
  expect_error(
    bootstrap$install_scRDSreport(profile = "full", species = "auto"),
    "needs input",
    fixed = TRUE
  )
  expect_error(
    bootstrap$install_scRDSreport(
      profile = "report_only", species = "auto", render = NA
    ),
    "render and upgrade"
  )
})

test_that("an already-ready installation is offline and manager independent", {
  testthat::local_mocked_bindings(
    .dependency_package_state = function(...) ready_dependency_state(),
    .dependency_system_state = function(...) ready_dependency_state(),
    .bootstrap_dependency_managers = function(...) stop("manager must not run"),
    .install_requested_dependencies = function(...) stop("installer must not run"),
    .package = "scRDSreport"
  )
  unused_library <- tempfile("scrdsreport-unused-library-")
  expect_false(dir.exists(unused_library))
  expect_no_error(
    install_dependencies(
      profile = "report_only", species = "auto", render = FALSE,
      lib = unused_library, strict = TRUE, verbose = FALSE
    )
  )
  expect_false(dir.exists(unused_library))
})

test_that("recommended dependency failures warn after required dependencies pass", {
  testthat::local_mocked_bindings(
    .dependency_package_state = function(package, ...) {
      if (identical(package, "quarto")) missing_dependency_state() else ready_dependency_state()
    },
    .dependency_system_state = function(...) ready_dependency_state(),
    .bootstrap_dependency_managers = function(...) invisible(TRUE),
    .install_requested_dependencies = function(...) "quarto: mock install failure",
    .package = "scRDSreport"
  )
  expect_warning(
    install_dependencies(
      profile = "report_only", species = "auto", render = TRUE,
      strict = TRUE, verbose = FALSE
    ),
    "recommended dependencies remain unavailable"
  )
})

test_that("recommended packages still missing after a quiet installer warning are reported", {
  testthat::local_mocked_bindings(
    .dependency_package_state = function(package, ...) {
      if (identical(package, "quarto")) missing_dependency_state() else ready_dependency_state()
    },
    .dependency_system_state = function(...) ready_dependency_state(),
    .bootstrap_dependency_managers = function(...) invisible(TRUE),
    .install_requested_dependencies = function(...) character(),
    .package = "scRDSreport"
  )
  expect_warning(
    install_dependencies(
      profile = "report_only", species = "auto", render = TRUE,
      strict = TRUE, verbose = FALSE
    ),
    "quarto"
  )
})
