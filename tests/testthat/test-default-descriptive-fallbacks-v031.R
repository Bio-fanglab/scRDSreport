test_that("default full configuration enables bounded descriptive outputs", {
  config <- report_config()
  expect_equal(config$differential$strategy, "auto")
  expect_equal(config$differential$fallback_grouping, c("group", "annotation", "cluster"))
  expect_equal(config$differential$max_contrasts, 6L)
  expect_true(config$differential$run_group_markers)
  expect_true(config$enrichment$descriptive_rank_summary)
  expect_true(config$trajectory$export_geometry_without_root)
  expect_true(config$communication$export_group_overview)
  expect_true(config$cnv$export_readiness)
})

test_that("automatic DE with an unreviewed design returns visible descriptive results", {
  object <- .release_v020_object(n_features = 60L, n_cells = 40L)
  object$group <- NULL
  inferred <- infer_sample_design(rep(c("WT_1", "WT_2", "KO_1", "KO_2"), each = 10L))
  object <- scRDSreport:::.attach_design(object, inferred)
  output <- tempfile("scrdsreport-default-de-")
  result <- scRDSreport:::.fa_module_differential(
    object, output,
    scRDSreport:::.fa_module_config(report_config(), "differential"),
    seed = 7L, verbose = FALSE
  )

  expect_equal(result$status, "partial")
  expect_equal(result$reason_code, "sample_design_needs_review_descriptive_completed")
  expect_false(result$details$inferential)
  expect_true(result$details$exploratory)
  rankings <- scRDSreport:::.fa_get_analysis_misc(result$object, "differential_rankings")
  expect_gt(nrow(rankings), 0L)
  expect_false(any(c("PValue", "FDR") %in% names(rankings)))
  expect_true(any(rankings$analysis_type == "exploratory_one_vs_rest_effect_size"))
  expect_true(any(vapply(result$artifacts, function(x) identical(x$type, "figure"), logical(1))))
})

test_that("annotation fallback produces one-vs-rest markers without pairwise explosion", {
  object <- .release_v020_object(n_features = 60L, n_cells = 120L)
  object$group <- NULL
  object$sample <- NULL
  object$celltype_manual <- NULL
  object$seurat_clusters <- rep(as.character(0:11), each = 10L)
  output <- tempfile("scrdsreport-cluster-markers-")
  cfg <- scRDSreport:::.fa_module_config(report_config(), "differential")
  result <- scRDSreport:::.fa_module_differential(
    object, output, cfg, seed = 8L, verbose = FALSE
  )

  expect_equal(result$status, "partial")
  expect_equal(result$details$grouping_source, "cluster_fallback")
  rankings <- scRDSreport:::.fa_get_analysis_misc(result$object, "differential_rankings")
  expect_setequal(unique(rankings$analysis_type), "exploratory_one_vs_rest_effect_size")
  expect_equal(length(unique(rankings$contrast)), 12L)
  expect_lte(length(unique(rankings$contrast)), cfg$max_marker_groups)
  expect_false(any(c("PValue", "FDR") %in% names(rankings)))
})

test_that("exploratory rankings feed descriptive gene-set summaries without P values", {
  rankings <- data.frame(
    feature = paste0("g", 1:8), SYMBOL = paste0("G", 1:8),
    logFC = c(2, 1.5, 1, -2, -1.5, -1, 0.5, -0.5),
    stratum = "all_cells", contrast = "B vs A",
    analysis_type = "exploratory_effect_size_only",
    stringsAsFactors = FALSE
  )
  sets <- list(up_set = c("G1", "G2", "G3"), down_set = c("G4", "G5", "G6"))
  result <- scRDSreport:::.fa_descriptive_gene_set_effects(rankings, sets)
  expect_equal(nrow(result), 2L)
  expect_true(all(result$analysis_type == "descriptive_gene_set_effect_summary"))
  expect_false(any(grepl("^p$|pvalue|fdr|padj", names(result), ignore.case = TRUE)))
  expect_gt(result$mean_logFC[result$gene_set == "up_set"], 0)
  expect_lt(result$mean_logFC[result$gene_set == "down_set"], 0)
})

test_that("descriptive contrasts name pooled CPM accurately", {
  counts <- Matrix::Matrix(
    matrix(c(10, 0, 20, 0, 0, 5, 0, 15), nrow = 2L,
           dimnames = list(c("g1", "g2"), paste0("s", 1:4))),
    sparse = TRUE
  )
  result <- scRDSreport:::.fa_exploratory_contrast(
    counts, c("A", "A", "B", "B"), "A", "B"
  )
  expect_true(all(c("pooled_cpm_group_a", "pooled_cpm_group_b") %in% names(result)))
  expect_false(any(c("mean_cpm_group_a", "mean_cpm_group_b") %in% names(result)))
})

test_that("composite analysis engines retain installed package versions", {
  testthat::local_mocked_bindings(
    .fa_package_version = function(package) {
      versions <- c(clusterProfiler = "4.16.0", GSVA = "2.2.0")
      if (package %in% names(versions)) unname(versions[[package]]) else NA_character_
    },
    .package = "scRDSreport"
  )
  expect_equal(
    scRDSreport:::.fa_engine_version("clusterProfiler/GSVA"),
    "clusterProfiler=4.16.0; GSVA=2.2.0"
  )
  expect_equal(scRDSreport:::.fa_engine_version("GSVA::gsva"), "2.2.0")
})

test_that("existing cell-cycle columns are exported without species resources", {
  object <- .release_v020_object(n_features = 40L, n_cells = 24L)
  object$S.Score <- seq(-1, 1, length.out = ncol(object))
  object$G2M.Score <- rev(object$S.Score)
  object$Phase <- rep(c("G1", "S", "G2M"), length.out = ncol(object))
  output <- tempfile("scrdsreport-existing-cycle-")
  result <- scRDSreport:::.fa_module_cell_cycle(
    object, output, list(), seed = 9L, verbose = FALSE, species = "unknown"
  )
  expect_equal(result$status, "completed")
  expect_equal(result$reason_code, "existing_cell_cycle_preserved")
  expect_true(any(grepl("cell_cycle_scores", vapply(result$artifacts, `[[`, character(1), "path"))))
})

test_that("CNV without a reference remains needs-input but exports readiness", {
  object <- .release_v020_object(n_features = 40L, n_cells = 24L)
  output <- tempfile("scrdsreport-cnv-readiness-")
  result <- scRDSreport:::.fa_module_cnv(
    object, output, list(export_readiness = TRUE),
    seed = 10L, verbose = FALSE, species = "unknown"
  )
  expect_equal(result$status, "needs_input")
  expect_equal(result$reason_code, "cnv_reference_missing")
  expect_true(any(grepl("cnv_readiness", vapply(result$artifacts, `[[`, character(1), "path"))))
  expect_match(result$message, "never guesses")
})

test_that("plans allow descriptive DE with one sample but never auto-select CNV references", {
  plan <- scRDSreport:::.build_analysis_plan(
    report_config(),
    context = list(
      species = "mouse", has_annotation = TRUE, has_clusters = TRUE,
      n_samples = 1L, n_cells = 500L
    )
  )
  differential <- scRDSreport:::.module_by_id(plan, "differential")
  cnv <- scRDSreport:::.module_by_id(plan, "cnv")
  expect_true(differential$eligible)
  expect_equal(differential$reason, "descriptive_fallback_ready")
  expect_false(cnv$eligible)
  expect_equal(cnv$reason, "cnv_reference_required")
})

test_that("unrooted trajectory is partial geometry and never pseudotime", {
  skip_if_not_installed("monocle3")
  object <- .release_v020_object(n_features = 60L, n_cells = 120L)
  output <- tempfile("scrdsreport-unrooted-")
  result <- scRDSreport:::.fa_module_pseudotime(
    object, output,
    list(root = NULL, max_cells = 120L, max_features = 60L, num_dim = 10L,
         run_graph_test = FALSE, export_geometry_without_root = TRUE),
    seed = 11L, verbose = FALSE
  )
  expect_equal(result$status, "partial")
  expect_equal(result$reason_code, "trajectory_geometry_completed_root_missing")
  paths <- vapply(result$artifacts, `[[`, character(1), "path")
  expect_true(any(grepl("trajectory_geometry_unrooted", paths)))
  expect_false(any(grepl("cell_pseudotime", paths)))
})

test_that("SingleR prediction artifacts use the actual fallback unit", {
  predictions <- data.frame(
    prediction_id = c("cell1", "cell2"), labels = c("T", "B"),
    stringsAsFactors = FALSE
  )
  expect_equal(scRDSreport:::.fa_prediction_row_unit(
    "auto", predictions,
    list(prediction_mode = "cell_level_fallback_missing_cluster_aggregation_dependency"),
    cluster_column = "seurat_clusters", cell_names = c("cell1", "cell2")
  ), "cell")
  expect_equal(scRDSreport:::.fa_prediction_row_unit(
    "auto", data.frame(prediction_id = c("0", "1")),
    list(prediction_mode = "cluster_level"),
    cluster_column = "seurat_clusters", cell_names = c("cell1", "cell2")
  ), "cluster")
})

test_that("CellChat creation failures retain context and diagnostics", {
  object <- .release_v020_object(n_features = 40L, n_cells = 24L)
  output <- tempfile("scrdsreport-cellchat-create-error-")
  dir.create(output)
  counts <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  testthat::local_mocked_bindings(
    .fa_pkg_available = function(package) TRUE,
    .fa_matrix = function(object, layer, assay) counts,
    .fa_species_resources = function(species, cfg) list(
      cellchat = "CellChatDB.mouse", orgdb = NULL,
      feature_keytype = NULL, symbol_column = "SYMBOL"
    ),
    .fa_pkg_object = function(package, name) list(database = name),
    .fa_orgdb = function(resources, cfg) NULL,
    .fa_feature_symbols = function(object, assay) {
      stats::setNames(rownames(counts), rownames(counts))
    },
    .fa_feature_mapping = function(features, symbols, orgdb, ...) {
      data.frame(feature = features, SYMBOL = features, stringsAsFactors = FALSE)
    },
    .fa_pkg_fun = function(package, name, exported = TRUE) {
      if (identical(package, "CellChat") && identical(name, "createCellChat")) {
        return(function(...) stop("synthetic object creation failure", call. = FALSE))
      }
      if (identical(package, "CellChat")) return(function(value, ...) value)
      NULL
    },
    .package = "scRDSreport"
  )
  result <- scRDSreport:::.fa_module_communication(
    object, output,
    list(annotation_column = "celltype_manual", min_cells = 2L,
         max_cells_per_group = 100L, export_group_overview = TRUE),
    seed = 12L, verbose = FALSE, species = "mouse"
  )
  expect_equal(result$status, "failed")
  expect_equal(result$reason_code, "cellchat_creation_failed")
  paths <- vapply(result$artifacts, `[[`, character(1), "path")
  expect_true(any(grepl("grouping_context", paths)))
  expect_true(any(grepl("cellchat_diagnostic", paths)))
  expect_true(all(file.exists(file.path(output, paths))))
})
