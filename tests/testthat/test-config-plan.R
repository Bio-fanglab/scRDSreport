test_that("report_config defaults to the complete full profile", {
  config <- report_config()
  expect_s3_class(config, "scRDSreport_config")
  expect_equal(config$profile, "full")
  expect_equal(names(config$modules), scRDSreport:::.report_module_ids())
  expect_true(all(config$modules))
  expect_equal(config$annotation$mode, "auto_if_missing")
  expect_equal(config$differential$strategy, "auto")
  expect_null(config$trajectory$root)
  expect_null(config$cnv$reference_groups)
})

test_that("module_options are validated and exposed to module runners", {
  config <- report_config(module_options = list(
    qc = list(filter = TRUE, max_percent_mt = 15),
    pseudotime = list(max_cells = 5000)
  ))
  expect_true(config$qc$filter)
  expect_equal(config$qc$max_percent_mt, 15)
  expect_equal(config$module_options$pseudotime$max_cells, 5000)
  merged <- report_config(
    differential = "wilcox",
    module_options = list(differential = list(max_contrasts = 3L))
  )
  expect_equal(merged$differential$strategy, "wilcox")
  expect_equal(merged$differential$max_contrasts, 3L)
  expect_error(
    report_config(module_options = list(unknown = list(enabled = TRUE))),
    "Unknown module option target"
  )
})

test_that("auto_if_missing preserves existing annotations and otherwise requires a matched reference", {
  config <- report_config()
  mouse_with_annotation <- scRDSreport:::.build_analysis_plan(
    config,
    context = list(
      species = "mouse", has_annotation = TRUE, has_clusters = TRUE,
      n_samples = 2L, n_cells = 500L
    )
  )
  expect_equal(
    scRDSreport:::.module_by_id(mouse_with_annotation, "celltype")$reason,
    "ready"
  )

  mouse_without_annotation <- scRDSreport:::.build_analysis_plan(
    config,
    context = list(
      species = "mouse", has_annotation = FALSE, has_clusters = TRUE,
      n_samples = 2L, n_cells = 500L
    )
  )
  expect_true(scRDSreport:::.module_by_id(mouse_without_annotation, "celltype")$eligible)

  rat_without_annotation <- scRDSreport:::.build_analysis_plan(
    config,
    context = list(
      species = "rat", has_annotation = FALSE, has_clusters = TRUE,
      n_samples = 2L, n_cells = 500L
    )
  )
  expect_equal(
    scRDSreport:::.module_by_id(rat_without_annotation, "celltype")$reason,
    "annotation_resources_unavailable"
  )

  rat_reference <- report_config(
    module_options = list(celltype = list(reference = structure(list(dummy = TRUE), class = "synthetic_reference")))
  )
  rat_with_reference <- scRDSreport:::.build_analysis_plan(
    rat_reference,
    context = list(
      species = "rat", has_annotation = FALSE, has_clusters = TRUE,
      n_samples = 2L, n_cells = 500L
    )
  )
  expect_true(scRDSreport:::.module_by_id(rat_with_reference, "celltype")$eligible)
})

test_that("profiles and explicit module selections retain every plan entry", {
  core <- report_config(profile = "core")
  expect_true(core$modules[["qc"]])
  expect_false(core$modules[["enrichment"]])
  expect_true(core$modules[["downloads"]])

  selected <- report_config(profile = "report_only", modules = c("enrichment"))
  expect_true(selected$modules[["enrichment"]])
  expect_true(selected$modules[["downloads"]])
  expect_false(selected$modules[["qc"]])

  plan <- scRDSreport:::.build_analysis_plan(core)
  expect_setequal(names(plan$modules), scRDSreport:::.report_module_ids())
  expect_equal(scRDSreport:::.module_by_id(plan, "enrichment")$status, "not_requested")
  expect_error(report_config(modules = "not_a_module"), "Unknown report module")
})

test_that("biological choices and resource limits are validated", {
  config <- report_config(
    annotation_mode = "manual",
    differential = "pseudobulk",
    trajectory_root = list(column = "state", value = "naive"),
    cnv_reference = c("T cells", "B cells"),
    resource_overrides = list(gtf = "genome.gtf", manual_markers = list(T = "CD3D")),
    limits = list(workers = 4, analysis_max_cells = 20000)
  )
  expect_equal(config$annotation$mode, "manual")
  expect_equal(config$trajectory$root$value, "naive")
  expect_equal(config$cnv$reference_groups, c("T cells", "B cells"))
  expect_equal(config$limits$workers, 4L)
  expect_equal(config$limits$analysis_max_cells, 20000)

  node_config <- report_config(trajectory_root = "Y_12")
  expect_equal(node_config$trajectory$root$type, "principal_node")
  expect_equal(node_config$trajectory$root$value, "Y_12")
  expect_error(report_config(limits = list(workers = 0)), "workers")
  expect_error(report_config(limits = list(typo = 1)), "Unknown limit")
  expect_error(report_config(cnv_reference = ""), "cannot be empty")
})

test_that("human and mouse resources are explicit and unknown never falls back to mouse", {
  human <- species_resources("Homo sapiens")
  expect_s3_class(human, "scRDSreport_species_resources")
  expect_true(human$supported)
  expect_equal(human$orgdb, "org.Hs.eg.db")
  expect_equal(human$kegg_code, "hsa")
  expect_equal(human$cellchat_db, "CellChatDB.human")
  expect_equal(human$ensembl_prefix, "^ENSG[0-9]")
  expect_equal(length(human$chromosomes), 25L)

  mouse <- species_resources("mouse")
  expect_true(mouse$supported)
  expect_equal(mouse$orgdb, "org.Mm.eg.db")
  expect_equal(mouse$kegg_code, "mmu")
  expect_equal(mouse$cellchat_db, "CellChatDB.mouse")
  expect_equal(mouse$txdb, "TxDb.Mmusculus.UCSC.mm10.knownGene")
  expect_equal(mouse$ensembl_prefix, "^ENSMUSG[0-9]")
  expect_equal(length(mouse$chromosomes), 22L)

  unknown <- species_resources("platypus")
  expect_false(unknown$supported)
  expect_null(unknown$orgdb)
  expect_null(unknown$kegg_code)
  expect_null(unknown$cellchat_db)
  expect_null(unknown$txdb)
  expect_false(any(grepl("mouse|Mm.eg|mm10", unlist(unknown), ignore.case = TRUE)))

  custom <- species_resources(
    "platypus",
    overrides = list(orgdb = "org.Custom.eg.db", gtf = "custom.gtf", kegg_code = "custom")
  )
  expect_true(custom$supported)
  expect_equal(custom$source, "user_override")
  expect_equal(custom$orgdb, "org.Custom.eg.db")
})

test_that("full plans contain every original report module and explain unmet prerequisites", {
  config <- report_config()
  plan <- scRDSreport:::.build_analysis_plan(
    config,
    context = list(
      species = "mouse", has_annotation = TRUE, has_clusters = TRUE,
      n_samples = 4L, n_cells = 1000L
    )
  )
  expect_s3_class(plan, "scRDSreport_analysis_plan")
  expect_equal(names(plan$modules), scRDSreport:::.report_module_ids())
  expect_true(all(vapply(plan$modules, `[[`, logical(1), "requested")))
  expect_identical(names(scRDSreport:::.module_by_id(plan, "qc")),
                   scRDSreport:::.module_schema_fields())
  expect_equal(scRDSreport:::.module_by_id(plan, "qc")$status, "planned")
  expect_true(scRDSreport:::.module_by_id(plan, "pseudotime")$eligible)
  expect_equal(scRDSreport:::.module_by_id(plan, "pseudotime")$reason, "root_candidate_only")
  expect_equal(scRDSreport:::.module_by_id(plan, "cnv")$reason, "cnv_reference_required")

  configured <- report_config(
    trajectory_root = "Y_1",
    cnv_reference = c("T cells", "B cells"),
    module_options = list(
      cnv = list(object_genome_assembly = "GRCm38/mm10")
    )
  )
  ready <- scRDSreport:::.build_analysis_plan(
    configured,
    context = list(
      species = "mouse", has_annotation = TRUE, has_clusters = TRUE,
      n_samples = 4L, n_cells = 1000L
    )
  )
  expect_true(scRDSreport:::.module_by_id(ready, "pseudotime")$eligible)
  expect_true(scRDSreport:::.module_by_id(ready, "cnv")$eligible)
})

test_that("CNV plans require genome-build confirmation before using a built-in TxDb", {
  empty <- scRDSreport:::.build_analysis_plan(
    report_config(module_options = list(cnv = list(
      reference_groups = character(), gene_order = character()
    ))),
    context = list(
      species = "mouse", has_annotation = TRUE, has_clusters = TRUE,
      n_samples = 2L, n_cells = 500L
    )
  )
  expect_equal(scRDSreport:::.module_by_id(empty, "cnv")$reason, "cnv_reference_required")

  unconfirmed <- scRDSreport:::.build_analysis_plan(
    report_config(cnv_reference = "Normal"),
    context = list(
      species = "mouse", has_annotation = TRUE, has_clusters = TRUE,
      n_samples = 2L, n_cells = 500L
    )
  )
  cnv <- scRDSreport:::.module_by_id(unconfirmed, "cnv")
  expect_false(cnv$eligible)
  expect_equal(cnv$reason, "cnv_genome_assembly_confirmation_required")

  confirmed <- scRDSreport:::.build_analysis_plan(
    report_config(
      cnv_reference = "Normal",
      module_options = list(cnv = list(object_genome_assembly = "mm10"))
    ),
    context = list(
      species = "mouse", has_annotation = TRUE, has_clusters = TRUE,
      n_samples = 2L, n_cells = 500L
    )
  )
  expect_true(scRDSreport:::.module_by_id(confirmed, "cnv")$eligible)
})

test_that("unknown species leave species-dependent modules visible but ineligible", {
  plan <- scRDSreport:::.build_analysis_plan(
    report_config(),
    context = list(
      species = "unknown", has_annotation = TRUE, has_clusters = TRUE,
      n_samples = 3L, n_cells = 500L
    )
  )
  expect_equal(scRDSreport:::.module_by_id(plan, "enrichment")$reason,
               "enrichment_resources_unavailable")
  expect_equal(scRDSreport:::.module_by_id(plan, "communication")$reason,
               "communication_context_only")
  expect_equal(scRDSreport:::.module_by_id(plan, "cell_cycle")$reason,
               "cell_cycle_resources_unavailable")
})

test_that("artifact records use the complete typed schema", {
  empty_table <- scRDSreport:::.artifact_table(scRDSreport:::.new_artifact_registry())
  expect_identical(names(empty_table), scRDSreport:::.artifact_schema_fields())
  expect_type(empty_table$rows, "double")
  expect_type(empty_table$complete, "logical")

  path <- tempfile(fileext = ".csv")
  writeLines(c("cell,value", "c1,1"), path)
  artifact <- scRDSreport:::.new_artifact(
    module = "qc",
    type = "table",
    path = path,
    label = "Cell QC",
    description = "One row per cell.",
    rows = 1,
    columns = 2,
    row_unit = "cell",
    column_unit = "field",
    column_dictionary = c(cell = "Cell identifier", value = "QC value"),
    preview = "tables/qc_preview.csv",
    embed_priority = 10
  )
  expect_s3_class(artifact, "scRDSreport_artifact")
  expect_identical(names(artifact), scRDSreport:::.artifact_schema_fields())
  expect_true(artifact$complete)
  expect_gt(artifact$bytes, 0)
  if (!is.na(artifact$sha256)) expect_match(artifact$sha256, "^[0-9a-f]{64}$")
  expect_equal(artifact$column_dictionary$cell, "Cell identifier")

  registry <- scRDSreport:::.new_artifact_registry()
  registry <- scRDSreport:::.register_artifact(registry, artifact)
  expect_equal(scRDSreport:::.artifact_ids(registry), artifact$id)
  table <- scRDSreport:::.artifact_table(registry)
  expect_identical(names(table), scRDSreport:::.artifact_schema_fields())
  expect_error(scRDSreport:::.register_artifact(registry, artifact), "already registered")
})

test_that("safe module runner captures warnings, artifacts, and optional failures", {
  plan <- scRDSreport:::.build_analysis_plan(report_config())
  qc <- scRDSreport:::.module_by_id(plan, "qc")
  path <- tempfile(fileext = ".txt")
  writeLines("ok", path)
  artifact <- scRDSreport:::.new_artifact(
    module = "qc", type = "table", path = path,
    label = "QC output", description = "Synthetic output"
  )

  success <- scRDSreport:::.run_module_safely(
    qc,
    function() {
      warning("captured warning")
      scRDSreport:::.module_result(value = 42, artifacts = artifact)
    },
    engine = "base::identity",
    seed = 123
  )
  expect_equal(success$record$status, "completed")
  expect_equal(success$value, 42)
  expect_match(success$record$warnings, "captured warning")
  expect_equal(success$record$artifact_ids, artifact$id)
  expect_true(is.numeric(success$record$timing$elapsed_seconds))

  failure <- scRDSreport:::.run_module_safely(
    qc,
    function() stop("module boom")
  )
  expect_equal(failure$record$status, "failed")
  expect_match(failure$record$error, "module boom")
  expect_null(failure$value)

  skipped_record <- qc
  skipped_record$eligible <- FALSE
  skipped_record$reason <- "missing_input"
  touched <- FALSE
  skipped <- scRDSreport:::.run_module_safely(skipped_record, function() touched <<- TRUE)
  expect_equal(skipped$record$status, "skipped")
  expect_false(touched)
})
