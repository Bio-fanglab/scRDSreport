test_that("public v0.2 configuration bridges module, root, and species resources", {
  config <- report_config(
    profile = "core",
    modules = list(include = "pseudotime", exclude = "celltype"),
    trajectory_root = list(type = "metadata_group", column = "state", value = "early"),
    resource_overrides = list(species = "human", gtf = "explicit.gtf"),
    module_options = list(qc = list(max_percent_mt = 12))
  )

  expect_s3_class(config, "scRDSreport_config")
  expect_true(config$modules[["pseudotime"]])
  expect_false(config$modules[["celltype"]])
  expect_true(config$modules[["downloads"]])
  expect_equal(config$qc$max_percent_mt, 12)

  pseudotime <- scRDSreport:::.fa_module_config(config, "pseudotime")
  expect_equal(pseudotime$root$type, "metadata_group")
  expect_equal(pseudotime$root$column, "state")
  expect_equal(pseudotime$root$value, "early")

  resources <- scRDSreport:::.resolve_species_resources(config, "auto")
  expect_equal(resources$species, "human")
  expect_equal(resources$gtf, "explicit.gtf")
  expect_equal(resources$orgdb, "org.Hs.eg.db")

  unknown <- species_resources("unknown")
  expect_false(unknown$supported)
  expect_null(unknown$orgdb)
  expect_null(unknown$txdb)
  expect_false(any(grepl("mouse|mm10|Mm.eg", unlist(unknown), ignore.case = TRUE)))

  rat <- species_resources("Rattus norvegicus")
  expect_equal(rat$species, "rat")
  expect_false(rat$supported)
  expect_null(rat$orgdb)
  expect_null(rat$kegg_code)
  expect_null(scRDSreport:::.fa_species_resources("rat")$orgdb)
  expect_null(scRDSreport:::.fa_species_resources("zebrafish")$kegg)
  expect_null(scRDSreport:::.fa_cell_cycle_genes("rat", list()))
  expect_length(scRDSreport:::.fa_default_tf_genes("zebrafish"), 0L)

  rat_override <- scRDSreport:::.fa_species_resources(
    "rat",
    list(orgdb = "org.Rn.eg.db", kegg_code = "rno")
  )
  expect_equal(rat_override$orgdb, "org.Rn.eg.db")
  expect_equal(rat_override$kegg, "rno")
  expect_equal(species_resources("Drosophila melanogaster")$species, "drosophila_melanogaster")
  expect_equal(species_resources("homo_sapiens")$species, "human")
  expect_equal(species_resources("Homo-Sapiens")$species, "human")
  expect_equal(species_resources("mmu")$species, "mouse")
  expect_equal(scRDSreport:::.fa_species("rattus-norvegicus"), "rat")
  expect_equal(scRDSreport:::.fa_species("danio_rerio"), "zebrafish")
})

test_that("explicit species provenance is not confused with automatic detection", {
  detected <- list(species = "mouse", confidence = "high", basis = "ENSMUSG feature IDs")

  explicit <- scRDSreport:::.species_selection_status("human", detected)
  expect_equal(explicit$selected, "human")
  expect_equal(explicit$detected, "mouse")
  expect_equal(explicit$confidence, "user")
  expect_equal(explicit$basis, "running species argument")

  configured <- scRDSreport:::.species_selection_status("auto", detected, "rat")
  expect_equal(configured$selected, "rat")
  expect_equal(configured$confidence, "user")
  expect_equal(configured$basis, "config resource_overrides$species")

  automatic <- scRDSreport:::.species_selection_status("auto", detected)
  expect_equal(automatic$selected, "mouse")
  expect_equal(automatic$confidence, "high")
  expect_equal(automatic$basis, "ENSMUSG feature IDs")
})

test_that("raw-barcode QC applies count thresholds before sparse mitochondrial filtering", {
  skip_if_not_installed("SeuratObject")
  counts <- matrix(
    0,
    nrow = 3L,
    ncol = 12L,
    dimnames = list(paste0("ENSMUSG", 1:3), paste0("cell", 1:12))
  )
  counts[2:3, ] <- 10
  counts[1, 1] <- 100
  object <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE)
  )
  feature_metadata <- data.frame(
    gene_symbols = c("mt-Nd1", "Actb", "Gapdh"),
    row.names = rownames(object)
  )
  object[["RNA"]] <- SeuratObject::AddMetaData(object[["RNA"]], feature_metadata)

  mitochondrial <- scRDSreport:::.mitochondrial_percent(object)
  expect_equal(mitochondrial$features, 1L)
  expect_gt(mitochondrial$percent[[1L]], 80)
  expect_equal(mitochondrial$percent[-1L], rep(0, 11L))

  filtered <- scRDSreport:::.prefilter_raw_barcodes(
    object,
    mode = "always",
    min_features = 1,
    min_counts = 1,
    max_percent_mt = 50,
    verbose = FALSE
  )
  expect_true(filtered$summary$applied)
  expect_equal(filtered$summary$mitochondrial_features, 1L)
  expect_equal(filtered$summary$cells_removed_mitochondrial, 1L)
  expect_equal(filtered$summary$cells_after, 11L)
  expect_false("cell1" %in% colnames(filtered$object))
  expect_s4_class(scRDSreport:::.layer_data(filtered$object, "RNA", "counts"), "sparseMatrix")
})

test_that("partial analyzed objects are not automatically barcode-filtered", {
  object <- .release_v020_object(n_cells = 20L)
  result <- scRDSreport:::.prefilter_raw_barcodes(
    object,
    mode = "auto",
    min_features = 10L,
    min_counts = 10L,
    max_percent_mt = 20,
    verbose = FALSE
  )
  expect_false(result$summary$applied)
  expect_equal(result$summary$cells_before, 20L)
  expect_equal(result$summary$cells_after, 20L)
  expect_match(result$summary$reason, "does not look like", ignore.case = TRUE)
})

test_that("unreviewed automatic sample groups never trigger inferential DE", {
  object <- .release_v020_object(n_features = 60L, n_cells = 40L)
  object$group <- NULL
  sample_ids <- rep(c("WT_1", "WT_2", "KO_1", "KO_2"), each = 10L)
  inferred <- infer_sample_design(sample_ids)
  expect_true(all(inferred$design$needs_review))
  object <- scRDSreport:::.attach_design(object, inferred)
  metadata <- scRDSreport:::.seurat_metadata(object)
  expect_true(all(metadata$.scRDSreport_design_needs_review))

  output <- tempfile("scrdsreport-unreviewed-de-")
  result <- scRDSreport:::.fa_module_differential(
    object = object,
    output = output,
    cfg = list(strategy = "auto", min_cells_per_stratum = 10L),
    seed = 1L,
    verbose = FALSE
  )

  expect_equal(result$status, "needs_input")
  expect_equal(result$reason_code, "sample_design_needs_review")
  expect_false(result$details$inferential)
  expect_true(result$details$design_needs_review)
  rankings <- scRDSreport:::.fa_get_analysis_misc(result$object, "differential_rankings")
  expect_false("PValue" %in% names(rankings))
  expect_false("FDR" %in% names(rankings))
  result_artifacts <- Filter(
    function(x) identical(x$type, "table") && grepl("de_", basename(x$path), fixed = TRUE),
    result$artifacts
  )
  expect_gt(length(result_artifacts), 0L)
  for (artifact in result_artifacts) {
    table <- utils::read.csv(file.path(output, artifact$path), check.names = FALSE)
    expect_false("PValue" %in% names(table) && any(is.finite(table$PValue)))
    expect_false("FDR" %in% names(table) && any(is.finite(table$FDR)))
  }
})

test_that("full orchestrator returns twelve statuses and typed, complete artifacts", {
  object <- .release_v020_object()
  output <- tempfile("scrdsreport-full-core-")
  dir.create(output)
  config <- report_config(
    profile = "report_only",
    modules = c("qc", "reduction", "cluster"),
    module_options = list(
      qc = list(
        min_features = 0, min_counts = 0, max_features = Inf,
        max_counts = Inf, max_percent_mt = Inf, filter = FALSE,
        max_plot_cells = 100L
      ),
      cluster = list(run_markers = FALSE)
    )
  )

  result <- scRDSreport:::.run_full_analysis(
    object = object,
    output = output,
    species_info = list(selected = "human"),
    config = config,
    seed = 17L,
    verbose = FALSE
  )

  expect_equal(result$schema_version, "2.0")
  expect_identical(names(result$modules), scRDSreport:::.report_module_ids())
  expect_equal(
    vapply(result$modules[c("qc", "reduction", "cluster", "downloads")], `[[`, character(1), "status"),
    c(qc = "completed", reduction = "completed", cluster = "completed", downloads = "completed")
  )
  expect_true(all(vapply(
    result$modules[setdiff(scRDSreport:::.report_module_ids(), c("qc", "reduction", "cluster", "downloads"))],
    function(x) identical(x$status, "skipped") && identical(x$reason_code, "disabled"),
    logical(1)
  )))

  required <- c(
    "artifact_id", "module", "type", "format", "path", "label",
    "description", "rows", "columns", "row_unit", "column_unit",
    "column_dictionary", "bytes", "sha256", "complete"
  )
  expect_gt(length(result$artifacts), 8L)
  expect_true(all(vapply(result$artifacts, function(x) all(required %in% names(x)), logical(1))))
  expect_true(all(vapply(result$artifacts, function(x) isTRUE(x$complete), logical(1))))
  expect_true(all(vapply(
    result$artifacts,
    function(x) file.exists(file.path(output, x$path)),
    logical(1)
  )))
  hashes <- vapply(result$artifacts, `[[`, character(1), "sha256")
  hashes <- hashes[!is.na(hashes)]
  if (length(hashes)) expect_true(all(grepl("^[0-9a-f]{64}$", hashes)))

  artifact_index <- result$artifacts[vapply(
    result$artifacts,
    function(x) identical(x$module, "downloads") && grepl("artifact_index", x$path),
    logical(1)
  )]
  expect_length(artifact_index, 1L)
})

test_that("trajectory_root public schema resolves explicit cells and metadata groups", {
  object <- .release_v020_object(n_cells = 20L)
  object$state <- rep(c("early", "late"), each = 10L)
  metadata <- scRDSreport:::.seurat_metadata(object)
  built <- list(cells = colnames(object), features = rownames(object))

  group_config <- report_config(
    trajectory_root = list(type = "metadata_group", column = "state", values = "early")
  )
  group_root <- scRDSreport:::.fa_resolve_trajectory_root(
    object,
    built,
    metadata,
    scRDSreport:::.fa_module_config(group_config, "pseudotime")
  )
  expect_true(group_root$supplied)
  expect_true(group_root$valid)
  expect_equal(group_root$type, "metadata_group")
  expect_setequal(group_root$root_cells, rownames(metadata)[metadata$state == "early"])

  cell_config <- report_config(
    trajectory_root = list(type = "cells", cells = c("cell2", "cell4", "absent"))
  )
  cell_root <- scRDSreport:::.fa_resolve_trajectory_root(
    object,
    built,
    metadata,
    scRDSreport:::.fa_module_config(cell_config, "pseudotime")
  )
  expect_true(cell_root$valid)
  expect_equal(cell_root$type, "cells")
  expect_setequal(cell_root$root_cells, c("cell2", "cell4"))

  node_config <- report_config(trajectory_root = "Y_7")
  node_root <- scRDSreport:::.fa_resolve_trajectory_root(
    object,
    built,
    metadata,
    scRDSreport:::.fa_module_config(node_config, "pseudotime")
  )
  expect_true(node_root$valid)
  expect_equal(node_root$root_pr_nodes, "Y_7")
})

test_that("missing TxDb resources stop CNV before any dense conversion", {
  object <- .release_v020_object(n_features = 60L, n_cells = 20L)
  input_counts <- scRDSreport:::.fa_matrix(object, "counts")
  expect_s4_class(input_counts, "sparseMatrix")
  input_dimensions <- dim(input_counts)

  missing_txdb <- "scRDSreport.DefinitelyMissing.TxDb"
  expect_null(scRDSreport:::.fa_txdb_object(missing_txdb))
  expect_null(scRDSreport:::.fa_gene_order(list(txdb = missing_txdb, orgdb = NULL)))

  result <- scRDSreport:::.fa_module_cnv(
    object = object,
    output = tempfile("scrdsreport-cnv-resource-"),
    cfg = list(
      reference_groups = "TypeA",
      annotation_column = "celltype_manual",
      txdb = missing_txdb,
      min_ordered_genes = 10L
    ),
    seed = 1L,
    verbose = FALSE,
    species = "unknown"
  )
  expect_true(result$status %in% c("skipped", "needs_input"))
  expect_true(result$reason_code %in% c("dependency_missing", "gene_order_missing"))
  output_counts <- scRDSreport:::.fa_matrix(result$object, "counts")
  expect_s4_class(output_counts, "sparseMatrix")
  expect_identical(dim(output_counts), input_dimensions)
})

test_that("Quarto template retains all modules, needs-input states, embedded downloads, and dictionaries", {
  template <- system.file("quarto", "report.qmd", package = "scRDSreport")
  expect_true(nzchar(template))
  text <- paste(readLines(template, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

  chapters <- c(
    "质量控制（QC）", "降维分析", "聚类分析", "细胞类型注释",
    "差异表达分析", "富集分析", "伪时间分析", "细胞通讯分析",
    "细胞周期分析", "转录因子分析", "拷贝数变异（CNV）分析",
    "完整结果与下载"
  )
  headings <- sub("^# +", "", grep("^# ", strsplit(text, "\n", fixed = TRUE)[[1L]], value = TRUE))
  expect_true(all(chapters %in% headings))
  expect_match(text, 'return\\("needs_input"\\)')
  expect_match(text, "needs_input = \"需要输入\"", fixed = TRUE)
  expect_match(text, "window.scRDSDownload=function", fixed = TRUE)
  expect_match(text, "script[data-scrds-file", fixed = TRUE)
  expect_match(text, "下载内置文件", fixed = TRUE)
  expect_match(text, "render_table_block <- function", fixed = TRUE)
  expect_match(text, "interpretation_panel(", fixed = TRUE)
  expect_match(text, "artifact_column_dictionary <- function", fixed = TRUE)
  expect_match(text, "展开查看每一列的含义", fixed = TRUE)
})

test_that("legacy and full-analysis artifacts merge with relative paths and SHA-256", {
  output <- tempfile("scrdsreport-manifest-v2-")
  dir.create(file.path(output, "tables"), recursive = TRUE)
  legacy_path <- file.path(output, "tables", "legacy.csv")
  writeLines(c("cell,value", "c1,1"), legacy_path)
  legacy <- scRDSreport:::.manifest_row(
    section = "metadata_analysis",
    label = "Legacy metadata",
    path = legacy_path,
    root = output,
    rows = 1L,
    columns = 2L
  )

  artifact <- scRDSreport:::.fa_write_table_artifact(
    data.frame(cell = c("c1", "c2"), score = c(1, 2)),
    output = output,
    module = "qc",
    name = "synthetic_qc",
    label = "Synthetic QC",
    description = "One row per cell.",
    row_unit = "cell",
    column_dictionary = list(score = "Synthetic test score.")
  )
  module_rows <- scRDSreport:::.artifacts_to_manifest_rows(list(artifact), output)
  merged <- scRDSreport:::.bind_rows_fill(legacy, module_rows)

  expect_equal(nrow(merged), 2L)
  expect_equal(legacy$path, "tables/legacy.csv")
  expect_false(any(grepl(paste0("^", output), merged$path)))
  expect_equal(module_rows$module, "qc")
  expect_equal(module_rows$rows, 2L)
  expect_match(module_rows$column_dictionary, "Synthetic test score", fixed = TRUE)
  expect_true(file.exists(file.path(output, module_rows$path)))
  expect_match(module_rows$sha256, "^[0-9a-f]{64}$")
  expect_match(legacy$sha256, "^[0-9a-f]{64}$")
})
