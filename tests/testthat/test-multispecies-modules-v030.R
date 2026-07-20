test_that("explicit gene sets take precedence and retain provenance", {
  resources <- scRDSreport:::.fa_species_resources("mouse", list())
  sets <- scRDSreport:::.fa_gene_sets(
    list(
      gene_sets = list(pathway_a = c("Actb", "Gapdh", "Myc")),
      auto_gene_sets = TRUE
    ),
    species = "mouse",
    resources = resources
  )
  expect_equal(unname(sets$pathway_a), c("Actb", "Gapdh", "Myc"))
  expect_equal(attr(sets, "source"), "user_gene_sets_or_gmt")
  expect_equal(attr(sets, "msigdbr_species"), "Mus musculus")
  expect_true(is.na(attr(sets, "msigdbr_db_species")))
  expect_equal(attr(sets, "collections"), "user_supplied")

  disabled <- scRDSreport:::.fa_gene_sets(
    list(auto_gene_sets = FALSE), species = "rat",
    resources = scRDSreport:::.fa_species_resources("rat", list())
  )
  expect_length(disabled, 0L)
  expect_equal(attr(disabled, "source"), "none")
  expect_equal(attr(disabled, "collections"), "H")

  human_disabled <- scRDSreport:::.fa_gene_sets(
    list(auto_gene_sets = FALSE), species = "human",
    resources = scRDSreport:::.fa_species_resources("human", list())
  )
  expect_equal(attr(human_disabled, "collections"), "H")
})

test_that("cell-cycle resources record native or ortholog mapping", {
  human <- scRDSreport:::.fa_cell_cycle_genes("human", list())
  expect_true(length(human$s.genes) > 20L)
  expect_true(length(human$g2m.genes) > 20L)
  expect_true(all(c("phase_set", "source_human_symbol", "target_symbol") %in% names(human$mapping)))
  expect_match(human$source, "native human")

  custom <- scRDSreport:::.fa_cell_cycle_genes(
    "zebrafish", list(s_genes = c("mcm5", "pcna"), g2m_genes = c("cdk1", "top2a"))
  )
  expect_equal(custom$source, "user")
  expect_equal(custom$s.genes, c("mcm5", "pcna"))
})

test_that("auto_if_missing resolves to preserve when an RDS annotation exists", {
  counts <- Matrix::Matrix(
    matrix(
      sample.int(5L, 240L, replace = TRUE) - 1L,
      nrow = 24L,
      dimnames = list(paste0("gene", seq_len(24L)), paste0("cell", seq_len(10L)))
    ),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts)
  object$celltype <- rep(c("T", "B"), each = 5L)
  object$sample <- rep(c("s1", "s2"), each = 5L)
  output <- tempfile("celltype-policy-")
  dir.create(output)
  result <- scRDSreport:::.fa_module_celltype(
    object, output, list(mode = "auto_if_missing"), seed = 11L,
    verbose = FALSE, species = "mouse"
  )
  expect_equal(result$status, "completed")
  expect_equal(result$reason_code, "annotation_preserved")
  expect_equal(result$details$requested_mode, "auto_if_missing")
  expect_equal(result$details$mode, "preserve")
  expect_equal(result$details$annotation_column, "celltype")
  expect_false(".scRDSreport_celltype_SingleR" %in% names(result$object[[]]))
})

test_that("invalid explicit annotation references never fall back silently", {
  expect_null(scRDSreport:::.fa_load_annotation_reference(
    "mouse",
    list(reference = "missing-reference.rds", allow_reference_download = TRUE)
  ))
  expect_null(scRDSreport:::.fa_load_annotation_reference(
    "mouse",
    list(reference = "celldex::TypoData", allow_reference_download = TRUE)
  ))
})

test_that("annotation reference loader preserves a loader failure", {
  testthat::local_mocked_bindings(
    .fa_pkg_fun = function(package, name, exported = TRUE) {
      function(...) stop("Cannot open lock file", call. = FALSE)
    },
    .package = "scRDSreport"
  )
  loaded <- scRDSreport:::.fa_load_annotation_reference(
    "mouse", list(reference = "celldex::ImmGenData")
  )
  expect_s3_class(loaded, "scRDSreport_reference_load_error")
  expect_equal(loaded$reference, "celldex::ImmGenData")
  expect_match(loaded$message, "lock file")
})

test_that("Entrez feature IDs map through the species OrgDb", {
  skip_if_not_installed("org.Mm.eg.db")
  orgdb <- scRDSreport:::.fa_orgdb(
    scRDSreport:::.fa_species_resources("mouse", list()), list()
  )
  features <- c("12566", "12567")
  mapping <- scRDSreport:::.fa_feature_mapping(
    features, stats::setNames(features, features), orgdb,
    preferred_keytype = "ENTREZID"
  )
  expect_equal(mapping$SYMBOL, c("Cdk2", "Cdk4"))
  expect_equal(mapping$ENTREZID, features)
  expect_equal(unique(mapping$mapping_keytype), "ENTREZID")
})

test_that("module resource overrides reach the runtime resource object", {
  config <- report_config(resource_overrides = list(
    orgdb = "org.Custom.eg.db",
    scientific_name = "Custom species",
    msigdbr_species = "Custom species",
    msigdbr_db_species = "HS"
  ))
  cfg <- scRDSreport:::.fa_module_config(config, "enrichment")
  resources <- scRDSreport:::.fa_species_resources("mouse", cfg)
  expect_equal(resources$orgdb, "org.Custom.eg.db")
  expect_equal(resources$scientific_name, "Custom species")
  expect_equal(resources$msigdbr_species, "Custom species")
  expect_equal(resources$msigdbr_db_species, "HS")
  expect_equal(resources$msigdbr_default_collection, "H")
  expect_true(resources$msigdbr_ortholog_projection)
})
