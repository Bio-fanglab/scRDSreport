test_that("common species registry exposes complete core metadata", {
  expected <- data.frame(
    species = c("human", "mouse", "rat", "zebrafish", "pig", "cattle", "chicken", "dog", "macaque"),
    scientific_name = c(
      "Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio",
      "Sus scrofa", "Bos taurus", "Gallus gallus", "Canis lupus familiaris",
      "Macaca mulatta"
    ),
    orgdb = c(
      "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Dr.eg.db",
      "org.Ss.eg.db", "org.Bt.eg.db", "org.Gg.eg.db", "org.Cf.eg.db",
      "org.Mmu.eg.db"
    ),
    kegg = c("hsa", "mmu", "rno", "dre", "ssc", "bta", "gga", "cfa", "mcc"),
    stringsAsFactors = FALSE
  )

  observed <- scRDSreport:::supported_species()
  expect_s3_class(observed, "data.frame")
  expect_equal(observed$species, expected$species)
  expect_equal(observed$scientific_name, expected$scientific_name)
  expect_equal(observed$orgdb, expected$orgdb)
  expect_equal(observed$kegg_code, expected$kegg)
  expect_true(all(!is.na(observed$taxonomy_id)))
  expect_true(all(!is.na(observed$ensembl_prefix)))
  expect_true(all(!is.na(observed$msigdbr_species)))
  expect_equal(observed$msigdbr_db_species, c("HS", "MM", rep("HS", 7L)))
  expect_equal(observed$msigdbr_default_collection, c("H", "MH", rep("H", 7L)))
  expect_equal(observed$msigdbr_ortholog_projection, c(FALSE, FALSE, rep(TRUE, 7L)))
  expect_equal(observed$capability_tier, c("full", "full", rep("core", 7L)))

  mouse <- species_resources("mouse")
  expect_equal(mouse$auto_annotation_reference, "celldex::MouseRNAseqData")
  expect_equal(mouse$msigdbr_default_collection, "MH")
  expect_true("celldex::ImmGenData" %in% mouse$reference_candidates)
  expect_match(mouse$annotation_context, "immune-focused")
})

test_that("scientific names and common aliases normalize to registry keys", {
  aliases <- c(
    "Homo sapiens" = "human",
    "Mus musculus" = "mouse",
    "Rattus norvegicus" = "rat",
    "Danio rerio" = "zebrafish",
    "Sus scrofa" = "pig",
    "Bos taurus" = "cattle",
    "Gallus gallus" = "chicken",
    "Canis lupus familiaris" = "dog",
    "Macaca mulatta" = "macaque",
    "swine" = "pig",
    "cow" = "cattle",
    "canine" = "dog",
    "rhesus macaque" = "macaque"
  )
  observed <- vapply(names(aliases), function(value) {
    species_resources(value)$species
  }, character(1))
  expect_equal(unname(observed), unname(aliases))
})

test_that("core profiles never borrow human or mouse CellChat and TxDb resources", {
  for (species in c("rat", "zebrafish", "pig", "cattle", "chicken", "dog", "macaque")) {
    resources <- species_resources(species)
    expect_true(resources$supported, info = species)
    expect_equal(resources$capability_tier, "core", info = species)
    expect_null(resources$cellchat_db, info = species)
    expect_null(resources$txdb, info = species)
    expect_null(resources$auto_annotation_reference, info = species)
    expect_false(resources$capabilities[["cellchat"]], info = species)
    expect_false(resources$capabilities[["cnv_coordinates"]], info = species)
    expect_true(resources$capabilities[["identifier_mapping"]], info = species)
    expect_true(resources$capabilities[["go_kegg"]], info = species)
    expect_true(resources$capabilities[["msigdb"]], info = species)
    expect_true(resources$capabilities[["cell_cycle"]], info = species)
    expect_true(resources$capabilities[["tf_catalog"]], info = species)
    expect_equal(resources$cell_cycle_strategy, "babelgene_human_seurat_cc_orthologs", info = species)
    expect_equal(resources$tf_catalog_strategy, "orgdb_go_tf_catalog", info = species)
  }

  rat <- species_resources("rat", overrides = list(
    cellchat_db = "RatCellChatDB.explicit",
    gtf = "mRatBN7.2.explicit.gtf"
  ))
  expect_equal(rat$source, "built_in_with_overrides")
  expect_equal(rat$cellchat_db, "RatCellChatDB.explicit")
  expect_true(rat$capabilities[["cellchat"]])
  expect_true(rat$capabilities[["cnv_coordinates"]])
})

test_that("flat and nested QC patterns stay synchronized", {
  zebrafish <- species_resources("zebrafish")
  expect_equal(zebrafish$patterns$mitochondrial, zebrafish$mitochondrial_pattern)
  expect_equal(zebrafish$patterns$ribosomal, zebrafish$ribosomal_pattern)
  expect_equal(zebrafish$patterns$hemoglobin, zebrafish$hemoglobin_pattern)

  custom <- species_resources("rat", overrides = list(
    patterns = list(mitochondrial = "^mito_", ignore_case = FALSE)
  ))
  expect_equal(custom$mitochondrial_pattern, "^mito_")
  expect_equal(custom$patterns$mitochondrial, "^mito_")
  expect_false(custom$pattern_ignore_case)
})

test_that("stable Ensembl prefixes detect all registered common species", {
  prefixes <- c(
    human = "ENSG", mouse = "ENSMUSG", rat = "ENSRNOG",
    zebrafish = "ENSDARG", pig = "ENSSSCG", cattle = "ENSBTAG",
    chicken = "ENSGALG", dog = "ENSCAFG", macaque = "ENSMMUG"
  )
  for (species in names(prefixes)) {
    features <- paste0(prefixes[[species]], sprintf("%011d", seq_len(30L)), ".1")
    result <- scRDSreport:::.detect_species(features)
    expect_equal(result$species, species, info = species)
    expect_equal(result$confidence, "high", info = species)
    expect_match(result$basis, prefixes[[species]], fixed = TRUE, info = species)
  }
})

test_that("automatic detection refuses mixed-species and empty feature sets", {
  mixed <- c(
    paste0("ENSG", sprintf("%011d", seq_len(30L))),
    paste0("ENSMUSG", sprintf("%011d", seq_len(30L)))
  )
  result <- scRDSreport:::.detect_species(mixed)
  expect_equal(result$species, "unknown")
  expect_match(result$basis, "mixed or ambiguous")

  empty <- scRDSreport:::.detect_species(character())
  expect_equal(empty$species, "unknown")
  expect_equal(empty$basis, "no feature IDs")

  uppercase <- scRDSreport:::.detect_species(c("TP53", "GAPDH", "ACTB", "CD3D"))
  expect_equal(uppercase$species, "unknown")
  expect_match(uppercase$basis, "cannot distinguish")
})

test_that("MSigDB override derivations stay internally consistent", {
  mouse_human_db <- species_resources(
    "mouse", overrides = list(msigdbr_db_species = "HS")
  )
  expect_equal(mouse_human_db$msigdbr_default_collection, "H")
  expect_true(mouse_human_db$msigdbr_ortholog_projection)
  expect_match(mouse_human_db$gene_sets_strategy, "msigdbr_hs_to_")
})
