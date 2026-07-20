.merge_resource_lists <- function(base, overrides) {
  if (!length(overrides)) return(base)
  utils::modifyList(base, overrides, keep.null = TRUE)
}

.normalize_species_name <- function(species) {
  if (!is.character(species) || length(species) != 1L || is.na(species) ||
      !nzchar(trimws(species))) {
    stop("species must be one non-empty character value.", call. = FALSE)
  }
  value <- tolower(trimws(species))
  value <- gsub("[ -]+", "_", value)
  aliases <- c(
    homo_sapiens = "human", hsapiens = "human", hs = "human", hsa = "human",
    mus_musculus = "mouse", mmusculus = "mouse", mm = "mouse", mmu = "mouse",
    rattus_norvegicus = "rat", rnorvegicus = "rat", rn = "rat", rno = "rat",
    danio_rerio = "zebrafish", drerio = "zebrafish", dr = "zebrafish", dre = "zebrafish",
    zebra_fish = "zebrafish",
    sus_scrofa = "pig", sscrofa = "pig", ssc = "pig", swine = "pig", porcine = "pig",
    bos_taurus = "cattle", btaurus = "cattle", bta = "cattle", cow = "cattle", bovine = "cattle",
    gallus_gallus = "chicken", ggallus = "chicken", gga = "chicken",
    canis_lupus_familiaris = "dog", canis_familiaris = "dog", cfa = "dog", canine = "dog",
    macaca_mulatta = "macaque", mmulatta = "macaque", mcc = "macaque",
    rhesus = "macaque", rhesus_macaque = "macaque", rhesus_monkey = "macaque",
    automatic = "auto", unspecified = "unknown"
  )
  if (value %in% names(aliases)) value <- unname(aliases[[value]])
  value
}

.species_selection_status <- function(requested, detected, configured = NULL) {
  requested <- .normalize_species_name(requested)
  conflict_for <- function(selected) {
    detected_species <- .normalize_species_name(detected$species %||% "unknown")
    known <- !detected_species %in% c("auto", "unknown")
    known && !identical(.normalize_species_name(selected), detected_species)
  }
  if (!identical(requested, "auto")) {
    conflict <- conflict_for(requested)
    return(list(
      requested = requested,
      detected = detected$species,
      selected = requested,
      confidence = "user",
      basis = "running species argument",
      conflict = conflict,
      message = if (conflict) {
        paste0("User-selected species '", requested, "' conflicts with feature-ID detection '", detected$species, "'. The explicit user choice controls resources; verify the input.")
      } else {
        "The explicit running() species argument controls organism-specific resources."
      }
    ))
  }
  if (!is.null(configured)) {
    configured <- .normalize_species_name(configured)
    conflict <- conflict_for(configured)
    return(list(
      requested = requested,
      detected = detected$species,
      selected = configured,
      confidence = "user",
      basis = "config resource_overrides$species",
      conflict = conflict,
      message = if (conflict) {
        paste0("Configured species '", configured, "' conflicts with feature-ID detection '", detected$species, "'. Verify the input and override.")
      } else {
        "The configured species override controls organism-specific resources."
      }
    ))
  }
  list(
    requested = requested,
    detected = detected$species,
    selected = detected$species,
    confidence = detected$confidence,
    basis = detected$basis,
    conflict = FALSE,
    message = paste0("Automatic species selection used: ", detected$basis, ".")
  )
}

.resource_capabilities <- function(resources) {
  populated <- function(name) {
    value <- resources[[name]]
    !is.null(value) && length(value) > 0L && !all(is.na(value))
  }
  c(
    identifier_mapping = populated("orgdb"),
    qc_patterns = all(vapply(
      c("mitochondrial_pattern", "ribosomal_pattern", "hemoglobin_pattern"),
      populated, logical(1)
    )),
    go_kegg = populated("orgdb") && populated("kegg_code"),
    msigdb = populated("msigdbr_species"),
    automatic_annotation = populated("auto_annotation_reference"),
    cellchat = populated("cellchat_db"),
    cnv_coordinates = populated("txdb") || populated("gtf"),
    cell_cycle = populated("cell_cycle_strategy") &&
      !identical(resources$cell_cycle_strategy, "user_supplied"),
    tf_catalog = populated("tf_catalog_strategy") &&
      !identical(resources$tf_catalog_strategy, "user_supplied")
  )
}

.species_profile <- function(species, scientific_name, taxonomy_id,
                             orgdb, kegg_code, ensembl_prefix,
                             genome_assembly, chromosomes, msigdbr_species,
                             msigdbr_db_species, msigdbr_default_collection,
                             mitochondrial_pattern, ribosomal_pattern,
                             hemoglobin_pattern, pattern_ignore_case = TRUE,
                             capability_tier = "core", cellchat_db = NULL,
                             txdb = NULL, auto_annotation_reference = NULL,
                             annotation_context = NULL,
                             reference_candidates = character(),
                             cell_cycle_strategy = "babelgene_human_seurat_cc_orthologs",
                             tf_catalog_strategy = "orgdb_go_tf_catalog") {
  resources <- list(
    schema_version = "1.1",
    species = species,
    scientific_name = scientific_name,
    taxonomy_id = as.integer(taxonomy_id),
    capability_tier = capability_tier,
    supported = TRUE,
    source = "built_in",
    reason = paste0(
      "Built-in ", scientific_name, " resource mapping (", capability_tier,
      " capability tier)."
    ),
    orgdb = orgdb,
    feature_keytype = "ENSEMBL",
    symbol_column = "SYMBOL",
    kegg_code = kegg_code,
    cellchat_db = cellchat_db,
    txdb = txdb,
    genome_assembly = genome_assembly,
    ensembl_prefix = ensembl_prefix,
    mitochondrial_pattern = mitochondrial_pattern,
    ribosomal_pattern = ribosomal_pattern,
    hemoglobin_pattern = hemoglobin_pattern,
    pattern_ignore_case = isTRUE(pattern_ignore_case),
    patterns = list(
      mitochondrial = mitochondrial_pattern,
      ribosomal = ribosomal_pattern,
      hemoglobin = hemoglobin_pattern,
      ignore_case = isTRUE(pattern_ignore_case)
    ),
    chromosomes = as.character(chromosomes),
    msigdbr_species = msigdbr_species,
    msigdbr_db_species = msigdbr_db_species,
    msigdbr_default_collection = msigdbr_default_collection,
    msigdbr_ortholog_projection = !identical(
      toupper(msigdbr_db_species),
      if (identical(species, "mouse")) "MM" else if (identical(species, "human")) "HS" else ""
    ),
    cell_cycle_strategy = cell_cycle_strategy,
    tf_catalog_strategy = tf_catalog_strategy,
    gene_sets_strategy = paste0(
      "msigdbr_", tolower(msigdbr_db_species), "_to_",
      gsub("[^a-z0-9]+", "_", tolower(scientific_name))
    ),
    auto_annotation_reference = auto_annotation_reference,
    annotation_context = annotation_context,
    reference_candidates = as.character(reference_candidates),
    gtf = NULL,
    gene_sets = NULL,
    manual_markers = NULL
  )
  resources$capabilities <- .resource_capabilities(resources)
  resources
}

.common_species_registry <- function() {
  list(
    human = .species_profile(
      species = "human", scientific_name = "Homo sapiens", taxonomy_id = 9606L,
      capability_tier = "full", orgdb = "org.Hs.eg.db", kegg_code = "hsa",
      ensembl_prefix = "^ENSG[0-9]", genome_assembly = "GRCh38/hg38",
      chromosomes = c(as.character(seq_len(22L)), "X", "Y", "MT"),
      msigdbr_species = "Homo sapiens", msigdbr_db_species = "HS",
      msigdbr_default_collection = "H",
      mitochondrial_pattern = "^MT-",
      ribosomal_pattern = "^RP[SL]", hemoglobin_pattern = "^HB[ABDEGQMZ]",
      pattern_ignore_case = FALSE, cellchat_db = "CellChatDB.human",
      txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
      auto_annotation_reference = "celldex::HumanPrimaryCellAtlasData",
      annotation_context = "general human reference",
      reference_candidates = c(
        "celldex::HumanPrimaryCellAtlasData", "celldex::BlueprintEncodeData",
        "celldex::MonacoImmuneData"
      ),
      cell_cycle_strategy = "seurat_cc_genes_human_symbols",
      tf_catalog_strategy = "orgdb_go_tf_catalog"
    ),
    mouse = .species_profile(
      species = "mouse", scientific_name = "Mus musculus", taxonomy_id = 10090L,
      capability_tier = "full", orgdb = "org.Mm.eg.db", kegg_code = "mmu",
      ensembl_prefix = "^ENSMUSG[0-9]", genome_assembly = "GRCm38/mm10",
      chromosomes = c(as.character(seq_len(19L)), "X", "Y", "MT"),
      msigdbr_species = "Mus musculus", msigdbr_db_species = "MM",
      msigdbr_default_collection = "MH",
      mitochondrial_pattern = "^mt-",
      ribosomal_pattern = "^Rp[sl]", hemoglobin_pattern = "^Hb[ab]",
      pattern_ignore_case = FALSE, cellchat_db = "CellChatDB.mouse",
      txdb = "TxDb.Mmusculus.UCSC.mm10.knownGene",
      auto_annotation_reference = "celldex::MouseRNAseqData",
      annotation_context = "general mouse reference; ImmGenData is an optional immune-focused alternative",
      reference_candidates = c("celldex::MouseRNAseqData", "celldex::ImmGenData"),
      cell_cycle_strategy = "babelgene_human_seurat_cc_orthologs",
      tf_catalog_strategy = "orgdb_go_tf_catalog"
    ),
    rat = .species_profile(
      species = "rat", scientific_name = "Rattus norvegicus", taxonomy_id = 10116L,
      orgdb = "org.Rn.eg.db", kegg_code = "rno", ensembl_prefix = "^ENSRNOG[0-9]",
      genome_assembly = "mRatBN7.2", chromosomes = c(as.character(seq_len(20L)), "X", "Y", "MT"),
      msigdbr_species = "Rattus norvegicus", msigdbr_db_species = "HS",
      msigdbr_default_collection = "H",
      mitochondrial_pattern = "^Mt-",
      ribosomal_pattern = "^Rp[sl]", hemoglobin_pattern = "^Hb[ab]"
    ),
    zebrafish = .species_profile(
      species = "zebrafish", scientific_name = "Danio rerio", taxonomy_id = 7955L,
      orgdb = "org.Dr.eg.db", kegg_code = "dre", ensembl_prefix = "^ENSDARG[0-9]",
      genome_assembly = "GRCz11", chromosomes = c(as.character(seq_len(25L)), "MT"),
      msigdbr_species = "Danio rerio", msigdbr_db_species = "HS",
      msigdbr_default_collection = "H",
      mitochondrial_pattern = "^mt-",
      ribosomal_pattern = "^rp[sl]", hemoglobin_pattern = "^hb"
    ),
    pig = .species_profile(
      species = "pig", scientific_name = "Sus scrofa", taxonomy_id = 9823L,
      orgdb = "org.Ss.eg.db", kegg_code = "ssc", ensembl_prefix = "^ENSSSCG[0-9]",
      genome_assembly = "Sscrofa11.1", chromosomes = c(as.character(seq_len(18L)), "X", "Y", "MT"),
      msigdbr_species = "Sus scrofa", msigdbr_db_species = "HS",
      msigdbr_default_collection = "H",
      mitochondrial_pattern = "^MT-",
      ribosomal_pattern = "^RP[SL]", hemoglobin_pattern = "^HB"
    ),
    cattle = .species_profile(
      species = "cattle", scientific_name = "Bos taurus", taxonomy_id = 9913L,
      orgdb = "org.Bt.eg.db", kegg_code = "bta", ensembl_prefix = "^ENSBTAG[0-9]",
      genome_assembly = "ARS-UCD1.3", chromosomes = c(as.character(seq_len(29L)), "X", "Y", "MT"),
      msigdbr_species = "Bos taurus", msigdbr_db_species = "HS",
      msigdbr_default_collection = "H",
      mitochondrial_pattern = "^MT-",
      ribosomal_pattern = "^RP[SL]", hemoglobin_pattern = "^HB"
    ),
    chicken = .species_profile(
      species = "chicken", scientific_name = "Gallus gallus", taxonomy_id = 9031L,
      orgdb = "org.Gg.eg.db", kegg_code = "gga", ensembl_prefix = "^ENSGALG[0-9]",
      genome_assembly = "GRCg7b", chromosomes = c(as.character(seq_len(39L)), "W", "Z", "MT"),
      msigdbr_species = "Gallus gallus", msigdbr_db_species = "HS",
      msigdbr_default_collection = "H",
      mitochondrial_pattern = "^MT-",
      ribosomal_pattern = "^RP[SL]", hemoglobin_pattern = "^HB"
    ),
    dog = .species_profile(
      species = "dog", scientific_name = "Canis lupus familiaris", taxonomy_id = 9615L,
      orgdb = "org.Cf.eg.db", kegg_code = "cfa", ensembl_prefix = "^ENSCAFG[0-9]",
      genome_assembly = "ROS_Cfam_1.0", chromosomes = c(as.character(seq_len(38L)), "X", "Y", "MT"),
      msigdbr_species = "Canis lupus familiaris", msigdbr_db_species = "HS",
      msigdbr_default_collection = "H",
      mitochondrial_pattern = "^MT-",
      ribosomal_pattern = "^RP[SL]", hemoglobin_pattern = "^HB"
    ),
    macaque = .species_profile(
      species = "macaque", scientific_name = "Macaca mulatta", taxonomy_id = 9544L,
      orgdb = "org.Mmu.eg.db", kegg_code = "mcc", ensembl_prefix = "^ENSMMUG[0-9]",
      genome_assembly = "Mmul_10", chromosomes = c(as.character(seq_len(20L)), "X", "Y", "MT"),
      msigdbr_species = "Macaca mulatta", msigdbr_db_species = "HS",
      msigdbr_default_collection = "H",
      mitochondrial_pattern = "^MT-",
      ribosomal_pattern = "^RP[SL]", hemoglobin_pattern = "^HB"
    )
  )
}

.human_species_resources <- function() .common_species_registry()$human

.mouse_species_resources <- function() .common_species_registry()$mouse

.empty_species_resources <- function(species) {
  resources <- list(
    schema_version = "1.1",
    species = species,
    scientific_name = NULL,
    taxonomy_id = NA_integer_,
    capability_tier = "unregistered",
    supported = FALSE,
    source = "none",
    reason = if (identical(species, "auto")) {
      "Species has not been resolved yet; no organism-specific resource was selected."
    } else {
      paste0("No built-in organism-specific resources are registered for '", species, "'.")
    },
    orgdb = NULL,
    feature_keytype = NULL,
    symbol_column = "SYMBOL",
    kegg_code = NULL,
    cellchat_db = NULL,
    txdb = NULL,
    genome_assembly = NULL,
    ensembl_prefix = NULL,
    mitochondrial_pattern = NULL,
    ribosomal_pattern = NULL,
    hemoglobin_pattern = NULL,
    pattern_ignore_case = TRUE,
    patterns = list(mitochondrial = NULL, ribosomal = NULL, hemoglobin = NULL, ignore_case = TRUE),
    chromosomes = character(),
    msigdbr_species = NULL,
    msigdbr_db_species = NULL,
    msigdbr_default_collection = NULL,
    msigdbr_ortholog_projection = NA,
    cell_cycle_strategy = "user_supplied",
    tf_catalog_strategy = "user_supplied",
    gene_sets_strategy = "user_supplied",
    auto_annotation_reference = NULL,
    annotation_context = NULL,
    reference_candidates = character(),
    gtf = NULL,
    gene_sets = NULL,
    manual_markers = NULL
  )
  resources$capabilities <- .resource_capabilities(resources)
  resources
}

.has_custom_species_resources <- function(resources) {
  fields <- c(
    "orgdb", "kegg_code", "cellchat_db", "txdb", "gtf", "gene_sets",
    "manual_markers", "ensembl_prefix", "msigdbr_species", "msigdbr_db_species",
    "msigdbr_default_collection", "auto_annotation_reference"
  )
  any(vapply(fields, function(name) {
    value <- resources[[name]]
    !is.null(value) && length(value) > 0L
  }, logical(1)))
}

.synchronize_resource_patterns <- function(resources, overrides) {
  nested <- overrides$patterns
  if (is.list(nested)) {
    nested_to_flat <- c(
      mitochondrial = "mitochondrial_pattern",
      ribosomal = "ribosomal_pattern",
      hemoglobin = "hemoglobin_pattern",
      ignore_case = "pattern_ignore_case"
    )
    for (name in names(nested_to_flat)) {
      flat <- nested_to_flat[[name]]
      if (is.null(overrides[[flat]]) && name %in% names(nested)) resources[[flat]] <- nested[[name]]
    }
  }
  resources$pattern_ignore_case <- isTRUE(resources$pattern_ignore_case)
  resources$patterns <- list(
    mitochondrial = resources$mitochondrial_pattern,
    ribosomal = resources$ribosomal_pattern,
    hemoglobin = resources$hemoglobin_pattern,
    ignore_case = resources$pattern_ignore_case
  )
  resources
}

.synchronize_resource_derivations <- function(resources, overrides) {
  database <- resources$msigdbr_db_species
  if (!is.null(database) && length(database)) {
    database <- toupper(as.character(database[[1L]]))
    if (is.null(overrides$msigdbr_default_collection)) {
      resources$msigdbr_default_collection <- if (identical(database, "MM")) "MH" else "H"
    }
    if (is.null(overrides$msigdbr_ortholog_projection)) {
      resources$msigdbr_ortholog_projection <- !(
        (identical(resources$species, "human") && identical(database, "HS")) ||
          (identical(resources$species, "mouse") && identical(database, "MM"))
      )
    }
    if (is.null(overrides$gene_sets_strategy)) {
      target <- resources$scientific_name %||% resources$species
      resources$gene_sets_strategy <- paste0(
        "msigdbr_", tolower(database), "_to_",
        gsub("[^a-z0-9]+", "_", tolower(target))
      )
    }
  }
  resources
}

#' List built-in common-species resource profiles
#'
#' Returns the species profiles that can be selected directly with
#' [species_resources()]. The `full` tier currently identifies profiles with
#' native automatic-annotation, CellChat, and genome-coordinate resources.
#' The `core` tier provides species-matched identifier mapping, QC patterns,
#' KEGG, MSigDB metadata, ortholog-aware cell-cycle mapping, and an OrgDb TF
#' catalogue strategy. Automatic annotation, CellChat, and CNV coordinates
#' remain unavailable until the user supplies a validated resource for that
#' same species and genome build. Registered capability does not imply that an
#' optional R package is installed; the returned table reports both.
#'
#' @return A data frame with one row per built-in species profile and resource
#'   identifiers or availability flags in columns.
#' @export
supported_species <- function() {
  registry <- .common_species_registry()
  scalar <- function(resource, name) {
    value <- resource[[name]]
    if (is.null(value) || !length(value) || is.na(value[[1L]])) NA_character_ else as.character(value[[1L]])
  }
  rows <- lapply(registry, function(resource) {
    data.frame(
      species = resource$species,
      scientific_name = resource$scientific_name,
      taxonomy_id = resource$taxonomy_id,
      capability_tier = resource$capability_tier,
      orgdb = scalar(resource, "orgdb"),
      kegg_code = scalar(resource, "kegg_code"),
      ensembl_prefix = scalar(resource, "ensembl_prefix"),
      genome_assembly = scalar(resource, "genome_assembly"),
      msigdbr_species = scalar(resource, "msigdbr_species"),
      msigdbr_db_species = scalar(resource, "msigdbr_db_species"),
      msigdbr_default_collection = scalar(resource, "msigdbr_default_collection"),
      msigdbr_ortholog_projection = isTRUE(resource$msigdbr_ortholog_projection),
      automatic_annotation = isTRUE(resource$capabilities[["automatic_annotation"]]),
      automatic_annotation_dependencies_installed = if (isTRUE(resource$capabilities[["automatic_annotation"]])) {
        requireNamespace("SingleR", quietly = TRUE) && requireNamespace("celldex", quietly = TRUE)
      } else FALSE,
      cellchat = isTRUE(resource$capabilities[["cellchat"]]),
      cellchat_dependency_installed = if (isTRUE(resource$capabilities[["cellchat"]])) {
        requireNamespace("CellChat", quietly = TRUE)
      } else FALSE,
      cnv_coordinates = isTRUE(resource$capabilities[["cnv_coordinates"]]),
      orgdb_installed = !is.null(resource$orgdb) && requireNamespace(resource$orgdb, quietly = TRUE),
      txdb_installed = !is.null(resource$txdb) && requireNamespace(resource$txdb, quietly = TRUE),
      msigdbr_installed = requireNamespace("msigdbr", quietly = TRUE),
      babelgene_installed = requireNamespace("babelgene", quietly = TRUE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Resolve organism-specific resources without cross-species fallback
#'
#' Built-in mappings are provided for common human, mouse, rat, zebrafish,
#' pig, cattle, chicken, dog, and macaque datasets. Only human and mouse have
#' built-in CellChat and TxDb resources. Other profiles deliberately leave
#' those fields empty instead of borrowing a database or genome build from a
#' different species. `"auto"`, `"unknown"`, and unregistered species return
#' an unsupported resource object unless matching resources are supplied
#' explicitly through `overrides`.
#'
#' @param species Species name or supported alias.
#' @param overrides Named list replacing resource fields such as `orgdb`,
#'   `kegg_code`, `cellchat_db`, `txdb`, `gtf`, or `gene_sets`.
#' @return A list of species resources with class `scRDSreport_species_resources`.
#' @export
species_resources <- function(species = "auto", overrides = list()) {
  species <- .normalize_species_name(species)
  if (!is.list(overrides) ||
      (length(overrides) && (is.null(names(overrides)) || any(!nzchar(names(overrides)))))) {
    stop("overrides must be a named list.", call. = FALSE)
  }

  registry <- .common_species_registry()
  registered <- species %in% names(registry)
  resources <- if (registered) registry[[species]] else .empty_species_resources(species)
  resources <- .merge_resource_lists(resources, overrides)
  resources$species <- species
  resources <- .synchronize_resource_patterns(resources, overrides)
  resources <- .synchronize_resource_derivations(resources, overrides)

  if (!registered && .has_custom_species_resources(resources)) {
    resources$supported <- TRUE
    resources$source <- "user_override"
    resources$capability_tier <- "custom"
    resources$reason <- paste0("Organism-specific resources for '", species, "' were supplied explicitly.")
  } else if (registered && length(overrides)) {
    resources$source <- "built_in_with_overrides"
    resources$reason <- paste0(
      "Built-in ", resources$scientific_name,
      " resource mapping with explicit user overrides."
    )
  }
  resources$capabilities <- .resource_capabilities(resources)

  structure(resources, class = c("scRDSreport_species_resources", "list"))
}

.resolve_species_resources <- function(config, species = "auto") {
  .validate_report_config(config)
  overrides <- config$resource_overrides
  override_species <- overrides$species
  if (!is.null(override_species)) {
    if (identical(species, "auto") || identical(species, "unknown")) {
      species <- override_species
    }
    overrides$species <- NULL
  }
  species_resources(species, overrides = overrides)
}
