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
    automatic = "auto", unspecified = "unknown"
  )
  if (value %in% names(aliases)) value <- unname(aliases[[value]])
  value
}

.species_selection_status <- function(requested, detected, configured = NULL) {
  requested <- .normalize_species_name(requested)
  if (!identical(requested, "auto")) {
    return(list(
      requested = requested,
      detected = detected$species,
      selected = requested,
      confidence = "user",
      basis = "running species argument"
    ))
  }
  if (!is.null(configured)) {
    return(list(
      requested = requested,
      detected = detected$species,
      selected = .normalize_species_name(configured),
      confidence = "user",
      basis = "config resource_overrides$species"
    ))
  }
  list(
    requested = requested,
    detected = detected$species,
    selected = detected$species,
    confidence = detected$confidence,
    basis = detected$basis
  )
}

.empty_species_resources <- function(species) {
  list(
    schema_version = "1.0",
    species = species,
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
    chromosomes = character(),
    cell_cycle_strategy = "user_supplied",
    tf_catalog_strategy = "user_supplied",
    gene_sets_strategy = "user_supplied",
    auto_annotation_reference = NULL,
    gtf = NULL,
    gene_sets = NULL,
    manual_markers = NULL
  )
}

.human_species_resources <- function() {
  list(
    schema_version = "1.0",
    species = "human",
    supported = TRUE,
    source = "built_in",
    reason = "Built-in human resource mapping.",
    orgdb = "org.Hs.eg.db",
    feature_keytype = "ENSEMBL",
    symbol_column = "SYMBOL",
    kegg_code = "hsa",
    cellchat_db = "CellChatDB.human",
    txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    genome_assembly = "GRCh38/hg38",
    ensembl_prefix = "^ENSG[0-9]",
    mitochondrial_pattern = "^MT-",
    ribosomal_pattern = "^RP[SL]",
    hemoglobin_pattern = "^HB[ABDEGQMZ]",
    pattern_ignore_case = FALSE,
    chromosomes = c(as.character(seq_len(22L)), "X", "Y", "MT"),
    cell_cycle_strategy = "seurat_cc_genes_human_symbols",
    tf_catalog_strategy = "human_curated_or_user_supplied",
    gene_sets_strategy = "msigdbr_homo_sapiens",
    auto_annotation_reference = "celldex::HumanPrimaryCellAtlasData",
    gtf = NULL,
    gene_sets = NULL,
    manual_markers = NULL
  )
}

.mouse_species_resources <- function() {
  list(
    schema_version = "1.0",
    species = "mouse",
    supported = TRUE,
    source = "built_in",
    reason = "Built-in mouse resource mapping.",
    orgdb = "org.Mm.eg.db",
    feature_keytype = "ENSEMBL",
    symbol_column = "SYMBOL",
    kegg_code = "mmu",
    cellchat_db = "CellChatDB.mouse",
    txdb = "TxDb.Mmusculus.UCSC.mm10.knownGene",
    genome_assembly = "GRCm38/mm10",
    ensembl_prefix = "^ENSMUSG[0-9]",
    mitochondrial_pattern = "^mt-",
    ribosomal_pattern = "^Rp[sl]",
    hemoglobin_pattern = "^Hb[ab]",
    pattern_ignore_case = FALSE,
    chromosomes = c(as.character(seq_len(19L)), "X", "Y", "MT"),
    cell_cycle_strategy = "seurat_cc_genes_titlecase_mouse_symbols",
    tf_catalog_strategy = "mouse_curated_or_user_supplied",
    gene_sets_strategy = "msigdbr_mus_musculus",
    auto_annotation_reference = "celldex::ImmGenData",
    gtf = NULL,
    gene_sets = NULL,
    manual_markers = NULL
  )
}

.has_custom_species_resources <- function(resources) {
  fields <- c(
    "orgdb", "kegg_code", "cellchat_db", "txdb", "gtf", "gene_sets",
    "manual_markers", "ensembl_prefix"
  )
  any(vapply(fields, function(name) {
    value <- resources[[name]]
    !is.null(value) && length(value) > 0L
  }, logical(1)))
}

#' Resolve organism-specific resources without implicit mouse fallback
#'
#' Built-in mappings are provided for human and mouse. `"auto"`, `"unknown"`,
#' and unregistered species return an unsupported resource object containing no
#' mouse database, marker, chromosome, or genome fallback. Other species can be
#' enabled only through explicit `overrides`.
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

  resources <- switch(
    species,
    human = .human_species_resources(),
    mouse = .mouse_species_resources(),
    .empty_species_resources(species)
  )
  resources <- .merge_resource_lists(resources, overrides)
  resources$species <- species

  if (!identical(species, "human") && !identical(species, "mouse") &&
      .has_custom_species_resources(resources)) {
    resources$supported <- TRUE
    resources$source <- "user_override"
    resources$reason <- paste0("Organism-specific resources for '", species, "' were supplied explicitly.")
  }

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
