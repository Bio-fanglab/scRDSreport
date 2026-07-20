# scRDSreport 0.3.0

- Expanded the built-in species registry from human and mouse to nine common
  species: human, mouse, rat, zebrafish, pig, cattle, chicken, dog, and rhesus
  macaque. `supported_species()` reports registered resources, capability tiers,
  and which optional dependencies are installed in the current R environment.
- Added species-specific aliases, taxonomy, Ensembl prefixes, QC patterns,
  chromosomes, genome assemblies, OrgDb packages, KEGG codes, and explicit
  conflict reporting when a user-selected species disagrees with feature-ID
  detection. Unknown species still never fall back to mouse.
- Made `annotation_mode = "auto_if_missing"` the default: existing RDS
  annotations are preserved; otherwise a trusted species-matched reference is
  attempted. Built-in references currently use HumanPrimaryCellAtlasData for
  human and the general MouseRNAseqData reference for mouse, with ImmGenData as
  an explicit immune-focused alternative. Generated labels use a separate
  column and uncertain SingleR labels remain missing.
- Human and mouse now use their HS/H and MM/MH MSigDB database/collection
  pairs, respectively. The other
  seven built-in nonhuman species use explicit human-to-target ortholog
  projection through msigdbr, with database version and available ortholog
  evidence exported instead of presenting the result as a native database.
- Added nonhuman cell-cycle mapping from Seurat human S/G2M genes through
  babelgene, and species OrgDb-derived TF catalogues from
  `GOALL:GO:0003700`. Exact cell-cycle mappings and resource provenance are
  downloadable.
- Clarified that `species = "mouse"` selects mouse OrgDb, KEGG, QC, MSigDB,
  cell-cycle, TF, CellChat, annotation, and genome-resource routes, while an
  already analyzed RDS still needs an explicit `profile = "full"` config to
  request new advanced modules.
- Added ENTREZID-aware feature mapping, conservative handling of ambiguous
  uppercase symbols, species-aware mitochondrial prefiltering for Ensembl-only
  inputs, and a recorded cell-level SingleR fallback when optional cluster
  aggregation support is unavailable.
- Documented scRDSreport as the RDS-to-report package stage of the FangLab
  `run_scrnaseq.sh`/`scRNAseq.qmd` workflow. FASTQ alignment, quantification,
  Cell Ranger, nf-core, and Nextflow remain upstream of the package boundary.

# scRDSreport 0.2.0

- Rebuilt the report around the complete twelve-chapter analysis structure of
  the FangLab single-cell Quarto workflow while keeping computation outside the
  document renderer.
- Added a fault-isolated full analysis plan covering QC, dimensionality,
  clustering, cell composition, differential analysis, enrichment, trajectory,
  CellChat, cell cycle, TF expression, CNV preparation/analysis, and downloads.
  Ineligible modules remain visible and explain exactly why they were skipped.
- Added explicit human and mouse resource profiles. Unknown species never fall
  back silently to mouse resources.
- Preserved input annotations by default. Optional annotation modes never
  overwrite an existing RDS annotation column.
- Added module-level parameters, engines, timings, warnings, errors, artifacts,
  SHA-256 hashes, row/column units, and per-column dictionaries to manifest v2.
- Redesigned the standalone Quarto report using the original QMD's numbered
  scientific-documentation layout, with module status, tabbed results,
  colorblind-safe figures, visible download controls, and structured
  interpretation notes below every table and figure.
- Sample-level differential inference now requires defensible biological
  replication; exploratory cell-level comparisons are labelled as such.
- CNV and trajectory modules require an explicit reference/root decision where
  automatic inference would create a biological claim.

# scRDSreport 0.1.0

- Initial public release.
- One-command RDS inspection, SCP completion, export, and Quarto reporting.
- Original and analysis expression matrices use feature-by-cell orientation.
- Existing RDS annotations are preserved without generated cell-type labels.
- UTF-8, DT downloads, ggsci reduction plots, and per-table data dictionaries.
- Conservative sample-group and replicate inference with review flags.
- Fixed Quarto rendering when numeric cluster identifiers are colored with the
  discrete SCI palette.
