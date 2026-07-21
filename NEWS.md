# scRDSreport 0.3.2

- Added `dependency_status()`, `install_dependencies()`, and
  `check_dependencies()` so a requested report profile can be installed and
  verified by actually loading every required namespace rather than checking
  package-directory presence alone.
- Split installation across CRAN, Bioconductor, and declared GitHub sources,
  retained the active Bioconductor repositories while installing SCP, and
  added explicit GitHub routes for CellChat and current Monocle 3.
- Added input-aware species dependency selection. `species = "auto"` now needs
  an input RDS before an OrgDb/TxDb is selected; `species = "all"` remains an
  explicit opt-in for all nine built-in organisms.
- Added separate checks for the external Quarto CLI, JAGS, and the HDF5 build
  prerequisite. Installing the R package named `quarto` is no longer presented
  as installing the CLI.
- Removed unused direct `Suggests` entries that are either not called by this
  package or already belong to an analysis engine's transitive dependency
  declaration, avoiding fragile `dependencies = TRUE` over-installation.
- Added a standalone bootstrap script for a fresh R session, including a
  writable user-library fallback and a final strict loadability check.
- Raised the declared baseline to R 4.2 so it matches the currently required
  SCP dependency generation; the newest Monocle 3 path has stricter limits
  documented by the installer.

# scRDSreport 0.3.1

- Made useful, scientifically bounded fallback output the default for advanced
  modules. When reliable biological replication is unavailable, differential
  analysis now exports bounded annotation/cluster one-versus-rest rankings and
  descriptive pooled-CPM/log2-fold-change tables and figures without P values
  or FDR; pooled CPM is named explicitly and is not presented as a sample mean.
- Added descriptive gene-set effect summaries that can consume exploratory
  differential rankings. These summaries report overlap, direction, and effect
  sizes only and are never presented as significance tests.
- Kept biological decisions explicit while still producing diagnostic output:
  trajectory analysis exports unrooted geometry and root candidates without an
  explicit root; communication exports grouping context and diagnostics when
  inference is unavailable; CNV exports readiness information without guessing
  a normal reference or producing a CNV signal.
- Added a genome-build gate before inferCNV can use a built-in TxDb, and retain
  readiness/input artifacts when inferCNV object creation or execution fails.
- Reused complete cell-cycle scores already present in an input RDS before
  requiring external species resources.
- Fixed annotation and object-audit tabs that could display generated HTML as
  literal source text. The left report table of contents now starts expanded.
- Improved the complete-file index with a wider explanation column and a fixed
  per-row scroll area, and made fallback table descriptions distinguish
  descriptive output from formal inference.

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
- Documented scRDSreport as an RDS-to-report R package that inspects supported
  single-cell objects, optionally completes analysis, exports data products,
  and renders the final Quarto report.

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
