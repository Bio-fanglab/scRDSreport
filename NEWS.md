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
