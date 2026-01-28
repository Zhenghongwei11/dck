# scPagwas Outputs (Source-of-Record Tables)

This directory stores the primary, source-of-record tables used for reporting and figure provenance.

## Runs

### Main run
- Outputs: `results/scpagwas/main/`
- Tables:
  - `merged_celltype_pvalue.csv` (cluster/celltype-level P values)
  - `singlecell_results.csv` (cell-level scores + background-corrected P/adjP/Z)

### Repeat run (same processed atlas)
- Outputs: `results/scpagwas/repeat_run/`
- Tables:
  - `merged_celltype_pvalue.csv`
  - `singlecell_results.csv`
  - `cluster_by_group_counts.tsv` (cell counts by cluster × group for this run)

## Evidence
- Software versions: `results/scpagwas/evidence/software_versions.txt`
- Tail logs are retained under each run directory for traceability.

## Notes
- Cluster IDs can change across runs because clustering is data-dependent. Comparisons should be made at the cell-type/state level using marker-based interpretation, not by assuming cluster numeric IDs are stable.
- The committed runs are both based on the same processed atlas (20,243 cells) and are intended to check numerical stability rather than to represent independent datasets.
