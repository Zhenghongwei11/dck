# Cell Context Definitions (Primary Dataset: GSE131882)

This document defines the analysis units used in downstream interpretation (cell type vs cell state) and records the current marker-based mapping used for reporting.

## Principle: treat clusters as data-dependent cell states

Numeric cluster labels (e.g., `6`, `7`, `15`) are not guaranteed to be stable across runs because clustering depends on the exact cells retained after QC, normalization, and neighborhood graph construction. For cross-run comparisons (main vs repeat run) and cross-dataset replication, we interpret clusters as *cell contexts* and compare them using marker-based cell type/state labels rather than assuming cluster IDs are comparable.

## Current “cell context” mapping (marker-based)

Source-of-record:
- `results/annotation/cluster_annotations_filled.tsv`
- Marker tables under `results/annotation/markers/`

High-level contexts observed among top-risk clusters:
- **Injured epithelial state(s)**: stress/injury immediate-early and inflammatory programs; interpret as a state that can overlay multiple nephron segments.
- **Loop of Henle / TAL program**: canonical TAL markers support a stable segment identity.
- **Collecting duct programs**: intercalated vs principal-like signatures; confirm with segment markers when needed.
- **Progenitor/repair-like epithelial state (provisional)**: PROM1+ signature; treat as a hypothesis to validate with additional markers and external references.

## Reporting guidance

When writing results:
- Use the cluster numeric ID only as a technical label within a single run.
- Prefer “cell context” names (e.g., *injured tubular epithelial state*) in the narrative.
- When comparing main vs repeat runs, compare the *context label* supported by markers and the direction/strength of enrichment, not the raw cluster IDs.
