# Pipeline overview

This document summarizes the end‑to‑end workflow and how `scripts/` connect.

## One‑command reproduction
- Full run: `make reproduce`
- Audits only: `make audit` (set `AUDIT_LEVEL=ci` for lightweight checks)
- GWAS preflight: `make preflight`
- Local-only reproduction (no cloud compute): `bash scripts/reproduce_local.sh`

## What `make reproduce` includes
- Download manifest-driven inputs via `make data` (see `data/manifest.tsv` and `docs/DATA_DOWNLOAD.md`).
- Unpack GSE131882 (GEO tar) into curated per-sample files (`make unpack-gse131882`).
- Preprocess both datasets:
  - `preprocess-gse131882`
  - `preprocess-gse195460`
- Run primary scPagwas twice for stability:
  - `analysis-primary` → `results/scpagwas/main/` and `results/scpagwas/repeat_run/`
- Run replication scPagwas:
  - `analysis-replication` → `results/scpagwas/gse195460/`
- Regenerate figure anchor tables from the current Seurat objects:
  - QC anchors: `make anchors-qc`
  - UMAP anchors: `make anchors-umap`
- Generate publication figures: `make figures`

## Principles
- Run from repository root with **relative paths only**.
- Key numerical outputs are written to `results/` as CSV/TSV “anchor tables”.
- Figures are generated from scripts and saved under `plots/publication/`.

## Cloud execution note (micromamba)
For cloud VMs we support running R inside a micromamba environment to avoid system-level dependency issues.

- Wrapper: `scripts/run_r_mamba.sh`
- Typical pattern:
  - `./scripts/run_r_mamba.sh scripts/02_preprocess.R`
  - `./scripts/run_r_mamba.sh scripts/03_analysis.R`

