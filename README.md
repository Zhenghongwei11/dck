# DKD genetics × kidney snRNA-seq integration (reviewer repo)

This repository contains code and lightweight result tables for a DKD genetics × kidney single-nucleus transcriptomic workflow (scPagwas integration + robustness checks), sufficient to regenerate the main figures from the provided “source-of-record” tables.

## What is included
- `scripts/`: reproducible analysis + figure scripts
- `results/`: derived tables used by figures (no raw data)
- `plots/publication/`: final figures (PDF + PNG)
- `data/manifest.tsv` + `docs/DATA_DOWNLOAD.md`: public-data download pointers for full re-runs

## What is not included
- Manuscript / submission files (submitted separately)
- Large raw/processed single-cell objects and full per-pathway regression dumps (cloud compute recommended)

## Quick start (regenerate figures from included tables)

```bash
EZHU_DISABLE_RENV=1 make figures FIGURES=2,3,4,5,6,7
```

Outputs are written to `plots/publication/` and `plots/publication/png/`.

## Full re-run (cloud recommended)

See `docs/GCP.md` and `docs/INSTALL_MICROMAMBA_SEURAT4.md`. For a full rebuild (downloads + preprocessing + scPagwas + figures), a practical profile is **16 vCPU / 128 GB RAM / ≥300 GB disk**.

## OpenGWAS credentials (MR / coloc)
Some optional stages (e.g., MR, colocalization) query OpenGWAS via API and require a JWT.

- Preferred: write the token to `~/.config/ezhu/opengwas_jwt` (and `chmod 600` the file)
- Alternative: export `OPENGWAS_JWT` in your shell
- Do not commit tokens to the repository.
