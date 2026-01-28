# DKD genetics × single-cell atlas pipeline

Reproducible DKD genetics × kidney single-cell/single-nucleus analysis workflow, with manifest-driven data acquisition and script-generated figures.

## Quick Links
- Data registry: `data/manifest.tsv`
- Data download instructions: `docs/DATA_DOWNLOAD.md`
- Environment notes: `docs/ENVIRONMENT.md`
- Figure index (what tables back each panel): `docs/FIGURE_INDEX.md`
- Cloud run guide (GCP): `docs/GCP.md`

## Repository Layout
- `data/`: external datasets (raw / processed / references), governed by `data/manifest.tsv`
- `scripts/`: numbered, reproducible pipeline scripts
- `results/`: intermediate outputs and “anchor tables” used by figures/manuscript
- `plots/`: figures; final versions in `plots/publication/`
- `docs/`: documentation and reviewer notes

## Reproduce (reviewer-friendly)

Raw datasets are **not** committed. Download/validate them first:

1. `make data` (or follow `docs/DATA_DOWNLOAD.md`)
2. `make preflight` (set `GWAS_ID=...` to choose the primary GWAS; default is `GCST005881`)

Then run the core pipeline:

3. `make preprocess`
4. `make analysis`
5. `make figures`

For a single command end-to-end run:

```bash
make reproduce
```

Figure 1 (workflow schematic) is intentionally not version-locked here; the prompt/spec is tracked at `docs/FIGURE1_GENERATION_PROMPT.md`.

## Cloud run (recommended)

Two practical VM profiles:

- Full rebuild (downloads + Seurat preprocessing + scPagwas + figures): **16 vCPU / 128 GB RAM / ≥300 GB disk**
- Re-run scPagwas only (single-cell objects already exist under `data/processed/`): **8 vCPU / 64 GB RAM / ≥150 GB disk**

To minimize reruns:
- Keep `data/raw/` and `data/processed/` on a persistent disk between sessions.
- The default `make data` excludes two very large CKDGen GWAS files (`GCST008058`, `GCST008064`). To include everything: `make data DOWNLOAD_EXCLUDE_IDS=`.

## OpenGWAS Credentials (MR / coloc)
Some optional stages (e.g., MR, colocalization) query OpenGWAS via API and require a JWT.

- Preferred: write the token to `~/.config/ezhu/opengwas_jwt` (and `chmod 600` the file)
- Alternative: export `OPENGWAS_JWT` in your shell
- Do not commit tokens to the repository; auditors should run with their own OpenGWAS accounts.
