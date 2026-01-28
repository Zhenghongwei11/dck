# Data Download Guide

This repository intentionally does **not** commit raw datasets under `data/raw/` (they are large and publicly available). Instead, we provide a machine-readable registry (`data/manifest.tsv`) plus scripts to download/validate the exact inputs used in this study.

## What you need

- `bash`, `curl`, `tar`
- A checksum tool:
  - Linux: `md5sum`
  - macOS: `md5`
- Sufficient disk space for GEO archives and processed single-cell objects.

## Option A (recommended): automated download from the manifest

From the repo root:

```bash
make data
```

This runs `scripts/01_download_data.R`, which:
- Reads `data/manifest.tsv`
- Downloads entries with `download_method=script` into their `local_path`
- Validates checksums where provided

Notes:
- `make data` excludes two very large CKDGen GWAS files by default (`GCST008058`, `GCST008064`). To include everything, run:
  - `make data DOWNLOAD_EXCLUDE_IDS=`
- To download only specific manifest IDs, set `EZHU_DOWNLOAD_IDS` (comma-separated), e.g.:
  - `EZHU_DOWNLOAD_IDS=GCST005881,GSE131882 make data`

After download, run a quick GWAS sanity check:

```bash
make preflight
```

Notes:
- Use `GWAS_ID=...` to choose the primary GWAS row from `data/manifest.tsv` (default is `GCST005881`):
  - `make preflight GWAS_ID=GCST005881`
  - `make analysis GWAS_ID=GCST005881`
- The preflight stage writes the harmonized file to `data/processed/gwas_harmonized__<GWAS_ID>.tsv.gz` and the analysis stage is configured to read the same file (see `Makefile`).

## Option B: manual download (exact paths)

If you prefer manual download (or your environment blocks scripted downloads), download each URL in `data/manifest.tsv` and place it at the exact `local_path`. The pipeline assumes these paths.

### 1) GWAS summary statistics (GWAS Catalog / EBI FTP)

Create the target folder:

```bash
mkdir -p data/raw/gwas
```

Download:

```bash
curl -L -o data/raw/gwas/GCST005881.h.tsv.gz \
  "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005881/harmonised/29703844-GCST005881-EFO_0000401.h.tsv.gz"

curl -L -o data/raw/gwas/GCST008058.dbgap.txt.gz \
  "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008058/20171016_MW_eGFR_overall_ALL_nstud61.dbgap.txt.gz"

curl -L -o data/raw/gwas/GCST008064.dbgap.txt.gz \
  "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008064/CKD_overall_ALL_JW_20180223_nstud30.dbgap.txt.gz"
```

Checksum verification (examples):

- Linux:
  - `md5sum data/raw/gwas/GCST005881.h.tsv.gz`
- macOS:
  - `md5 data/raw/gwas/GCST005881.h.tsv.gz`

Expected MD5 values are recorded in `data/manifest.tsv`.

Optional (legacy, underpowered; not used in the current primary analysis):

```bash
curl -L -o data/raw/gwas/GCST90435703.tsv.gz \
  "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435703/GCST90435703.tsv.gz"
```

### 2) Single-cell data (GEO FTP)

Create folders:

```bash
mkdir -p data/raw/singlecell/GSE131882
mkdir -p data/raw/singlecell/GSE195460/h5
```

GSE131882 (tar archive):

```bash
curl -L -o data/raw/singlecell/GSE131882/GSE131882_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131882/suppl/GSE131882_RAW.tar"
```

Unpack into the expected structure:

```bash
make unpack-gse131882
```

GSE195460 (H5 matrices; download the files listed in `data/manifest.tsv`):

```bash
# Example (Control1); repeat for other *_filtered_feature_bc_matrix.h5 files in the manifest
curl -L -o data/raw/singlecell/GSE195460/h5/GSE195460_Control1_filtered_feature_bc_matrix.h5 \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE195nnn/GSE195460/suppl/GSE195460_Control1_filtered_feature_bc_matrix.h5"
```

## After download: reproduce the analysis

Minimum (preprocessing + core scPagwas + figures):

```bash
make preprocess
make analysis
make figures
```

Notes:
- The MR/colocalization stages are optional and may require an OpenGWAS token; see `README.md` for the JWT setup.
- Figure 1 is a schematic generated externally; Figures 2–7 are produced by scripts under `scripts/figures/`.
