# Cloud Install Recipe (micromamba + Seurat 4 + scPagwas)

This document captures a **known-good** installation + execution path for running the full scPagwas pipeline on a fresh Linux VM without getting stuck in CRAN source compilation (e.g., `httpuv`, `data.table`, `ggpubr`).

## Why this exists

- **Seurat installed from CRAN** often triggers a cascade of source builds (`httpuv`, `data.table`, etc.) that can fail depending on system libraries and compiler toolchain.
- **This repo contains a micromamba environment file** that installs the heavy R stack from **conda-forge binaries** (Seurat 4 + dependencies) to avoid those failures.
- **This repo uses `renv` for local work**, but on cloud runs we intentionally disable `renv` so that the conda environment’s library paths are used.

### scPagwas compatibility note (Seurat v5)

`scPagwas` is not compatible with SeuratObject v5 because `GetAssayData(slot=...)` is defunct in v5. The supported path for this project is:
- Seurat < 5.0.0
- SeuratObject < 5.0.0

If you see errors like `GetAssayData(): slot argument ... is defunct`, your environment has drifted to Seurat v5; recreate a clean micromamba env from `env/ezhu-r-seurat4.yml`.

## Source-of-truth files

- Environment (Seurat 4, binaries): `env/ezhu-r-seurat4.yml`
- R runner wrapper (disables `renv` + disables `.Rprofile`): `scripts/run_r_mamba.sh`
- scPagwas installer (pinned upstream commit): `scripts/00_install_scpagwas.R`
- One-shot end-to-end runner: `scripts/run_cloud_full.sh`

## 0) VM recommendation

- RAM: **96 GiB** (single-cell integration + scPagwas stability)
- CPU: 1 core is the safest default (`EZHU_N_CORES=1`)
- Disk: ensure enough space for `data/raw/` downloads and `data/processed/` artifacts

## 1) Clone the repo (SSH)

```bash
This project is intended to be run from the local repository checkout (this folder). If you are working from a shared archive (e.g., a review bundle), extract it first and then `cd` into the extracted directory.
git checkout main
git pull --ff-only
```

## 2) Install micromamba (one-time)

```bash
mkdir -p ~/bin
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
install -m 0755 bin/micromamba ~/bin/micromamba
rm -rf bin
export PATH="$HOME/bin:$PATH"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$(micromamba shell hook -s bash)"
```

## 3) Create the known-good R environment (one-time)

Optional (fast) sanity-check that every dependency name in `env/ezhu-r-seurat4.yml` exists on conda channels:

```bash
./scripts/audit/check_env_conda_packages.py --file env/ezhu-r-seurat4.yml
```

```bash
micromamba env create -f env/ezhu-r-seurat4.yml -y || \
  micromamba env update -f env/ezhu-r-seurat4.yml -n ezhu-r-seurat4
```

## 4) Run the pipeline (recommended: no `activate`)

Avoid `micromamba activate ...` in strict shells because activation scripts can reference unset variables and fail under `set -u`.

Use the wrapper (it uses `micromamba run -n ...` and sets `EZHU_DISABLE_RENV=1`):

```bash
export EZHU_MAMBA_ENV=ezhu-r-seurat4
export RSCRIPT_BIN=./scripts/run_r_mamba.sh

# Stability defaults
export EZHU_SEED=20251227
export EZHU_INTEGRATION_METHOD=seurat
export EZHU_N_CORES=1
export EZHU_DOWNLOAD_TIMEOUT=7200

nohup ./scripts/run_cloud_full.sh > logs/cloud_full_driver.log 2>&1 &
tail -n 120 logs/cloud_full_driver.log
```

Note: `scPagwas` is installed from GitHub and is not part of the base conda env. If you created a fresh env and run steps manually (instead of the full runner), install it once into the env-scoped user library:

```bash
export EZHU_MAMBA_ENV=ezhu-r-seurat4-clean
./scripts/run_r_mamba.sh scripts/00_install_scpagwas.R
```

### Optional: a single helper script

If you want a single helper that (i) creates/updates the env from `env/ezhu-r-seurat4.yml`, then (ii) starts the one-shot run:

```bash
chmod +x scripts/run_seurat4_mamba_full.sh
nohup ./scripts/run_seurat4_mamba_full.sh > logs/run_seurat4_mamba_full.log 2>&1 &
tail -n 120 logs/run_seurat4_mamba_full.log
```

## 5) Verify the environment quickly

**Do not** run plain `Rscript` from the repo root unless you also pass `--vanilla` and disable `renv`, otherwise `.Rprofile` may activate `renv` and hide the conda libraries.

Good checks:

```bash
# Using micromamba directly:
micromamba run -n ezhu-r-seurat4 env EZHU_DISABLE_RENV=1 Rscript --vanilla -e 'stopifnot(requireNamespace("Seurat", quietly=TRUE))'
micromamba run -n ezhu-r-seurat4 env EZHU_DISABLE_RENV=1 Rscript --vanilla -e 'stopifnot(requireNamespace("scPagwas", quietly=TRUE))'

# Or using the repo wrapper:
./scripts/run_r_mamba.sh -e 'stopifnot(requireNamespace("Seurat", quietly=TRUE)); cat("Seurat OK\\n")'
./scripts/run_r_mamba.sh -e 'stopifnot(requireNamespace("scPagwas", quietly=TRUE)); cat("scPagwas OK\\n")'
```

Recommended version guard (run via the wrapper so `renv` cannot interfere):

```bash
export EZHU_MAMBA_ENV=ezhu-r-seurat4-clean
./scripts/run_r_mamba.sh -e 'cat("Seurat:", as.character(packageVersion("Seurat")), "\n");
  cat("SeuratObject:", as.character(packageVersion("SeuratObject")), "\n");
  stopifnot(packageVersion("Seurat") < "5.0.0");
  stopifnot(packageVersion("SeuratObject") < "5.0.0")'
```

## 6) Common pitfalls (and fixes)

### Pitfall A: `requireNamespace("ggpubr") is not TRUE` even after `micromamba install`

Cause: `Rscript` loaded `.Rprofile` → activated `renv` → `.libPaths()` points to `renv/library/` and hides the conda libraries.

Fix:
- Use `./scripts/run_r_mamba.sh ...` (recommended), or
- Use `micromamba run -n ezhu-r-seurat4 env EZHU_DISABLE_RENV=1 Rscript --vanilla ...`

### Pitfall B: `micromamba activate` fails under `set -u` (unbound variables)

Fix:
- Prefer `micromamba run -n ...` (recommended), or
- Temporarily disable nounset around activation:
  - `set +u; micromamba activate ezhu-r-seurat4; set -u`

### Pitfall C: scPagwas GitHub 404

This repo installs scPagwas from:
- `sulab-wmu/scPagwas@8dd7975e262a4400928601389142b06cb78367a0`

Do **not** use `maialab/scPagwas` (404 / wrong repo).

### Pitfall D: Bioconductor version mismatch after conda upgrades

If the conda environment upgrades `r-base` (e.g., from R 4.4 → 4.5), BiocManager may error like:

> `Bioconductor version '3.20' requires R version '4.4'`

Fix options:
- Re-run `./scripts/run_r_mamba.sh scripts/00_install_scpagwas.R` (it auto-selects a compatible Bioconductor version), or
- Set an explicit Bioconductor version for the installer:
  - `EZHU_BIOC_VERSION=3.22 ./scripts/run_r_mamba.sh scripts/00_install_scpagwas.R`

### Pitfall D2: `rtracklayer` fails to compile from source (modern GCC)

This only applies if you explicitly enable the **VCF backend** for colocalization (set `EZHU_COLOC_USE_VCF=1`). The default coloc backend uses the OpenGWAS API and does **not** require `rtracklayer`/`VariantAnnotation`.

Symptom (during coloc VCF backend installs, e.g. `gwasvcf` / `VariantAnnotation`):

> `error: too many arguments to function 'free' ...` (or other low-level UCSC/rtracklayer compile errors)

Cause:
- R/Bioc is attempting to compile `rtracklayer` (and related Bioconductor packages) from source.
- On some modern toolchains (e.g., GCC 15.x), this can fail due to upstream C code compatibility.

Fix (recommended; avoid source compilation):
1. Ensure you are using the correct conda channels (do **not** use a non-existent channel like `bioc`):
   - Use `conda-forge` and `bioconda`.
2. Install Bioconductor packages as **prebuilt conda binaries**:

```bash
export EZHU_MAMBA_ENV=ezhu-r-seurat4-matrix161
micromamba install -y -n "$EZHU_MAMBA_ENV" -c conda-forge -c bioconda \
  bioconductor-variantannotation bioconductor-genomicranges \
  bioconductor-rtracklayer bioconductor-rsamtools bioconductor-rhtslib
```

Then install the remaining coloc R deps via the wrapper:

```bash
EZHU_MAMBA_ENV=$EZHU_MAMBA_ENV ./scripts/run_r_mamba.sh scripts/00_install_coloc_deps.R
```

### Pitfall E: `GetAssayData(slot=...) ... is now defunct` (SeuratObject v5)

Symptom (typically during `make analysis` / scPagwas):

> `Error: The slot argument of GetAssayData() was deprecated ... and is now defunct`

Cause: the R environment drifted to **Seurat v5 / SeuratObject v5**, but scPagwas requires **Seurat v4-era APIs**.

Hard requirement for this project:
- `Seurat < 5.0.0`
- `SeuratObject < 5.0.0`

Fix:
1. **Create a clean env** from `env/ezhu-r-seurat4.yml` (recommended to use a new env name, not mutate an existing one):

```bash
export EZHU_MAMBA_ENV=ezhu-r-seurat4-clean
micromamba env create -f env/ezhu-r-seurat4.yml -n $EZHU_MAMBA_ENV -y || \
  micromamba env update -f env/ezhu-r-seurat4.yml -n $EZHU_MAMBA_ENV -y
```

2. Verify versions:

```bash
EZHU_MAMBA_ENV=$EZHU_MAMBA_ENV ./scripts/run_r_mamba.sh -e 'cat("Seurat:", as.character(packageVersion("Seurat")), "\n"); cat("SeuratObject:", as.character(packageVersion("SeuratObject")), "\n")'
```

3. **Re-run preprocessing** to regenerate Seurat RDS artifacts in a compatible format:

```bash
EZHU_MAMBA_ENV=$EZHU_MAMBA_ENV make preprocess N_CORES=1 RSCRIPT=./scripts/run_r_mamba.sh
```

Notes:
- This repo’s installer `scripts/00_install_scpagwas.R` intentionally does **not** install/upgrade Seurat via CRAN to avoid accidentally pulling v5.
- Always run R scripts via `./scripts/run_r_mamba.sh` so user-level packages cannot shadow conda packages.

### Pitfall F: Conda packages installed, but R still loads the “wrong” ones

Symptom: you installed packages via micromamba, but `Rscript` still cannot `library(...)`, or loads unexpected versions.

Cause: user-level libraries (e.g. `R_LIBS_USER`) can shadow conda libraries.

Fix:
- Use the wrapper `./scripts/run_r_mamba.sh` (it sets `EZHU_DISABLE_RENV=1`, disables `.Rprofile`, and isolates `R_LIBS_USER`).
- Avoid plain `Rscript` from the repo root.

Extra note: even calling the **conda env’s** `Rscript` directly can still read the repo `.Rprofile` if you run it from the repo root. Always use:
- `./scripts/run_r_mamba.sh ...` (recommended), or
- `micromamba run -n <env> env EZHU_DISABLE_RENV=1 Rscript --vanilla ...`

Recovery: if you ever ran `install.packages("SeuratObject")` (or `Seurat`, or `Matrix`) and your env “drifted”:
- The wrapper uses an env-scoped user library under `~/.cache/ezhu/r_mamba_userlib/<ENV_NAME>/`.
- If those core packages exist there, delete them so they cannot shadow the micromamba stack:
  - `rm -rf ~/.cache/ezhu/r_mamba_userlib/<ENV_NAME>/SeuratObject ~/.cache/ezhu/r_mamba_userlib/<ENV_NAME>/Seurat ~/.cache/ezhu/r_mamba_userlib/<ENV_NAME>/Matrix`

### Pitfall G: Anchors export fails on Seurat v5 “split layers”

Seurat v5 can store `RNA` assay matrices split across layers like `counts.GSM...` / `data.GSM...`.

Status: the repo scripts were updated to handle split layers without calling `JoinLayers()`:
- QC anchors: `scripts/02d_export_qc_anchors.R`
- UMAP + feature anchors: `scripts/02c_export_umap_anchors.R`

Fix:
- `git pull --ff-only` and re-run:
  - `make anchors-qc anchors-umap RSCRIPT=./scripts/run_r_mamba.sh`

### Pitfall H: GSE131882 unpack checksum mismatch (`.rds.gz` vs `.rds`)

Cause: the GEO raw tar contains gzip-compressed RDS files; the index/checksum table may record MD5 for the uncompressed `.rds`.

Fix:
- Ensure you are on the latest `main` and rerun:
  - `make unpack-gse131882 RSCRIPT=./scripts/run_r_mamba.sh`

### Pitfall I: Multiple cloud runs collide (duplicate processes)

Symptom: logs look interleaved, stages rerun unexpectedly, or output files change during a run.

Fix:
- Start exactly one driver process:
  - `nohup ./scripts/run_cloud_full.sh > logs/cloud_full_driver.log 2>&1 &`
- Tail logs, do not re-run the driver while it’s still running:
  - `tail -f logs/preprocess.log`

### Pitfall J: Threading / cores cause instability

For stability on large Seurat objects, prefer single-core + single-thread BLAS:

```bash
export EZHU_N_CORES=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLAS_NUM_THREADS=1
```

### Pitfall K: `invalid class “Graph” object: superclass "mMatrix" not defined`

Symptom (typically during `make preprocess` clustering / neighbor graph construction):

> `invalid class "Graph" object: superclass "mMatrix" not defined`

Cause: `Matrix` is too new (>= 1.7). In Matrix 1.7 the S4 class `mMatrix` was removed, but Seurat v4 Graph objects still reference it.

Fix:
- Ensure the micromamba env uses the pinned `r-matrix` version from `env/ezhu-r-seurat4.yml`, then recreate or update the env:

```bash
micromamba env update -f env/ezhu-r-seurat4.yml -n $EZHU_MAMBA_ENV -y
```

Quick check (via wrapper):

```bash
./scripts/run_r_mamba.sh -e 'library(Matrix); cat("Matrix:", as.character(packageVersion("Matrix")), "\n"); stopifnot(packageVersion("Matrix") < "1.7-0"); stopifnot(methods::isClass("mMatrix")); stopifnot(methods::getClassDef("mMatrix")@package == "Matrix")'
```

If you still see this error even with `Matrix < 1.7`, it usually means the environment is inconsistent (e.g., mixed libraries or a broken Matrix install). Update the repo and recreate a clean env; do not “monkey patch” `mMatrix` in R:
- `git pull --ff-only`
- `micromamba env create -f env/ezhu-r-seurat4.yml -n ezhu-r-seurat4-clean -y`
- `make preprocess RSCRIPT=./scripts/run_r_mamba.sh`

### Pitfall L: `Missing package: hdf5r` while preprocessing GSE195460

Cause: GSE195460 is stored as 10x `.h5` and requires `hdf5r` plus HDF5 system libs.

Fix:
- Update/recreate the micromamba env from `env/ezhu-r-seurat4.yml` (it includes `r-hdf5r`, `hdf5`, `zlib`, `libaec`), or install into the current env:

```bash
micromamba install -y -n $EZHU_MAMBA_ENV -c conda-forge r-hdf5r hdf5 zlib libaec
```

### Pitfall M: Compilation errors like `cannot find -lz` or `cannot find -llzma`

Cause: during `install_github` or CRAN source installs (e.g., `igraph`, `svglite`), the R compiler cannot find system compression libraries even if they exist in the conda prefix.

Fix:
- Install development headers and libraries explicitly via micromamba:
```bash
micromamba install -y -n $EZHU_MAMBA_ENV -c conda-forge zlib xz libxml2
```
- If issues persist, install the R-specific binary version of the problematic package:
```bash
micromamba install -y -n $EZHU_MAMBA_ENV -c conda-forge r-igraph r-svglite
```

### Pitfall N: CellChat error: `Cell labels cannot contain '0'!`

Cause: CellChat ident validation fails if any cluster/cell-type label is a pure numeric string containing '0' (common with default `seurat_clusters`).

Fix:
- Prepend a character prefix to all identity labels before creating the CellChat object:
```r
# Example fix in R
obj$cell_context <- paste0("C", as.character(obj$seurat_clusters))
```

### Pitfall O: Ad-hoc `micromamba install` causing Seurat v5 upgrade

Cause: running `micromamba install <new-package>` without pinning existing core dependencies can cause conda to solve for the latest available versions, inadvertently upgrading `r-seurat` from 4.x to 5.x.

Fix:
- Always include explicit version pins for core stack members in every install command:
```bash
micromamba install -y -n $EZHU_MAMBA_ENV -c conda-forge \
  r-base=4.3.3 r-seurat=4.4.0 r-matrix=1.6.4 <new-package-name>
```

## 7) Audit-friendly record keeping

The pipeline writes:
- Stage metadata: `results/metadata/*`
- Logs: `logs/*.log`

To preserve an “audit bundle” in GitHub:
- The full runner `scripts/run_cloud_full.sh` writes an audit bundle on exit under `docs/audit_runs/<run_id>/` (logs + environment metadata + artifact hashes).
- You can also generate a bundle manually:
  - `make audit-bundle`

Commit (recommended):
- `docs/audit_runs/<run_id>/`
- `results/figures/`, `results/scpagwas/`, `results/annotation/`
- `plots/publication/`
- `docs/MANUSCRIPT_*.md` and `docs/MANUSCRIPT_SUBMISSION.docx`
