# Running on Google Cloud (Terminal-first)

This project does not require Jupyter. The recommended execution mode is terminal-first on a Google Cloud compute environment with sufficient RAM.

## Option 0: Use Cloud Shell to verify your project (no local install)

1. Open Google Cloud Console and select your billing account / organization.
2. Click the `>_` icon (Cloud Shell).
3. Run:

```bash
gcloud auth list
gcloud projects list
gcloud config get-value project
```

To select a project:

```bash
gcloud config set project YOUR_PROJECT_ID
```

If you need a new project:
- Create it in the console (recommended for beginners), then re-run `gcloud projects list`.

## Option 1: Compute Engine VM (recommended)

1. Create a VM with enough memory for single-cell workloads (choose a high-memory machine type as needed).
2. SSH into the VM.
3. Install prerequisites:
   - `git`, `make`
   - R (matching the version in `renv.lock` if possible)
4. Clone the repository:

```bash
git clone https://github.com/<ORG>/<REPO>.git
cd <REPO>
```

5. Restore R packages and run:

```bash
make setup
make audit AUDIT_LEVEL=ci
make reproduce
```

### Recommended (micromamba) execution

For Seurat/scPagwas, a micromamba-managed R environment is often more stable than system R + renv.

Use the wrapper `scripts/run_r_mamba.sh` for heavy steps (it sets `EZHU_DISABLE_RENV=1`):

```bash
export EZHU_MAMBA_ENV=ezhu-r   # adjust if you use a different env name
export EZHU_SEED=20251227
export EZHU_INTEGRATION_METHOD=seurat

# One-command end-to-end run (downloads + preprocess + scPagwas + anchors + figures + audit)
nohup make reproduce N_CORES=8 RSCRIPT=./scripts/run_r_mamba.sh > logs/reproduce.log 2>&1 &
tail -n 120 logs/reproduce.log
```

If you also want the Word export on the VM:

```bash
make manuscript-docx
```

## Data policy (important)

For third‑party reproducibility, avoid depending on private storage (e.g., personal Drive buckets).
- Register all datasets in `data/manifest.tsv`.
- Prefer authoritative public sources; use paper supplements only when needed.
- Use integrity fields (sha256/md5/size) whenever possible.

## Long-running jobs (avoid session drops)

For multi-hour jobs, do not rely on a single interactive SSH session.

- Prefer `tmux` on the VM: start a session, launch `nohup` jobs, and detach safely.
- Redirect all noisy output to `logs/*.log` and inspect with `tail -n 80 ...` instead of streaming progress bars in chat.

## One-shot full rerun script

If you want a single command that runs the full corrected pipeline end-to-end (download → preprocess → scPagwas → anchors → annotation → figures → audits), use:

```bash
micromamba env create -f env/ezhu-r-seurat4.yml
export EZHU_MAMBA_ENV=ezhu-r-seurat4
export RSCRIPT_BIN=./scripts/run_r_mamba.sh

nohup ./scripts/run_cloud_full.sh > logs/cloud_full_driver.log 2>&1 &
tail -n 120 logs/cloud_full_driver.log
```

Defaults are tuned for stability:
- `EZHU_INTEGRATION_METHOD=seurat`
- `EZHU_N_CORES=1`
- `EZHU_DOWNLOAD_TIMEOUT=7200`

Override as needed, e.g.:

```bash
export EZHU_N_CORES=1
export EZHU_SEED=20251227
export RSCRIPT_BIN=./scripts/run_r_mamba.sh
nohup ./scripts/run_cloud_full.sh > logs/cloud_full_driver.log 2>&1 &
```

## CCC (CellChat) on cloud (Stage 3.4)

Cell–cell communication requires a processed Seurat object (not committed to Git). Run CCC on a VM that already has the Seurat RDS available, then commit only the exported **anchor tables** under `results/ccc/`.

Example (run on the VM, from repo root):

```bash
# Avoid renv auto-activation if you run under micromamba:
export EZHU_DISABLE_RENV=1

# Point to the processed Seurat object on the VM disk:
export EZHU_SEURAT_RDS=/path/to/singlecell_gse131882_seurat.rds

# Optional: reduce runtime / memory by downsampling (reproducible with EZHU_SEED):
export EZHU_SEED=20260102
export EZHU_CCC_MAX_CELLS=20000
export EZHU_CCC_MIN_CELLS=50

nohup make ccc > logs/ccc_cellchat.log 2>&1 &
tail -n 80 logs/ccc_cellchat.log
```

Expected outputs (to commit):
- `results/ccc/<run_id>/cellchat_group_counts.tsv`
- `results/ccc/<run_id>/cellchat_interactions.tsv`
- `results/ccc/<run_id>/cellchat_pathway_summary.tsv`

Large intermediate objects are written under `data/processed/ccc/` by design and must not be committed.

## Post-analysis: Subsampling Stability (Robustness)

If you need to run the 200x subsampling stability analysis after the main scPagwas run:

1. **Environment**: Ensure `r-readr` is installed in the mamba env.
2. **Column Mapping**: The script `scripts/12_scpagwas_subsampling_stability.R` has been patched to handle the column names in modern scPagwas outputs (`cell`, `cluster`, `scPagwas.gPAS.score`).
3. **Execution**:
   ```bash
   export EZHU_SUBSAMPLING_N=200
   export EZHU_SUBSAMPLING_CELL_SCORES=results/scpagwas/repeat_run/cell_scores.tsv
   export EZHU_SUBSAMPLING_SKIP_PLOTS=1  # Avoid X11 errors on cloud
   ./scripts/run_r_mamba.sh scripts/12_scpagwas_subsampling_stability.R
   ```

Expected output: `results/scpagwas/robustness/bootstrap_cluster_stability.tsv`.

## Cost control (recommended)

Compute Engine costs are billed independently of GenAI App Builder promotional credits.

- Keep VMs in `STOPPED` state when idle (stopping ends compute charges; disks may still be billed).
- Set a small budget + alert in the Billing console before running large jobs.
