# Environment & Setup (WIP)

Goal: enable a third party to reproduce results from scratch with minimal ambiguity.

## Recommendations
- Prefer a single, explicit environment manager for end-to-end runs. In practice this repo supports:
  - **Micromamba** (recommended for full analysis): use the recorded YAML to get a known-good Seurat/scPagwas stack.
  - **renv** (optional): only enable if you intend to maintain a complete `renv.lock` and run `renv::restore()` in CI/clean machines.
- Record `sessionInfo()` for every major analysis stage and write it to `results/`.
- If cloud notebooks are used (e.g., Colab), export executed notebooks and record runtime details (CPU/GPU type, RAM, image/version).

## Cloud execution note (micromamba)

Some single-cell dependencies (notably Seurat/scPagwas combinations) can be brittle on system R installs. For GCP runs we support a micromamba-based R environment:

- Wrapper: `scripts/run_r_mamba.sh`
- It disables `renv` auto-activation by setting `EZHU_DISABLE_RENV=1` and `R_PROFILE_USER=/dev/null`.

For lightweight tasks (e.g., manuscript refresh / DOCX export) you can also disable `renv` locally:

```bash
EZHU_DISABLE_RENV=1 make figures
```

### “Known-good” environment (Seurat 4 + binaries)

If you see failures like:
- `ERROR: dependency 'httpuv' is not available` (Seurat/Shiny toolchain compilation)
- `data.table` compilation failures

Use the recorded, known-good micromamba environment file:

- `env/ezhu-r-seurat4.yml` (R 4.3.3 + Seurat 4.4.0 + prebuilt deps)

On a fresh VM:

```bash
# Install micromamba (one-time)
mkdir -p ~/bin
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
install -m 0755 bin/micromamba ~/bin/micromamba
rm -rf bin
export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook -s bash)"

# Create env (one-time)
micromamba env create -f env/ezhu-r-seurat4.yml

# Use it for the pipeline
export EZHU_MAMBA_ENV=ezhu-r-seurat4
export RSCRIPT_BIN=./scripts/run_r_mamba.sh
```

Then install scPagwas from the pinned upstream repo (this repo DOES NOT use `maialab/scPagwas`):

```bash
$RSCRIPT_BIN scripts/00_install_scpagwas.R
```

## Quick start
- Full analysis (recommended): follow the micromamba section above, then run `make reproduce` with `RSCRIPT_BIN=./scripts/run_r_mamba.sh`.
- renv-based local restore (optional): `make setup` then `make reproduce` (requires a complete `renv.lock` and may require compilation).
