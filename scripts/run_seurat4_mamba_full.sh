#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

mkdir -p logs

MAMBA_BIN="${MAMBA_BIN:-}"
if [[ -z "${MAMBA_BIN}" ]]; then
  if command -v micromamba >/dev/null 2>&1; then
    MAMBA_BIN="micromamba"
  else
    MAMBA_BIN="$HOME/bin/micromamba"
  fi
fi

if [[ ! -x "${MAMBA_BIN}" ]] && ! command -v "${MAMBA_BIN}" >/dev/null 2>&1; then
  echo "micromamba not found. Install it first (see docs/INSTALL_MICROMAMBA_SEURAT4.md)." >&2
  exit 127
fi

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$("${MAMBA_BIN}" shell hook -s bash)"

ENV_NAME="${EZHU_MAMBA_ENV:-ezhu-r-seurat4}"

echo "[run-seurat4] repo_root=$REPO_ROOT"
echo "[run-seurat4] MAMBA_BIN=$MAMBA_BIN"
echo "[run-seurat4] ENV_NAME=$ENV_NAME"

echo "[run-seurat4] ensure env from env/ezhu-r-seurat4.yml"
if "${MAMBA_BIN}" env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  "${MAMBA_BIN}" env update -f env/ezhu-r-seurat4.yml -n "${ENV_NAME}" -y > logs/mamba_env_update.log 2>&1
else
  "${MAMBA_BIN}" env create -f env/ezhu-r-seurat4.yml -y > logs/mamba_env_create.log 2>&1
fi

# Run everything via the wrapper to avoid renv/.Rprofile interference.
export EZHU_MAMBA_ENV="${ENV_NAME}"
export RSCRIPT_BIN="${RSCRIPT_BIN:-./scripts/run_r_mamba.sh}"

export EZHU_SEED="${EZHU_SEED:-20251227}"
export EZHU_INTEGRATION_METHOD="${EZHU_INTEGRATION_METHOD:-seurat}"
export EZHU_N_CORES="${EZHU_N_CORES:-1}"
export EZHU_DOWNLOAD_TIMEOUT="${EZHU_DOWNLOAD_TIMEOUT:-7200}"

echo "[run-seurat4] starting full pipeline via scripts/run_cloud_full.sh"
nohup ./scripts/run_cloud_full.sh > logs/cloud_full_driver.log 2>&1 &
echo $! > logs/cloud_full_driver.pid
echo "[run-seurat4] pid=$(cat logs/cloud_full_driver.pid)"
echo "[run-seurat4] tail -n 120 logs/cloud_full_driver.log"

