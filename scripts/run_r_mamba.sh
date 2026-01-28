#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${EZHU_MAMBA_ENV:-ezhu-r}"

MAMBA_BIN="${MAMBA_BIN:-}"
if [[ -z "${MAMBA_BIN}" ]]; then
  if command -v micromamba >/dev/null 2>&1; then
    MAMBA_BIN="micromamba"
  else
    MAMBA_BIN="${HOME}/bin/micromamba"
  fi
fi

if [[ ! -x "${MAMBA_BIN}" ]] && ! command -v "${MAMBA_BIN}" >/dev/null 2>&1; then
  echo "micromamba not found. Set MAMBA_BIN or ensure micromamba is on PATH." >&2
  exit 127
fi

export EZHU_DISABLE_RENV=1
export R_PROFILE_USER="${R_PROFILE_USER:-/dev/null}"
export R_ENVIRON_USER="${R_ENVIRON_USER:-/dev/null}"
# Use an env-scoped user library so GitHub/CRAN installs (e.g. scPagwas) persist
# without polluting other environments, and to avoid accidental shadowing.
export R_LIBS_USER="${R_LIBS_USER:-$HOME/.cache/ezhu/r_mamba_userlib/${ENV_NAME}}"
mkdir -p "$R_LIBS_USER" >/dev/null 2>&1 || true

# Hard guardrails: these core packages must come from the micromamba env (Seurat v4 stack).
for pkg in Seurat SeuratObject Matrix; do
  if [[ -d "${R_LIBS_USER}/${pkg}" ]]; then
    echo "ERROR: R user library contains '${pkg}' at: ${R_LIBS_USER}/${pkg}" >&2
    echo "This will shadow the micromamba-installed versions and can break the pipeline (e.g., SeuratObject v5)." >&2
    echo "Fix: remove it (rm -rf '${R_LIBS_USER}/${pkg}') and re-run, or recreate the micromamba env." >&2
    exit 2
  fi
done

exec "${MAMBA_BIN}" run -n "${ENV_NAME}" Rscript --vanilla "$@"
