#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# This script reproduces the local (non-cloud) stages:
# - Publication figures (from source-of-record tables)
# - CI-level audits (no raw data required)
#
# It intentionally does NOT re-run cloud-heavy steps (preprocess/scPagwas).

export EZHU_DISABLE_RENV="${EZHU_DISABLE_RENV:-1}"
export EZHU_SEED="${EZHU_SEED:-20251227}"

echo "[local] repo_root=$REPO_ROOT"
echo "[local] EZHU_DISABLE_RENV=$EZHU_DISABLE_RENV"
echo "[local] EZHU_SEED=$EZHU_SEED"

echo "[local] generating figures"
make figures

echo "[local] exporting causal candidates (no heavy data required)"
make causal

echo "[local] running audits (ci)"
make audit AUDIT_LEVEL=ci

echo "[local] done"
