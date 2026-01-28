#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

EXIT_CODE="${1:-0}"

GIT_SHA="$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")"
UTC_TS="$(date -u +"%Y%m%dT%H%M%SZ" 2>/dev/null || date +"%Y%m%dT%H%M%S")"
RUN_ID="${EZHU_RUN_ID:-${UTC_TS}_${GIT_SHA}}"

AUDIT_ROOT="${EZHU_AUDIT_ROOT:-docs/audit_runs}"
RUN_DIR="${EZHU_AUDIT_RUN_DIR:-${AUDIT_ROOT}/${RUN_ID}}"

mkdir -p "$RUN_DIR"/{meta,logs,results,plots,docs,env}

{
  echo "run_id=${RUN_ID}"
  echo "exit_code=${EXIT_CODE}"
  echo "git_sha=${GIT_SHA}"
  echo "utc_ts=${UTC_TS}"
} > "$RUN_DIR/meta/run_meta.txt"

{
  echo "repo_root=$REPO_ROOT"
  echo "pwd=$(pwd)"
  echo
  echo "== git =="
  git status --porcelain=v1 || true
  echo
  git log -n 1 --oneline || true
} > "$RUN_DIR/meta/git_state.txt"

{
  echo "== os-release =="
  if [[ -f /etc/os-release ]]; then
    cat /etc/os-release
  else
    uname -a || true
  fi
  echo
  echo "== cpu =="
  if command -v lscpu >/dev/null 2>&1; then
    lscpu || true
  else
    sysctl -a 2>/dev/null | grep -E "machdep\\.cpu|hw\\.memsize" || true
  fi
  echo
  echo "== memory =="
  if command -v free >/dev/null 2>&1; then
    free -h || true
  else
    vm_stat 2>/dev/null || true
  fi
  echo
  echo "== disk =="
  df -h || true
} > "$RUN_DIR/meta/machine.txt"

{
  echo "== selected env vars =="
  env | grep -E "^(EZHU_|RSCRIPT_BIN=|EZHU_MAMBA_ENV=|MAMBA_|R_PROFILE_USER=|R_LIBS=|R_LIBS_USER=|OMP_NUM_THREADS=|OPENBLAS_NUM_THREADS=|MKL_NUM_THREADS=|VECLIB_MAXIMUM_THREADS=|BLAS_NUM_THREADS=)" | sort || true
} > "$RUN_DIR/meta/env_vars.txt"

if [[ -f "env/ezhu-r-seurat4.yml" ]]; then
  cp -f "env/ezhu-r-seurat4.yml" "$RUN_DIR/env/ezhu-r-seurat4.yml"
fi
if [[ -f "renv.lock" ]]; then
  cp -f "renv.lock" "$RUN_DIR/env/renv.lock"
fi

if command -v micromamba >/dev/null 2>&1 && [[ -n "${EZHU_MAMBA_ENV:-}" ]]; then
  {
    echo "== micromamba env list =="
    micromamba env list || true
    echo
    echo "== micromamba list (${EZHU_MAMBA_ENV}) =="
    micromamba list -n "${EZHU_MAMBA_ENV}" || true
    echo
    echo "== micromamba env export (${EZHU_MAMBA_ENV}) =="
    micromamba env export -n "${EZHU_MAMBA_ENV}" || true
  } > "$RUN_DIR/meta/micromamba_${EZHU_MAMBA_ENV}.txt"
fi

RSCRIPT_BIN="${RSCRIPT_BIN:-}"
if [[ -z "$RSCRIPT_BIN" ]]; then
  if [[ -x "./scripts/run_r_mamba.sh" ]]; then
    RSCRIPT_BIN="./scripts/run_r_mamba.sh"
  else
    RSCRIPT_BIN="Rscript"
  fi
fi

{
  echo "== R sessionInfo() =="
  if "$RSCRIPT_BIN" -e 'sessionInfo()' 2>/dev/null; then
    true
  else
    echo "(failed to run sessionInfo via RSCRIPT_BIN=$RSCRIPT_BIN)"
  fi
} > "$RUN_DIR/meta/sessionInfo.txt"

if ls logs/*.log >/dev/null 2>&1; then
  cp -f logs/*.log "$RUN_DIR/logs/" || true
fi
if [[ -f logs/cloud_full_driver.log ]]; then
  cp -f logs/cloud_full_driver.log "$RUN_DIR/logs/" || true
fi

hash_cmd=""
if command -v sha256sum >/dev/null 2>&1; then
  hash_cmd="sha256sum"
elif command -v shasum >/dev/null 2>&1; then
  hash_cmd="shasum -a 256"
fi

if [[ -n "$hash_cmd" ]]; then
  (
    cd "$RUN_DIR"
    find . -type f ! -name "manifest.sha256" -print0 \
      | sort -z \
      | xargs -0 $hash_cmd
  ) > "$RUN_DIR/manifest.sha256" 2>/dev/null || true

  # Also record checksums for the *pipeline artifacts* in-place (they are committed
  # in the repo at the same git SHA), without duplicating large files into the audit bundle.
  (
    cd "$REPO_ROOT"
    {
      find results/figures results/scpagwas results/annotation results/causal results/metadata plots/publication docs \
        -type f \
        \( \
          -path "docs/DECISION_LOG.md" -o \
          -path "docs/ENVIRONMENT.md" -o \
          -path "docs/INSTALL_MICROMAMBA_SEURAT4.md" -o \
          -path "results/figures/*" -o \
          -path "results/scpagwas/*" -o \
          -path "results/annotation/*" -o \
          -path "results/causal/*" -o \
          -path "results/metadata/*" -o \
          -path "plots/publication/*" \
        \) \
        -print0 2>/dev/null \
        | sort -z \
        | xargs -0 $hash_cmd
    } 2>/dev/null
  ) > "$RUN_DIR/meta/artifacts.sha256" 2>/dev/null || true
fi

echo "OK: audit bundle written to $RUN_DIR"
