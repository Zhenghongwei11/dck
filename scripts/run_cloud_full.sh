#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

mkdir -p logs

# Defaults tuned for stability on large-memory VMs.
export EZHU_SEED="${EZHU_SEED:-20251227}"
export EZHU_INTEGRATION_METHOD="${EZHU_INTEGRATION_METHOD:-seurat}"
export EZHU_N_CORES="${EZHU_N_CORES:-1}"
export EZHU_DOWNLOAD_TIMEOUT="${EZHU_DOWNLOAD_TIMEOUT:-7200}"

GIT_SHA="$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")"
UTC_TS="$(date -u +"%Y%m%dT%H%M%SZ" 2>/dev/null || date +"%Y%m%dT%H%M%S")"
export EZHU_RUN_ID="${EZHU_RUN_ID:-${UTC_TS}_${GIT_SHA}}"

if [[ "${EZHU_BUNDLE_ON_EXIT:-1}" == "1" ]] && [[ -x "./scripts/audit/bundle_run.sh" ]]; then
  _ezhu_bundle_on_exit() {
    local exit_code="$?"
    ./scripts/audit/bundle_run.sh "$exit_code" || true
  }
  trap _ezhu_bundle_on_exit EXIT
fi

RSCRIPT_BIN="${RSCRIPT_BIN:-}"
if [[ -z "${RSCRIPT_BIN}" ]]; then
  if [[ -x "./scripts/run_r_mamba.sh" ]]; then
    RSCRIPT_BIN="./scripts/run_r_mamba.sh"
  else
    RSCRIPT_BIN="Rscript"
  fi
fi

if [[ "${RSCRIPT_BIN}" == "./scripts/run_r_mamba.sh" ]]; then
  export EZHU_MAMBA_ENV="${EZHU_MAMBA_ENV:-ezhu-r-seurat4}"
fi

echo "[cloud-full] repo_root=$REPO_ROOT"
echo "[cloud-full] RSCRIPT_BIN=$RSCRIPT_BIN"
if [[ -n "${EZHU_MAMBA_ENV:-}" ]]; then
  echo "[cloud-full] EZHU_MAMBA_ENV=$EZHU_MAMBA_ENV"
fi
echo "[cloud-full] EZHU_RUN_ID=$EZHU_RUN_ID"
echo "[cloud-full] EZHU_SEED=$EZHU_SEED"
echo "[cloud-full] EZHU_INTEGRATION_METHOD=$EZHU_INTEGRATION_METHOD"
echo "[cloud-full] EZHU_N_CORES=$EZHU_N_CORES"
echo "[cloud-full] EZHU_DOWNLOAD_TIMEOUT=$EZHU_DOWNLOAD_TIMEOUT"

echo "[cloud-full] Matrix compatibility check (Seurat v4 requires mMatrix)"
if "$RSCRIPT_BIN" -e 'ok<-requireNamespace("Matrix",quietly=TRUE)&&requireNamespace("methods",quietly=TRUE); if(!ok) quit(status=1); suppressPackageStartupMessages(library(Matrix)); v<-as.character(utils::packageVersion("Matrix")); m_ok<-methods::isClass("mMatrix") && !is.null(methods::getClassDef("mMatrix")) && methods::getClassDef("mMatrix")@package=="Matrix"; ok2<-utils::compareVersion(v,"1.7-0")<0 && m_ok; quit(status=ifelse(ok2,0,3))' >/dev/null 2>&1; then
  echo "[cloud-full] Matrix ok" | tee logs/matrix_check.log
else
  echo "[cloud-full] ERROR: Matrix too new (>=1.7) or mMatrix missing; Seurat v4 Graph objects will fail." | tee logs/matrix_check.log
  echo "[cloud-full] Fix: recreate/update micromamba env from env/ezhu-r-seurat4.yml (pins a known-good r-matrix build)." | tee -a logs/matrix_check.log
  exit 3
fi

echo "[cloud-full] runtime deps (Seurat + scPagwas + hdf5r check)"
if "$RSCRIPT_BIN" -e 'ok<-requireNamespace("Seurat", quietly=TRUE)&&requireNamespace("SeuratObject", quietly=TRUE)&&requireNamespace("scPagwas", quietly=TRUE)&&requireNamespace("AnnotationDbi", quietly=TRUE)&&requireNamespace("org.Hs.eg.db", quietly=TRUE)&&requireNamespace("hdf5r", quietly=TRUE); if(!ok) quit(status=1); so<-as.character(utils::packageVersion("SeuratObject")); s<-as.character(utils::packageVersion("Seurat")); ok2<-utils::compareVersion(so,"5.0.0")<0 && utils::compareVersion(s,"5.0.0")<0; quit(status=ifelse(ok2,0,2))' >/dev/null 2>&1; then
  echo "[cloud-full] runtime deps ok" | tee logs/install_runtime_deps.log
else
  # If Seurat/SeuratObject are too new (v5+), scPagwas will fail at runtime.
  if "$RSCRIPT_BIN" -e 'if(!requireNamespace("SeuratObject",quietly=TRUE)||!requireNamespace("Seurat",quietly=TRUE)) quit(status=1); so<-as.character(utils::packageVersion("SeuratObject")); s<-as.character(utils::packageVersion("Seurat")); quit(status=ifelse(utils::compareVersion(so,"5.0.0")>=0||utils::compareVersion(s,"5.0.0")>=0,2,0))' >/dev/null 2>&1; then
    true
  else
    echo "[cloud-full] ERROR: Seurat v5 detected; scPagwas requires Seurat/SeuratObject v4 (<5.0.0)." | tee logs/install_runtime_deps.log
    echo "[cloud-full] Fix: recreate micromamba env from env/ezhu-r-seurat4.yml (clean env) then rerun." | tee -a logs/install_runtime_deps.log
    exit 2
  fi

  echo "[cloud-full] installing scPagwas runtime deps via scripts/00_install_scpagwas.R" | tee -a logs/install_runtime_deps.log
  "$RSCRIPT_BIN" scripts/00_install_scpagwas.R >> logs/install_runtime_deps.log 2>&1
fi

echo "[cloud-full] audit (ci)"
make audit AUDIT_LEVEL=ci RSCRIPT="$RSCRIPT_BIN" > logs/audit_ci.log 2>&1

echo "[cloud-full] data download"
make data RSCRIPT="$RSCRIPT_BIN" > logs/data_download.log 2>&1

echo "[cloud-full] unpack gse131882"
make unpack-gse131882 RSCRIPT="$RSCRIPT_BIN" > logs/unpack_gse131882.log 2>&1

echo "[cloud-full] preflight"
make preflight RSCRIPT="$RSCRIPT_BIN" > logs/preflight.log 2>&1

echo "[cloud-full] preprocess (gse131882 + gse195460)"
make preprocess N_CORES="$EZHU_N_CORES" RSCRIPT="$RSCRIPT_BIN" > logs/preprocess.log 2>&1

echo "[cloud-full] anchors (qc + umap)"
make anchors-qc anchors-umap RSCRIPT="$RSCRIPT_BIN" > logs/anchors.log 2>&1

echo "[cloud-full] scPagwas (primary main + repeat + replication)"
make analysis N_CORES="$EZHU_N_CORES" RSCRIPT="$RSCRIPT_BIN" > logs/analysis.log 2>&1

echo "[cloud-full] annotation (markers)"
make annotation RSCRIPT="$RSCRIPT_BIN" > logs/annotation.log 2>&1

echo "[cloud-full] replication concordance"
make replication-concordance RSCRIPT="$RSCRIPT_BIN" > logs/replication_concordance.log 2>&1

echo "[cloud-full] tcm enrichment"
make tcm RSCRIPT="$RSCRIPT_BIN" > logs/tcm.log 2>&1

echo "[cloud-full] causal candidates"
make causal RSCRIPT="$RSCRIPT_BIN" > logs/causal.log 2>&1

echo "[cloud-full] colocalization (optional)"
if [[ "${EZHU_ENABLE_COLOC:-0}" == "1" ]]; then
  make coloc RSCRIPT="$RSCRIPT_BIN" > logs/coloc.log 2>&1
else
  echo "[cloud-full] skip coloc: set EZHU_ENABLE_COLOC=1 to enable" | tee logs/coloc.log
fi

echo "[cloud-full] figures"
make figures RSCRIPT="$RSCRIPT_BIN" > logs/figures.log 2>&1

echo "[cloud-full] CCC (CellChat) (optional)"
export EZHU_SEURAT_RDS="${EZHU_SEURAT_RDS:-data/processed/singlecell_gse131882_seurat.rds}"
if "$RSCRIPT_BIN" -e 'ok<-requireNamespace("CellChat", quietly=TRUE); quit(status=ifelse(ok,0,1))' >/dev/null 2>&1; then
  make ccc RSCRIPT="$RSCRIPT_BIN" > logs/ccc.log 2>&1
else
  echo "[cloud-full] skip ccc: CellChat not installed in this R environment" | tee logs/ccc.log
fi

echo "[cloud-full] audit (full)"
make audit AUDIT_LEVEL=full RSCRIPT="$RSCRIPT_BIN" > logs/audit_full.log 2>&1

echo "[cloud-full] done"
