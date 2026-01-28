#!/usr/bin/env bash
# scripts/03_analysis_multi_gwas.sh
#
# Purpose:
# - Run scPagwas against one or more *control GWAS* IDs (incremental runs).
# - Writes outputs under: results/scpagwas/gwas_controls/<GWAS_ID>/
#
# Safety:
# - Requires explicit EZHU_GWAS_IDS to avoid accidentally re-running the primary GWAS.
# - Writes harmonized GWAS to a per-ID file (no shared default output path).
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

MANIFEST="${EZHU_MANIFEST:-data/manifest.tsv}"
RSCRIPT="${RSCRIPT:-./scripts/run_r_mamba.sh}"
N_CORES="${EZHU_N_CORES:-1}"

SINGLECELL_TAG="${EZHU_SINGLECELL_TAG:-gse131882}"
SINGLECELL_PROCESSED="${EZHU_SINGLECELL_PROCESSED:-data/processed/singlecell_gse131882_seurat.rds}"

# Default: run both celltype + single-cell layers for control GWAS.
# In practice, some scPagwas versions still enter the single-cell scoring stage even when `singlecell=FALSE`,
# which can yield confusing failures late in the run. Override with EZHU_SCPAGWAS_CONTROLS_SINGLECELL=0 if desired.
SCPAGWAS_CONTROLS_SINGLECELL="${EZHU_SCPAGWAS_CONTROLS_SINGLECELL:-1}"
SCPAGWAS_CONTROLS_CELLTYPE="${EZHU_SCPAGWAS_CONTROLS_CELLTYPE:-1}"
PRINT_HASHES="${EZHU_PRINT_HASHES:-1}"
SMOKE="${EZHU_SCPAGWAS_SMOKE:-0}"
# For smoke runs, limit the harmonized GWAS size to reduce SnpToGene/Link_pathway_blocks runtime.
if [[ "$SMOKE" == "1" && -z "${EZHU_GWAS_PREFLIGHT_MAX_VARIANTS:-}" ]]; then
  export EZHU_GWAS_PREFLIGHT_MAX_VARIANTS="200000"
fi

if [[ ! -f "$MANIFEST" ]]; then
  echo "ERROR: manifest not found: $MANIFEST" >&2
  exit 2
fi
if [[ ! -x "$RSCRIPT" ]] && ! command -v "$RSCRIPT" >/dev/null 2>&1; then
  echo "ERROR: RSCRIPT not found/executable: $RSCRIPT" >&2
  exit 2
fi

available_ids="$(awk -F'\t' 'NR>1 && $2=="gwas" {print $1}' "$MANIFEST" | tr '\n' ' ' | sed -e 's/[[:space:]]\\+/ /g' -e 's/ $//')"
if [[ -z "${EZHU_GWAS_IDS:-}" ]]; then
  echo "ERROR: EZHU_GWAS_IDS is required to run control GWAS analyses." >&2
  echo "Example: EZHU_GWAS_IDS=GCST008058,GCST008064 $0" >&2
  echo "Available GWAS IDs in $MANIFEST: $available_ids" >&2
  exit 2
fi

hash_cmd=""
if command -v sha256sum >/dev/null 2>&1; then
  hash_cmd="sha256sum"
elif command -v shasum >/dev/null 2>&1; then
  hash_cmd="shasum -a 256"
fi

main_merged="${EZHU_MAIN_MERGED_CELLTYPE_PVALUE:-results/scpagwas/repeat_run/merged_celltype_pvalue.csv}"
if [[ ! -f "$main_merged" && -f "results/scpagwas/main/merged_celltype_pvalue.csv" ]]; then
  main_merged="results/scpagwas/main/merged_celltype_pvalue.csv"
fi

SKIP_EXISTING="${EZHU_SKIP_EXISTING:-1}"

IFS=',' read -ra IDS <<< "${EZHU_GWAS_IDS}"
for raw_id in "${IDS[@]}"; do
  ID="$(echo "$raw_id" | xargs)"
  [[ -z "$ID" ]] && continue

  echo "== scPagwas control GWAS: $ID =="

  HARMONIZED_FILE="data/processed/gwas_harmonized__${ID}.tsv.gz"
  OUTPUT_DIR="results/scpagwas/gwas_controls/${ID}"
  PREFIX="${SINGLECELL_TAG}_${ID,,}"

  mkdir -p "data/processed" "results/metadata" "$OUTPUT_DIR"

  merged_out="${OUTPUT_DIR}/merged_celltype_pvalue.csv"
  if [[ ! -f "$merged_out" ]]; then
    merged_out="$(ls -1 "${OUTPUT_DIR}"/*Merged_celltype_pvalue.csv "${OUTPUT_DIR}"/*merged_celltype_pvalue.csv 2>/dev/null | head -n 1 || true)"
  fi

  if [[ "$SKIP_EXISTING" == "1" && -n "${merged_out}" && -f "${merged_out}" ]]; then
    echo "Found existing merged output: $merged_out"
    if [[ -n "$hash_cmd" && -f "$main_merged" ]]; then
      main_hash="$($hash_cmd "$main_merged" | awk '{print $1}')"
      ctrl_hash="$($hash_cmd "$merged_out" | awk '{print $1}')"
      if [[ "$PRINT_HASHES" == "1" ]]; then
        echo "  sha(primary)=$main_hash"
        echo "  sha(control )=$ctrl_hash"
      fi
      if [[ "$main_hash" == "$ctrl_hash" && "${EZHU_ALLOW_IDENTICAL_CONTROLS:-0}" != "1" ]]; then
        echo "ERROR: existing control output is identical to primary GWAS; refusing to proceed." >&2
        echo "Delete ${OUTPUT_DIR} and rerun, or set EZHU_ALLOW_IDENTICAL_CONTROLS=1 (not recommended)." >&2
        exit 3
      fi
    fi
    echo "Skipping $ID (EZHU_SKIP_EXISTING=1)."
    continue
  fi

  echo "Preflight GWAS -> $HARMONIZED_FILE"
  export EZHU_GWAS_MANIFEST_ID="$ID"
  export EZHU_GWAS_HARMONIZED_OUT="$HARMONIZED_FILE"
  export EZHU_GWAS_PREFLIGHT_QC_OUT="results/metadata/01_preflight_gwas__qc__${ID}.tsv"
  echo "  EZHU_GWAS_MANIFEST_ID=$EZHU_GWAS_MANIFEST_ID"
  echo "  EZHU_GWAS_HARMONIZED_OUT=$EZHU_GWAS_HARMONIZED_OUT"
  if [[ -n "${EZHU_GWAS_PREFLIGHT_MAX_VARIANTS:-}" ]]; then
    echo "  EZHU_GWAS_PREFLIGHT_MAX_VARIANTS=$EZHU_GWAS_PREFLIGHT_MAX_VARIANTS"
  fi
  _run_preflight() {
    "$RSCRIPT" scripts/01_preflight_gwas.R
  }
  _run_preflight

  if [[ ! -f "$HARMONIZED_FILE" ]]; then
    echo "ERROR: harmonized GWAS file not found after preflight: $HARMONIZED_FILE" >&2
    exit 1
  fi
  if [[ ! -s "results/metadata/01_preflight_gwas__qc__${ID}.tsv" ]]; then
    echo "ERROR: missing/empty preflight QC output: results/metadata/01_preflight_gwas__qc__${ID}.tsv" >&2
    exit 1
  fi

  # Fail fast: scPagwas expects rsIDs that match the LD reference panel.
  # NOTE: don't use a short-circuiting gzip|head/awk pipeline here; it can raise SIGPIPE (exit 141)
  # when the consumer exits early, and our script runs with `set -euo pipefail`.
  rs_count=""
  if command -v python3 >/dev/null 2>&1; then
    rs_count="$(python3 - "$HARMONIZED_FILE" <<'PY'
import gzip, sys
path = sys.argv[1]
c = 0
with gzip.open(path, "rt") as f:
    header = f.readline()
    for i, line in enumerate(f, start=1):
        if i > 200:
            break
        cols = line.rstrip("\n").split("\t")
        if len(cols) >= 3 and cols[2].startswith("rs"):
            c += 1
print(c)
PY
)"
  elif command -v gzip >/dev/null 2>&1; then
    # Fallback: avoid SIGPIPE by scanning the full stream (slower, but safe).
    rs_count="$(gzip -dc "$HARMONIZED_FILE" | awk -F'\t' 'NR>1 && NR<=201 && $3 ~ /^rs/ {c++} END{print c+0}')"
  fi

  if [[ -n "${rs_count:-}" && "${rs_count:-0}" == "0" ]]; then
    echo "ERROR: harmonized GWAS rsid column does not contain any 'rs...' IDs in the first 200 variants: $HARMONIZED_FILE" >&2
    echo "Fix: ensure the raw GWAS has an rsID column (often 'rsid' or 'rs_id'), or choose a different GWAS." >&2
    echo "If you *must* fill missing IDs (not recommended), set EZHU_PREFLIGHT_FILL_MISSING_RSID=1." >&2
    exit 1
  fi

  rows_out="$(awk -F'\t' '$1=="rows_out"{print $2; exit}' "results/metadata/01_preflight_gwas__qc__${ID}.tsv" 2>/dev/null || true)"
  rows_out="${rows_out:-}"
  if [[ -z "$rows_out" ]]; then rows_out="0"; fi
  echo "  preflight rows_out=$rows_out"

  # In SMOKE mode we may cap the number of harmonized variants. For smaller GWAS,
  # the first chunk can occasionally yield 0 passing variants (e.g., missing AF/SE),
  # which then crashes scPagwas downstream. Fail fast and retry once without the cap.
  if [[ "$SMOKE" == "1" && "$rows_out" == "0" && -n "${EZHU_GWAS_PREFLIGHT_MAX_VARIANTS:-}" ]]; then
    echo "WARN: preflight produced 0 variants under EZHU_GWAS_PREFLIGHT_MAX_VARIANTS=$EZHU_GWAS_PREFLIGHT_MAX_VARIANTS" >&2
    echo "      Retrying preflight once without EZHU_GWAS_PREFLIGHT_MAX_VARIANTS..." >&2
    unset EZHU_GWAS_PREFLIGHT_MAX_VARIANTS
    _run_preflight

    if [[ ! -s "results/metadata/01_preflight_gwas__qc__${ID}.tsv" ]]; then
      echo "ERROR: missing/empty preflight QC output after retry: results/metadata/01_preflight_gwas__qc__${ID}.tsv" >&2
      exit 1
    fi
    rows_out="$(awk -F'\t' '$1=="rows_out"{print $2; exit}' "results/metadata/01_preflight_gwas__qc__${ID}.tsv" 2>/dev/null || true)"
    rows_out="${rows_out:-}"
    if [[ -z "$rows_out" ]]; then rows_out="0"; fi
    echo "  preflight rows_out(after retry)=$rows_out"
    if [[ "$rows_out" == "0" ]]; then
      echo "ERROR: preflight still produced 0 passing variants; scPagwas cannot run." >&2
      echo "Fix: verify the raw GWAS file schema/contents and the AF/SE columns in data/raw/gwas/${ID}*" >&2
      exit 1
    fi
  fi

  echo "Run scPagwas -> $OUTPUT_DIR (prefix=$PREFIX)"
  export EZHU_GWAS_HARMONIZED="$HARMONIZED_FILE"
  export EZHU_GWAS_MANIFEST_ID="$ID"
  export EZHU_SINGLECELL_TAG="$SINGLECELL_TAG"
  export EZHU_SINGLECELL_PROCESSED="$SINGLECELL_PROCESSED"
  export EZHU_SCPAGWAS_DIR="$OUTPUT_DIR"
  export EZHU_SCPAGWAS_PREFIX="$PREFIX"
  export EZHU_N_CORES="$N_CORES"
  export EZHU_SCPAGWAS_SINGLECELL="$SCPAGWAS_CONTROLS_SINGLECELL"
  export EZHU_SCPAGWAS_CELLTYPE="$SCPAGWAS_CONTROLS_CELLTYPE"
  "$RSCRIPT" scripts/03_analysis.R
  # Re-scan after the run because the output file is created by R.
  merged_out="${OUTPUT_DIR}/merged_celltype_pvalue.csv"
  if [[ ! -f "$merged_out" ]]; then
    merged_out="$(ls -1 "${OUTPUT_DIR}"/*Merged_celltype_pvalue.csv "${OUTPUT_DIR}"/*merged_celltype_pvalue.csv 2>/dev/null | head -n 1 || true)"
  fi
  if [[ -z "${merged_out}" || ! -f "${merged_out}" ]]; then
    echo "ERROR: could not find merged celltype p-value output under: $OUTPUT_DIR" >&2
    exit 1
  fi

  if [[ -n "$hash_cmd" && -f "$main_merged" ]]; then
    main_hash="$($hash_cmd "$main_merged" | awk '{print $1}')"
    ctrl_hash="$($hash_cmd "$merged_out" | awk '{print $1}')"
    if [[ "$PRINT_HASHES" == "1" ]]; then
      echo "  sha(primary)=$main_hash"
      echo "  sha(control )=$ctrl_hash"
    fi
    if [[ "$main_hash" == "$ctrl_hash" && "${EZHU_ALLOW_IDENTICAL_CONTROLS:-0}" != "1" ]]; then
      echo "ERROR: control GWAS produced an identical merged_celltype_pvalue to the primary GWAS." >&2
      echo "  primary: $main_merged ($main_hash)" >&2
      echo "  control: $merged_out ($ctrl_hash)" >&2
      echo "This usually means the control GWAS input did not take effect." >&2
      echo "Set EZHU_ALLOW_IDENTICAL_CONTROLS=1 to bypass (not recommended)." >&2
      exit 3
    fi
  fi

  echo "OK: $ID completed"
done
