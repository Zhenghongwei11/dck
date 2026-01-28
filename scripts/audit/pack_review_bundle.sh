#!/usr/bin/env bash
set -euo pipefail

file_size_bytes() {
  # macOS/BSD: stat -f%z; GNU: stat -c%s
  if stat -f%z "$1" >/dev/null 2>&1; then
    stat -f%z "$1"
    return 0
  fi
  if stat -c%s "$1" >/dev/null 2>&1; then
    stat -c%s "$1"
    return 0
  fi
  echo "ERROR: unable to determine file size for $1 (stat)" >&2
  return 1
}

file_sha256() {
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$1" | awk '{print $1}'
    return 0
  fi
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
    return 0
  fi
  echo "ERROR: neither shasum nor sha256sum is available" >&2
  return 1
}

root_dir() {
  local script_dir
  script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  # scripts/audit -> repo root
  (cd "${script_dir}/../.." && pwd)
}

ROOT="$(root_dir)"
OUT_DIR="${ROOT}/docs/review_bundle"
ZIP_PATH="${OUT_DIR}/ezhu_review_bundle.zip"
MANIFEST_PATH="${OUT_DIR}/FILE_MANIFEST.tsv"
BUILD_INFO_PATH="${OUT_DIR}/BUILD_INFO.txt"

mkdir -p "${OUT_DIR}"

tmp="$(mktemp -d)"
trap 'rm -rf "${tmp}"' EXIT

stage="${tmp}/ezhu_review_bundle"
mkdir -p "${stage}"

if command -v git >/dev/null 2>&1 && git -C "${ROOT}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  # Prefer a git snapshot so the bundle can be built even from a sparse checkout.
  (cd "${ROOT}" && git archive --format=tar HEAD) | tar -x -C "${stage}"
else
  rsync -a \
    --exclude '.git/' \
    --exclude '.DS_Store' \
    --exclude 'renv/library/' \
    --exclude 'data/raw/' \
    --exclude 'data/processed/' \
    --exclude 'logs/' \
    --exclude 'trash/' \
    --exclude 'docs/review_bundle/' \
    "${ROOT}/" \
    "${stage}/"
fi

# Never include a previously-generated review bundle inside itself.
rm -rf "${stage}/docs/review_bundle" || true

# Keep at most one audit run to keep the bundle size practical.
# Default: keep the lexicographically-latest run directory under docs/audit_runs.
audit_keep="${EZHU_REVIEW_BUNDLE_AUDIT_RUN:-}"
if [ -d "${stage}/docs/audit_runs" ]; then
  if [ -z "${audit_keep}" ]; then
    audit_keep="$(ls -1 "${stage}/docs/audit_runs" 2>/dev/null | sort | tail -n 1 || true)"
  fi

  if [ -n "${audit_keep}" ] && [ -d "${stage}/docs/audit_runs/${audit_keep}" ]; then
    find "${stage}/docs/audit_runs" -mindepth 1 -maxdepth 1 -type d ! -name "${audit_keep}" -exec rm -rf {} + || true
  else
    rm -rf "${stage}/docs/audit_runs" || true
  fi
fi

{
  echo "build_utc=$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  echo "platform=$(uname -a)"
  if command -v git >/dev/null 2>&1 && git -C "${ROOT}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "git_commit=$(git -C "${ROOT}" rev-parse HEAD 2>/dev/null || true)"
    echo "git_status_porcelain_count=$(git -C "${ROOT}" status --porcelain 2>/dev/null | wc -l | tr -d ' ')"
  fi
} > "${BUILD_INFO_PATH}"

# Write a SHA-256 manifest for every file inside the staged bundle.
{
  echo -e "sha256\tsize_bytes\trelpath"
  (cd "${stage}" && find . -type f -print0) | while IFS= read -r -d '' f; do
    rel="${f#./}"
    size="$(file_size_bytes "${stage}/${rel}")"
    sha="$(file_sha256 "${stage}/${rel}")"
    echo -e "${sha}\t${size}\t${rel}"
  done
} > "${MANIFEST_PATH}"

# Embed manifest + build info inside the archive as well.
mkdir -p "${stage}/docs/review_bundle"
cp -f "${MANIFEST_PATH}" "${stage}/docs/review_bundle/FILE_MANIFEST.tsv"
cp -f "${BUILD_INFO_PATH}" "${stage}/docs/review_bundle/BUILD_INFO.txt"

rm -f "${ZIP_PATH}"
(cd "${tmp}" && zip -qry "${ZIP_PATH}" "ezhu_review_bundle")

echo "Wrote: ${ZIP_PATH}"
echo "Wrote: ${MANIFEST_PATH}"
echo "Wrote: ${BUILD_INFO_PATH}"
