#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

AUDIT_LEVEL="${AUDIT_LEVEL:-full}"
make reproduce AUDIT_LEVEL="$AUDIT_LEVEL"

