# Preflight Checks (GWAS / scPagwas)

This repository uses a small set of preflight checks to eliminate common integration failures before running scPagwas.

## Pitfall 1 — GWAS genome build mismatch (hg37 vs hg38)

Policy:
- The GWAS row in `data/manifest.tsv` must declare `build=hg37` or `build=hg38` in the `notes` column.
- The scPagwas run must use the matching `block_annotation` (hg37 vs hg38).

## Pitfall 2 — GWAS column schema mismatch

scPagwas examples expect a schema equivalent to:

`chrom pos rsid se beta maf`

The preflight script will:
- detect common column aliases
- derive `beta=log(OR)` if only `OR` is available
- enforce basic QC (SE>0, MAF in (0,0.5])
- output a harmonized file at `data/processed/gwas_harmonized.tsv.gz`

Run:
- `make preflight`

## Pitfall 3 — LD panel vs GWAS ancestry mismatch

Policy:
- The GWAS row must declare `ancestry=...` (e.g., EUR/EAS/AFR) in `notes`.
- If you use an LD panel, its ancestry must be recorded and discussed when it differs from the GWAS ancestry.

At minimum, record the decision in `docs/DECISION_LOG.md`.

