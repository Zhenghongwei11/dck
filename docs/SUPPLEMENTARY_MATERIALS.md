# Supplementary Materials (Draft)

This document enumerates the supplementary tables and notes referenced by the manuscript and maps each item to the repository “source-of-record” outputs under `results/`.

## Supplementary Tables

**Table S1. scPagwas cluster-level enrichment (primary atlas; main vs repeat).**  
Source tables: `results/scpagwas/main/merged_celltype_pvalue.csv`, `results/scpagwas/repeat_run/merged_celltype_pvalue.csv`.

**Table S2. scPagwas cell-level scores (primary atlas; main vs repeat).**  
Source tables: `results/scpagwas/main/cell_scores.tsv`, `results/scpagwas/repeat_run/cell_scores.tsv`.

**Table S3. scPagwas cluster-level enrichment (replication atlas; GSE195460).**  
Source table: `results/scpagwas/gse195460/rep_gse195460_Merged_celltype_pvalue.csv`.

**Table S13. scPagwas subsampling stability (bootstrap).**  
Source tables: `results/scpagwas/robustness/bootstrap_cluster_stability.tsv`, `results/scpagwas/robustness/bootstrap_replicate_correlations.tsv`.

**Table S4. MR screening validation (IVW/WM/Egger; heterogeneity and pleiotropy diagnostics).**  
Source table: `results/causal/mr_validation.tsv`.

**Table S5. MR estimates (IVW only; full list).**  
Source table: `results/causal/mr_summary_ivw.tsv`.

**Table S6. MR estimates (all methods; full list).**  
Source table: `results/causal/mr_summary_all_methods.tsv`.

**Table S7. MR dataset provenance and instrument strength summary.**  
Source tables: `results/causal/mr_dataset_provenance.tsv`, `results/causal/mr_instrument_strength.tsv`.

**Table S8. Feasibility and SNP overlap analysis for colocalization (attempt-level diagnostics).**  
Source table: `results/causal/coloc/coloc_attempts.tsv`.

**Table S9. Colocalization summary (strict overlap threshold; `EZHU_COLOC_MIN_OVERLAP=50`).**  
Source table: `results/causal/coloc_summary_overlap50.tsv`.

**Table S10. Colocalization summary (exploratory overlap threshold; `EZHU_COLOC_MIN_OVERLAP=25`).**  
Source table: `results/causal/coloc_summary_overlap25.tsv`.

Caption template for Tables S8–S10:  
“Colocalization was performed in loci defined as ±500 kb around the strongest exposure instrument when available (fallback: gene coordinate ±500 kb), using OpenGWAS regional associations on GRCh37. Exposure eQTLs were sourced from OpenGWAS and are not kidney-injury-state-specific; therefore, heterogeneous colocalization support across loci should be interpreted in the context of tissue/measurement mismatch and exposure-definition limits. Tables S9 and S10 differ only by the minimum overlap threshold (`EZHU_COLOC_MIN_OVERLAP=50` vs `25`). SNP-level locus tables are saved under `results/causal/coloc/locus_*.tsv.gz` and serve as the source of record for auditability.”

**Table S15. scPagwas multi-GWAS control comparison (Main vs Control GWAS).**  
Source tables: `results/scpagwas/gwas_controls/controls_summary.tsv`, `results/scpagwas/gwas_controls/controls_top_contexts.tsv`.

## Supplementary Figures (diagnostics)

**Figure S6. scPagwas subsampling stability (bootstrap within sample).**  
Outputs: `plots/publication/figureS6_scpagwas_subsampling_stability.pdf` and `plots/publication/png/figureS6_scpagwas_subsampling_stability.png`.

## Supplementary Note 1: Interpreting colocalization under tissue mismatch

The optional colocalization layer is implemented as a *locus-level consistency check* between an exposure eQTL resource and renal outcome GWAS. In this repository, loci are defined as ±500 kb around the strongest exposure instrument when available (fallback: gene coordinate ±500 kb; see `results/causal/coloc/coloc_attempts.tsv`), and the overlap threshold is varied across strict vs exploratory summaries (`EZHU_COLOC_MIN_OVERLAP=50` vs `25`).

Exposure eQTLs are sourced from OpenGWAS and are not kidney-injury-state-specific (many widely used OpenGWAS eQTL resources are systemic and/or blood-derived). Therefore, a low posterior probability for a shared causal variant (PP.H4) should be interpreted as **non-supportive evidence under the tested exposure definition**, rather than as a definitive refutation of kidney-cell-type-specific regulatory mechanisms.

Accordingly, the manuscript reports colocalization as a supportive analysis and explicitly separates it from the primary genetics × kidney cell-context localization results.

For audit stability, the “source of record” for each tested locus is the saved SNP-level locus table under `results/causal/coloc/locus_*.tsv.gz`. If any live OpenGWAS API query changes over time (or fails transiently), the repository can deterministically restore the same summary tables from these saved locus tables via `make coloc-rebuild`.
