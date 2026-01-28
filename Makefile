SHELL := /bin/bash

RSCRIPT ?= Rscript
AUDIT_LEVEL ?= full
N_CORES ?= 1
FIGURES ?=
FIGURE_FORMAT ?= pdf,png
SEED ?= 20251227
GWAS_ID ?= GCST005881
GWAS_HARMONIZED ?= data/processed/gwas_harmonized__$(GWAS_ID).tsv.gz
DOWNLOAD_EXCLUDE_IDS ?= GCST008058,GCST008064
SCPAGWAS_PRIMARY_MAIN_PREFIX ?= gse131882_main
SCPAGWAS_PRIMARY_REPEAT_PREFIX ?= gse131882_repeat
SCPAGWAS_REP_PREFIX ?= rep_gse195460

.PHONY: help setup data preprocess analysis figures audit audit-bundle reproduce
.PHONY: figure2 figure3 figure4 figure5 figure6 figure7
.PHONY: singlecell-universe
.PHONY: causal
.PHONY: mr-setup
.PHONY: mr-provenance
.PHONY: coloc-setup
.PHONY: coloc
.PHONY: coloc-rebuild
.PHONY: reproduce-local
.PHONY: stage00
.PHONY: preflight
.PHONY: ccc
.PHONY: replication-concordance
.PHONY: replication-metrics
.PHONY: anchors-qc
.PHONY: anchors-umap
.PHONY: unpack-gse131882
.PHONY: preprocess-gse131882
.PHONY: preprocess-gse195460
.PHONY: analysis-primary
.PHONY: analysis-replication
.PHONY: install-scpagwas
.PHONY: annotation

help:
	@echo "Targets:"
	@echo "  make setup       - restore R packages via renv"
	@echo "  make stage00     - environment + run metadata"
	@echo "  make data        - acquire/validate external data (manifest-driven)"
	@echo "  make preflight   - validate GWAS format/build/ancestry (set GWAS_ID=... to choose)"
	@echo "  make preprocess  - preprocessing stage"
	@echo "  make analysis    - analysis stage"
	@echo "    - Optional: N_CORES=8 to set EZHU_N_CORES"
	@echo "  make figures     - publication figures stage (Figures 2–7; Figure 1 is external)"
	@echo "    - Optional: FIGURES=2,3,4 (default: all available)"
	@echo "    - Optional: FIGURE_FORMAT=pdf,png (default: pdf,png)"
	@echo "  make ccc         - cell-cell communication (CellChat; requires Seurat RDS)"
	@echo "  make replication-concordance - summarize primary vs replication matching"
	@echo "  make replication-metrics - compute replication concordance metrics"
	@echo "  make anchors-qc  - export QC anchor tables from a Seurat RDS"
	@echo "  make causal      - causal inference stage (exports candidates; MR/SMR scaffold)"
	@echo "  make mr-setup    - install MR deps (TwoSampleMR/ieugwasr)"
	@echo "  make mr-provenance - write MR dataset/provenance tables"
	@echo "  make coloc-setup - install colocalization deps (coloc + ieugwasr; optional gwasvcf)"
	@echo "  make coloc       - optional colocalization follow-up for top MR hits"
	@echo "  make coloc-rebuild - rebuild coloc summaries from saved locus tables"
	@echo "  make audit       - run audits (AUDIT_LEVEL=full|ci)"
	@echo "  make reproduce   - run setup->data->preprocess->analysis->figures->audit"

setup:
	$(RSCRIPT) -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv"); renv::restore(prompt = FALSE)'

stage00:
	$(RSCRIPT) scripts/00_setup.R

data:
	EZHU_DOWNLOAD_EXCLUDE_IDS=$(DOWNLOAD_EXCLUDE_IDS) $(RSCRIPT) scripts/01_download_data.R

unpack-gse131882:
	EZHU_SEED=$(SEED) $(RSCRIPT) scripts/01b_unpack_gse131882.R

preflight:
	EZHU_GWAS_MANIFEST_ID=$(GWAS_ID) EZHU_GWAS_HARMONIZED_OUT=$(GWAS_HARMONIZED) $(RSCRIPT) scripts/01_preflight_gwas.R

preprocess-gse131882: unpack-gse131882
	EZHU_SEED=$(SEED) EZHU_SINGLECELL_TAG=gse131882 EZHU_SINGLECELL_INPUT_DIR=data/raw/singlecell/GSE131882/files EZHU_SINGLECELL_PROCESSED=data/processed/singlecell_gse131882_seurat.rds EZHU_N_CORES=$(N_CORES) $(RSCRIPT) scripts/02_preprocess.R

preprocess-gse195460:
	EZHU_SEED=$(SEED) EZHU_SINGLECELL_TAG=gse195460 EZHU_SINGLECELL_INPUT_DIR=data/raw/singlecell/GSE195460/h5 EZHU_SINGLECELL_PROCESSED=data/processed/singlecell_gse195460_seurat_v4.rds EZHU_N_CORES=$(N_CORES) $(RSCRIPT) scripts/02_preprocess.R

preprocess: preprocess-gse131882 preprocess-gse195460

analysis-primary:
	EZHU_GWAS_MANIFEST_ID=$(GWAS_ID) EZHU_GWAS_HARMONIZED=$(GWAS_HARMONIZED) EZHU_SEED=$(SEED) EZHU_SINGLECELL_TAG=gse131882 EZHU_SINGLECELL_PROCESSED=data/processed/singlecell_gse131882_seurat.rds EZHU_SCPAGWAS_DIR=results/scpagwas/main EZHU_SCPAGWAS_PREFIX=$(SCPAGWAS_PRIMARY_MAIN_PREFIX) EZHU_N_CORES=$(N_CORES) $(RSCRIPT) scripts/03_analysis.R
	EZHU_GWAS_MANIFEST_ID=$(GWAS_ID) EZHU_GWAS_HARMONIZED=$(GWAS_HARMONIZED) EZHU_SEED=$$(($(SEED)+1)) EZHU_SINGLECELL_TAG=gse131882 EZHU_SINGLECELL_PROCESSED=data/processed/singlecell_gse131882_seurat.rds EZHU_SCPAGWAS_DIR=results/scpagwas/repeat_run EZHU_SCPAGWAS_PREFIX=$(SCPAGWAS_PRIMARY_REPEAT_PREFIX) EZHU_N_CORES=$(N_CORES) $(RSCRIPT) scripts/03_analysis.R

analysis-replication:
	EZHU_GWAS_MANIFEST_ID=$(GWAS_ID) EZHU_GWAS_HARMONIZED=$(GWAS_HARMONIZED) EZHU_SEED=$(SEED) EZHU_SINGLECELL_TAG=gse195460 EZHU_SINGLECELL_PROCESSED=data/processed/singlecell_gse195460_seurat_v4.rds EZHU_SCPAGWAS_DIR=results/scpagwas/gse195460 EZHU_SCPAGWAS_PREFIX=$(SCPAGWAS_REP_PREFIX) EZHU_N_CORES=$(N_CORES) $(RSCRIPT) scripts/03_analysis.R

analysis: install-scpagwas preflight analysis-primary analysis-replication

install-scpagwas:
	$(RSCRIPT) scripts/00_install_scpagwas.R

figures:
	EZHU_FIGURES=$(FIGURES) EZHU_FIGURE_FORMAT=$(FIGURE_FORMAT) $(RSCRIPT) scripts/07_figures.R

figure2:
	$(MAKE) figures FIGURES=2

figure3:
	$(MAKE) figures FIGURES=3

figure4:
	$(MAKE) figures FIGURES=4

figure5:
	$(MAKE) figures FIGURES=5

figure6:
	$(MAKE) figures FIGURES=6

figure7:
	$(MAKE) figures FIGURES=7

scpagwas-gwas-controls:
	./scripts/03_analysis_multi_gwas.sh

scpagwas-gwas-controls-summary:
	$(RSCRIPT) scripts/11_scpagwas_gwas_controls_summary.R

ccc-diff:
	$(RSCRIPT) scripts/03b_cellchat_case_control_diff.R

ccc:
	$(RSCRIPT) scripts/03b_cellchat.R

anchors-qc:
	EZHU_SEED=$(SEED) EZHU_SEURAT_RDS=data/processed/singlecell_gse131882_seurat.rds $(RSCRIPT) scripts/02d_export_qc_anchors.R

anchors-umap:
	EZHU_SEED=$(SEED) EZHU_SEURAT_RDS=data/processed/singlecell_gse131882_seurat.rds $(RSCRIPT) scripts/02c_export_umap_anchors.R

annotation:
	EZHU_SEED=$(SEED) EZHU_SCPAGWAS_RDS=data/processed/singlecell_gse131882_scpagwas_$(SCPAGWAS_PRIMARY_REPEAT_PREFIX).rds $(RSCRIPT) scripts/04_annotation.R

replication-concordance:
	$(RSCRIPT) scripts/03c_replication_concordance.R

replication-metrics:
	$(RSCRIPT) scripts/10_replication_metrics.R

singlecell-universe:
	EZHU_DISABLE_RENV=1 $(RSCRIPT) scripts/02b_export_gene_universe_gse131882.R

causal:
	$(RSCRIPT) scripts/06_causal_inference.R

mr-setup:
	$(RSCRIPT) scripts/00_install_mr_deps.R

mr-provenance:
	$(RSCRIPT) scripts/09_mr_provenance.R

coloc-setup:
	$(RSCRIPT) scripts/00_install_coloc_deps.R

coloc: coloc-setup
	$(RSCRIPT) scripts/06b_colocalization.R

coloc-rebuild:
	$(RSCRIPT) scripts/10_coloc_rebuild_from_locus_tables.R

audit:
	EZHU_AUDIT_LEVEL=$(AUDIT_LEVEL) $(RSCRIPT) scripts/audit/run_audit.R

audit-bundle:
	./scripts/audit/bundle_run.sh 0

reproduce: setup stage00 data preprocess anchors-qc anchors-umap analysis annotation figures audit

reproduce-local: figures
	$(MAKE) audit AUDIT_LEVEL=ci
