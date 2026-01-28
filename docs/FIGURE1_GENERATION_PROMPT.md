---
title: "Figure 1 schematic generation prompt"
date: "2026-01-28"
status: "needs regeneration (workflow updated)"
target_image: "plots/publication/png/figure1.jpeg"
notes: "Prompt archived for transparency; keep numbers consistent with the result tables."
---

# Prompt

Create a clean scientific workflow schematic for a DKD genetics × kidney single-cell integration study. Balance visuals and short labels (no long paragraphs).

Hard constraints:
- Do NOT include “Figure 1” / “Fig. 1” / any figure numbering inside the artwork.
- No checkmarks/ticks/emojis/badges.
- White background; clean flat 2D look; no sci-fi UI; no glass; no heavy gradients.

Main title inside artwork:
Study overview and evidence tiers

Layout:
Three modules left-to-right with arrows:
Module 1 (Inputs) → Module 2 (Primary evidence) → Module 3 (Secondary/supporting)
Use solid arrows for primary evidence and dashed arrows for secondary/supporting. Include a small legend (solid vs dashed).

Module 1 (Inputs):
- Manhattan plot icon + label: “GWAS: GCST005881 (DKD in diabetes; GRCh37)”
- Two kidney atlas cards: “GSE131882 (primary)” and “GSE195460 (replication)”

Module 2 (Primary evidence):
Header: “Single-cell integration + scPagwas localization”
Show a Seurat integration icon feeding into a UMAP thumbnail and then a small heatmap/regression icon labeled “scPagwas”.
UMAP callouts (small, one line each):
- Podocyte (C12): FDR=5.48e−4
- Proximal tubule S3-like (C3): FDR=5.48e−4
- Collecting duct principal-like (C0): FDR=0.00125
- Injury-associated epithelial (C13): FDR=0.00161
QC strip (single line):
N=20,243; exclude 4,043/92,186 (4.39%)

Module 3 (Secondary/supporting) as four compact sub-panels:
A) “Cross-GWAS robustness” with a small dot-matrix icon labeled: DKD meta (GENIE); CKD (CKDGen); eGFR (CKDGen). Small annotation: “rank concordance: 0.42, 0.41, −0.04”.
B) “Screening MR (OpenGWAS)” with a small forest icon and gene labels: PAX8-AS1; JAG1; H3-3A; ID3 (tag: screening).
C) “Colocalization (coloc.abf)” with an overlap icon and “max PP.H4≈0.034” (tag: consistency).
D) “CellChat (exploratory)” with a network icon and pathway tags: COLLAGEN/LAMININ; VEGF (tag: exploratory).

Color rules:
Mostly grayscale; use muted orange accents for the primary DKD-localized contexts, muted blue accents for eGFR-related elements, and neutral outlines elsewhere.

Do not introduce any additional numbers or claims beyond those specified above.
Do not claim MR/coloc proves causality.

