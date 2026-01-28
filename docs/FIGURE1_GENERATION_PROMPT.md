---
title: "Figure 1 generation prompt"
date: "2026-01-28"
model: "Gemini (image generation)"
status: "needs regeneration (workflow updated)"
target_image: "plots/publication/png/figure1.jpeg"
notes: "Prompt archived for transparency; keep numbers consistent with audited result tables."
---

# Prompt

Create a BioRender-style (vector) study-design / workflow schematic for a DKD genetics × kidney single-cell integration study.

Hard constraints:
- Do NOT include “Figure 1” / “Fig. 1” / any figure numbering inside the artwork.
- Do NOT include any meta terms like “Q1”, “journal-grade”, “prompt”, or “AI”.
- No checkmarks/ticks/emojis/badges.
- No watermark/branding (do not write “BioRender” or “Created with …” anywhere).
- Clean flat 2D look; no sci-fi UI; no glass; no heavy gradients.
- Do NOT place a big global title in the figure. (The figure title will be handled by the manuscript caption.)
- Keep labels short (no paragraphs); keep typography consistent.

Layout:
Three columns left-to-right with arrows:
Inputs → Primary localization (scPagwas) → Secondary follow-up (screening/consistency)
Use solid arrows for primary evidence and dashed arrows for secondary/supporting. Include a small legend (solid vs dashed).

Column A (Inputs):
- Manhattan plot icon + label: “GWAS: GCST005881 (DKD in diabetes; GRCh37)”
- Two kidney atlas cards: “GSE131882 (primary)” and “GSE195460 (replication)”

Column B (Primary localization):
Header (short): “Single-cell integration + scPagwas”
Show a Seurat integration icon feeding into a UMAP thumbnail and then a small heatmap/regression icon labeled “scPagwas”.
UMAP callouts (small, one line each):
- Podocyte (C12): FDR=5.48e−4
- Proximal tubule S3-like (C3): FDR=5.48e−4
- Collecting duct principal-like (C0): FDR=0.00125
- Injury-associated epithelial (C13): FDR=0.00161
QC strip (single line):
N=20,243; exclude 4,043/92,186 (4.39%)

Column C (Secondary follow-up) as four compact stacked tiles:
Tile 1: “Cross-GWAS robustness” with a small dot-matrix icon labeled: DKD meta (GENIE); CKD (CKDGen); eGFR (CKDGen). Small annotation: “rank concordance: 0.42, 0.41, −0.04”.
Tile 2: “Screening MR (OpenGWAS)” with a small forest icon and gene labels: PAX8-AS1; JAG1; H3-3A; ID3.
Tile 3: “Colocalization (coloc.abf)” with an overlap icon and “max PP.H4≈0.034”.
Tile 4: “CellChat (exploratory)” with a network icon and pathway tags: COLLAGEN/LAMININ; VEGF.

Color rules:
- Avoid a cold monochrome look.
- Use soft neutral panel backgrounds (very light warm gray) with subtle colored header bars.
- Accents: warm orange for DKD-localized contexts; muted teal/blue for eGFR-related elements; neutral charcoal for outlines/text.
- Keep overall palette cohesive and publication-friendly; no neon colors.

Do not introduce any additional numbers or claims beyond those specified above.
- The phrase “Multi-Omics”
- The phrase “Causal Verification”
- Any claim that MR/coloc proves causality
