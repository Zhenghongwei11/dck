#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

INPUT_FILE <- "results/scpagwas/main/merged_celltype_pvalue.csv"
REPORT_FILE <- "docs/audit_runs/current/figure3_audit_report.md"
dir.create(dirname(REPORT_FILE), recursive = TRUE, showWarnings = FALSE)

cat("Starting audit of", INPUT_FILE, "\n")

if (!file.exists(INPUT_FILE)) {
  stop("Error: File not found at ", INPUT_FILE)
}

df <- read_csv(INPUT_FILE, show_col_types = FALSE)

# 1. Row Count
total_rows <- nrow(df)

# 2. Duplicate Check
num_duplicates <- sum(duplicated(df))

# 3. P-value Integrity
p_value_errors <- df %>% filter(pvalue < 0 | pvalue > 1)
num_p_errors <- nrow(p_value_errors)

# FDR Integrity
fdr_errors <- df %>% filter(fdr < 0 | fdr > 1)
num_fdr_errors <- nrow(fdr_errors)

# FDR Consistency
fdr_consistency_errors <- df %>% filter(fdr < pvalue - 1e-9)
num_fdr_consistency_errors <- nrow(fdr_consistency_errors)

# 4. Generate Report
sink(REPORT_FILE)
cat("# Figure 3 Data Audit Report\n\n")
cat("**Date**:", as.character(Sys.time()), "\n")
cat("**Source File**:", INPUT_FILE, "\n\n")

cat("## Summary\n")
cat("- **Total Rows**:", total_rows, "\n")
cat("- **Duplicate Rows**:", num_duplicates, ifelse(num_duplicates == 0, "✅ PASS", "❌ FAIL"), "\n")
cat("- **P-value Range Errors**:", num_p_errors, ifelse(num_p_errors == 0, "✅ PASS", "❌ FAIL"), "\n")
cat("- **FDR Range Errors**:", num_fdr_errors, ifelse(num_fdr_errors == 0, "✅ PASS", "❌ FAIL"), "\n")
cat("- **FDR < P-value Inconsistencies**:", num_fdr_consistency_errors, ifelse(num_fdr_consistency_errors == 0, "✅ PASS", "❌ FAIL"), "\n\n")

if (num_duplicates > 0) {
  cat("## Duplicate Details\n")
  print(df[duplicated(df), ])
  cat("\n\n")
}

if (num_p_errors > 0) {
  cat("## P-value Error Details\n")
  print(p_value_errors)
  cat("\n\n")
}

sink()

cat("Audit complete. Report generated at", REPORT_FILE, "\n")

if (num_duplicates > 0 || num_p_errors > 0) {
  quit(status = 1)
}
