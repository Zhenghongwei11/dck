#!/usr/bin/env Rscript
# scripts/11_scpagwas_gwas_controls_summary.R
#
# Purpose:
# - Summarize scPagwas cluster-level results across one primary GWAS and multiple control GWAS.
#
# Outputs (source-of-record tables for the manuscript):
# - results/scpagwas/gwas_controls/controls_summary.tsv
# - results/scpagwas/gwas_controls/controls_top_contexts.tsv
# Note: the publication figure is generated separately (see `scripts/figures/figure6_scpagwas_multi_gwas_robustness.R`).

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
  if (dir.exists(repo_root)) setwd(repo_root)
}

source("scripts/utils/repro.R")
ezhu_set_repo_root()
ezhu_activate_renv()

options(stringsAsFactors = FALSE)

read_merged_celltype <- function(path, gwas_source) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)

  df <- utils::read.csv(path, check.names = FALSE)
  names(df) <- tolower(names(df))

  # scPagwas merged celltype output typically has: celltype, pvalue, fdr
  cluster <- NULL
  if ("celltype" %in% names(df)) {
    cluster <- as.character(df$celltype)
  } else if ("label" %in% names(df)) {
    cluster <- as.character(df$label)
  } else if ("x" %in% names(df)) {
    cluster <- as.character(df$x)
  } else {
    cluster <- rownames(df)
  }
  cluster <- trimws(as.character(cluster))

  if (!"pvalue" %in% names(df)) stop("Missing pvalue column in: ", path, call. = FALSE)
  pvalue <- suppressWarnings(as.numeric(df$pvalue))
  if (any(!is.finite(pvalue))) stop("Non-numeric pvalue values in: ", path, call. = FALSE)

  fdr <- NULL
  if ("fdr" %in% names(df)) {
    fdr <- suppressWarnings(as.numeric(df$fdr))
  } else {
    fdr <- stats::p.adjust(pvalue, method = "BH")
  }

  out <- data.frame(
    cluster = cluster,
    pvalue = pvalue,
    fdr = fdr,
    gwas_source = as.character(gwas_source),
    stringsAsFactors = FALSE
  )
  out
}

pick_merged_celltype_file <- function(dir_path) {
  candidate1 <- file.path(dir_path, "merged_celltype_pvalue.csv")
  if (file.exists(candidate1)) return(candidate1)
  candidates <- list.files(dir_path, pattern = "Merged_celltype_pvalue\\.csv$", full.names = TRUE, ignore.case = TRUE)
  if (length(candidates) > 0) return(candidates[1])
  NA_character_
}

controls_root <- "results/scpagwas/gwas_controls"
output_dir <- controls_root
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

main_result_path <- Sys.getenv("EZHU_MAIN_MERGED_CELLTYPE_PVALUE", unset = "results/scpagwas/repeat_run/merged_celltype_pvalue.csv")
if (!file.exists(main_result_path)) main_result_path <- "results/scpagwas/main/merged_celltype_pvalue.csv"
if (!file.exists(main_result_path)) stop("Primary merged_celltype_pvalue.csv not found under results/scpagwas/{repeat_run,main}/", call. = FALSE)

message("Reading primary GWAS result: ", main_result_path)
df_main <- read_merged_celltype(main_result_path, gwas_source = "Main")

control_dirs <- list.dirs(controls_root, recursive = FALSE)
control_dirs <- control_dirs[basename(control_dirs) != ""]
control_dirs <- control_dirs[dir.exists(control_dirs)]

results_list <- list(Main = df_main)
file_map <- list(Main = main_result_path)

for (d in control_dirs) {
  gwas_id <- basename(d)
  res_file <- pick_merged_celltype_file(d)
  if (is.na(res_file) || !file.exists(res_file)) {
    warning("No merged celltype pvalue CSV found in: ", d, call. = FALSE)
    next
  }
  message("Reading control GWAS: ", gwas_id, " from ", res_file)
  results_list[[gwas_id]] <- read_merged_celltype(res_file, gwas_source = gwas_id)
  file_map[[gwas_id]] <- res_file
}

# Sanity check: controls must not be byte-identical to the primary file.
# This catches accidental copying / miswired inputs.
allow_identical <- identical(Sys.getenv("EZHU_ALLOW_IDENTICAL_CONTROLS", unset = ""), "1")
if (!allow_identical) {
  main_md5 <- as.character(tools::md5sum(main_result_path))
  for (nm in setdiff(names(file_map), "Main")) {
    ctrl_path <- file_map[[nm]]
    if (!file.exists(ctrl_path)) next
    ctrl_md5 <- as.character(tools::md5sum(ctrl_path))
    if (identical(main_md5, ctrl_md5)) {
      stop(
        "Control GWAS output is identical to primary merged_celltype_pvalue.csv.\n",
        "  primary: ", main_result_path, " (md5=", main_md5, ")\n",
        "  control:  ", ctrl_path, " (md5=", ctrl_md5, ")\n",
        "This indicates the control GWAS did not take effect.\n",
        "Set EZHU_ALLOW_IDENTICAL_CONTROLS=1 to bypass (not recommended).",
        call. = FALSE
      )
    }
  }
}

all_res <- do.call(rbind, results_list)
if (is.null(all_res) || nrow(all_res) == 0) stop("No results combined.", call. = FALSE)

# Attach human-readable cluster aliases when available.
alias_path <- Sys.getenv("EZHU_CLUSTER_ALIAS_MAP", unset = "results/figures/cluster_alias_map.tsv")
if (file.exists(alias_path)) {
  alias <- utils::read.delim(alias_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  names(alias) <- tolower(names(alias))
  if (all(c("cluster", "alias") %in% names(alias))) {
    alias$cluster <- trimws(as.character(alias$cluster))
    alias$alias <- trimws(as.character(alias$alias))
    all_res$cluster <- trimws(as.character(all_res$cluster))
    all_res <- merge(all_res, alias[, c("cluster", "alias")], by = "cluster", all.x = TRUE, sort = FALSE)
  }
}
if (!"alias" %in% names(all_res)) all_res$alias <- NA_character_
all_res$alias <- trimws(as.character(all_res$alias))
all_res$cluster_label <- ifelse(is.na(all_res$alias) | !nzchar(all_res$alias), all_res$cluster, all_res$alias)
all_res$alias <- NULL

summary_file <- file.path(output_dir, "controls_summary.tsv")
utils::write.table(all_res, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote summary: ", summary_file)

# Identify top contexts in the primary GWAS, then extract them across all GWAS.
main_sig <- df_main[df_main$fdr < 0.05, , drop = FALSE]
if (nrow(main_sig) == 0) {
  message("No significant clusters in primary GWAS (FDR < 0.05); using top 10 by p-value.")
  main_sig <- df_main[order(df_main$pvalue), , drop = FALSE]
  main_sig <- main_sig[seq_len(min(10L, nrow(main_sig))), , drop = FALSE]
}
top_clusters <- unique(trimws(as.character(main_sig$cluster)))
comparison <- all_res[all_res$cluster %in% top_clusters, , drop = FALSE]

top_contexts_file <- file.path(output_dir, "controls_top_contexts.tsv")
utils::write.table(comparison, top_contexts_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote top-context comparison: ", top_contexts_file)

ezhu_write_stage_metadata(
  "11_scpagwas_gwas_controls_summary",
  params = list(
    main_result_path = main_result_path,
    controls_root = controls_root,
    outputs = list(
      controls_summary = summary_file,
      controls_top_contexts = top_contexts_file
    )
  ),
  seed = NA_integer_
)
