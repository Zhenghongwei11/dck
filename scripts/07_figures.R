source("scripts/utils/repro.R")
ezhu_set_repo_root()
ezhu_activate_renv()

available <- c("2", "3", "4", "5", "6", "7")

seed_env <- Sys.getenv("EZHU_SEED", unset = "")
seed <- suppressWarnings(as.integer(seed_env))
if (is.na(seed)) seed <- 20251227L
if (!nzchar(seed_env)) Sys.setenv(EZHU_SEED = as.character(seed))
set.seed(seed)

parse_fig_list_local <- function(x) {
  x <- trimws(as.character(x))
  if (!nzchar(x)) return(character())
  if (tolower(x) %in% c("all", "*")) return(available)
  parts <- unlist(strsplit(x, "[,;\\s]+"))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

fig_env <- trimws(Sys.getenv("EZHU_FIGURES", unset = ""))
figs <- parse_fig_list_local(fig_env)
if (length(figs) == 0) figs <- available

if ("1" %in% figs) {
  message("Note: Figure 1 is generated externally (graphical abstract/prompt) and is not produced by scripts.")
  figs <- setdiff(figs, "1")
}

unknown <- setdiff(figs, available)
if (length(unknown) > 0) {
  warning("Ignoring unknown figure IDs: ", paste(unknown, collapse = ", "))
  figs <- intersect(figs, available)
}

if (length(figs) == 0) stop("No valid figures selected. Use EZHU_FIGURES=2,3,4,5,6,7 or EZHU_FIGURES=all.")

fig_scripts <- list(
  "2" = "scripts/figures/figure2_atlas_qc_and_cluster_composition.R",
  "3" = "scripts/figures/figure3_scpagwas_cluster_pvalues_main_vs_repeat.R",
  "4" = "scripts/figures/figure4_top_markers_text.R",
  "5" = "scripts/figures/figure5_mr_top_hits.R",
  "6" = "scripts/figures/figure6_scpagwas_multi_gwas_robustness.R",
  "7" = "scripts/figures/figure7_cellchat_ccc.R"
)

ezhu_write_stage_metadata(
  "07_figures",
  params = list(
    figures = figs,
    formats = Sys.getenv("EZHU_FIGURE_FORMAT", unset = "pdf,png")
  ),
  seed = seed
)

for (fig in figs) {
  script_path <- fig_scripts[[fig]]
  if (!file.exists(script_path)) stop("Missing figure script: ", script_path)
  message("Running Figure ", fig, " via ", script_path)
  script_env <- new.env(parent = globalenv())
  source(script_path, local = script_env)
}

message("Figures finished: ", paste(figs, collapse = ", "))
