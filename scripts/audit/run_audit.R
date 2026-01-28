args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = FALSE)
  if (dir.exists(repo_root)) setwd(repo_root)
}

source("scripts/utils/manifest.R")
source("scripts/audit/audit_utils.R")
source("scripts/audit/check_citations.R")
source("scripts/audit/check_renv.R")
source("scripts/audit/check_manifest.R")
source("scripts/audit/check_files.R")
source("scripts/audit/check_scripts.R")
source("scripts/audit/check_figure_index.R")
source("scripts/audit/check_results.R")
source("scripts/audit/check_figures.R")

args <- commandArgs(trailingOnly = TRUE)
level <- audit_parse_level(args)
audit_assert(level %in% c("full", "ci"), "Unknown audit level: ", level)

audit_info("Audit level: ", level)

audit_check_renv()

manifest <- ezhu_read_manifest("data/manifest.tsv")
audit_check_citations()
audit_check_scripts()
audit_check_manifest(manifest)
audit_check_files(manifest, level = level)
audit_check_results(level = level)
audit_check_figure_index(level = level)
audit_check_figures(level = level)

audit_info("All audits passed.")
