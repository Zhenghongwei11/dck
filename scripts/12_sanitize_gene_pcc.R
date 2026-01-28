args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
  if (dir.exists(repo_root)) setwd(repo_root)
}

source("scripts/utils/repro.R")
ezhu_set_repo_root()

options(stringsAsFactors = FALSE, width = 120)

p_floor <- suppressWarnings(as.numeric(Sys.getenv("EZHU_PVALUE_FLOOR", unset = "1e-300")))
if (!is.finite(p_floor) || p_floor <= 0) p_floor <- 1e-300

root <- Sys.getenv("EZHU_SCPAGWAS_ROOT", unset = "results/scpagwas")
root <- trimws(root)
if (!nzchar(root)) root <- "results/scpagwas"

if (!dir.exists(root)) stop("Missing scPagwas results directory: ", root, call. = FALSE)

files <- list.files(root, pattern = "gene_PCC\\.csv$", recursive = TRUE, full.names = TRUE)
files <- sort(files)
if (length(files) == 0) {
  message("No gene_PCC.csv files found under: ", root)
  quit(status = 0)
}

sanitize_one <- function(path) {
  df <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0) return(invisible(FALSE))

  needed <- c("PCC", "pvalue", "adj_logp", "adj_pvalue", "weight_pcc")
  if (!all(needed %in% names(df))) return(invisible(FALSE))

  df$PCC <- suppressWarnings(as.numeric(df$PCC))
  df$pvalue <- suppressWarnings(as.numeric(df$pvalue))
  df$adj_pvalue <- suppressWarnings(as.numeric(df$adj_pvalue))

  # Underflow guard: p=0 is not a real p-value; treat as numeric underflow.
  df$pvalue[is.finite(df$pvalue) & df$pvalue <= 0] <- p_floor
  df$adj_pvalue[is.finite(df$adj_pvalue) & df$adj_pvalue <= 0] <- p_floor

  # Clamp non-finite adjusted p-values to [p_floor, 1].
  df$adj_pvalue[!is.finite(df$adj_pvalue)] <- NA_real_
  df$adj_pvalue[df$adj_pvalue < p_floor] <- p_floor
  df$adj_pvalue[df$adj_pvalue > 1] <- 1

  # Recompute to avoid Inf/NaN in outputs (e.g., -log10(0)).
  df$adj_logp <- -log10(df$adj_pvalue)
  df$adj_logp[!is.finite(df$adj_logp)] <- NA_real_
  df$adj_logp[df$adj_logp < 0] <- 0

  df$weight_pcc <- df$PCC * df$adj_logp
  df$weight_pcc[!is.finite(df$weight_pcc)] <- NA_real_

  utils::write.csv(df, file = path, row.names = FALSE, quote = TRUE)
  invisible(TRUE)
}

n_ok <- 0L
n_total <- length(files)
for (path in files) {
  ok <- tryCatch(sanitize_one(path), error = function(e) FALSE)
  if (isTRUE(ok)) n_ok <- n_ok + 1L
}

message("Sanitized ", n_ok, " / ", n_total, " gene_PCC.csv files under ", root, " (p_floor=", format(p_floor, scientific = TRUE), ").")

