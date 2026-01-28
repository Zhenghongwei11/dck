#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(purrr))

default_report <- file.path("docs", "audit_runs", "current", "data_audit_log.md")
report_file <- Sys.getenv("EZHU_DATA_AUDIT_REPORT", unset = default_report)
report_file <- trimws(report_file)
if (!nzchar(report_file)) report_file <- default_report
REPORT_FILE <- normalizePath(report_file, winslash = "/", mustWork = FALSE)

dir.create(dirname(REPORT_FILE), recursive = TRUE, showWarnings = FALSE)

cat("# 全量数据层审计报告 (Carpet Bombing)\n\n", file = REPORT_FILE)
cat("生成时间:", as.character(Sys.time()), "\n\n", file = REPORT_FILE, append = TRUE)

# 递归寻找所有数据文件
data_files <- list.files("results", pattern = "\\.(csv|tsv)$", recursive = TRUE, full.names = TRUE)

cat("找到", length(data_files), "个数据文件进行审计。\n\n", file = REPORT_FILE, append = TRUE)

audit_file <- function(f) {
  cat("正在审计:", f, "\n")
  
  # 根据后缀读取
  df <- tryCatch({
    if (grepl("\\.csv$", f)) read_csv(f, show_col_types = FALSE)
    else read_tsv(f, show_col_types = FALSE)
  }, error = function(e) return(NULL))
  
  if (is.null(df) || nrow(df) == 0) return(list(file = f, status = "EMPTY/ERROR"))
  
  # 1. 重复行检查
  dupes <- sum(duplicated(df))
  
  # 2. P值检查 (如果列名包含 pvalue)
  p_cols <- grep("pvalue|p.value|p_val", names(df), ignore.case = TRUE, value = TRUE)
  p_errors <- 0
  if (length(p_cols) > 0) {
    for (pc in p_cols) {
      if (is.numeric(df[[pc]])) {
        p_errors <- p_errors + sum(df[[pc]] < 0 | df[[pc]] > 1, na.rm = TRUE)
      }
    }
  }
  
  # 3. 异常值检查 (Inf / NaN)
  numeric_cols <- sapply(df, is.numeric)
  inf_count <- 0
  if (any(numeric_cols)) {
    inf_count <- sum(sapply(df[, numeric_cols], function(x) sum(is.infinite(x) | is.nan(x))), na.rm = TRUE)
  }
  
  return(list(
    file = f,
    rows = nrow(df),
    dupes = dupes,
    p_errors = p_errors,
    inf_count = inf_count
  ))
}

results <- map(data_files, audit_file)

# 写入报告表格
cat("| 文件路径 | 行数 | 重复行 | P值错误 | 异常值(Inf/NaN) | 状态 |\n", file = REPORT_FILE, append = TRUE)
cat("| :--- | :--- | :--- | :--- | :--- | :--- |\n", file = REPORT_FILE, append = TRUE)

for (res in results) {
  status <- if (res$dupes == 0 && res$p_errors == 0 && res$inf_count == 0) "✅ PASS" else "❌ FAIL"
  if (identical(res$status, "EMPTY/ERROR")) {
     cat(sprintf("| %s | - | - | - | - | ⚠️ READ ERROR |\n", res$file), file = REPORT_FILE, append = TRUE)
  } else {
     cat(sprintf("| %s | %d | %d | %d | %d | %s |\n", res$file, res$rows, res$dupes, res$p_errors, res$inf_count, status), file = REPORT_FILE, append = TRUE)
  }
}

cat("\n审计结束。\n")
