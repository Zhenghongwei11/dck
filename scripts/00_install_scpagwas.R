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

options(error = function() {
  traceback(40)
  quit(status = 1, save = "no")
})

options(expressions = max(500000, getOption("expressions", 5000)))

repos <- c(CRAN = "https://cloud.r-project.org")
disable_renv <- identical(Sys.getenv("EZHU_DISABLE_RENV", unset = ""), "1")

cran_pkgs <- c(
  "Matrix",
  "Rcpp",
  "R.methodsS3",
  "R.oo",
  "R.utils",
  "data.table",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "ggsci",
  "ggtext",
  "ggthemes",
  "gridExtra",
  "glmnet",
  "irlba",
  "reshape",
  "reshape2",
  "SOAR",
  "bigreadr",
  "biganalytics",
  "bigmemory",
  "bigstatsr",
  "RMTstat"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = repos)
}

bioc_pkgs <- c("IRanges", "GenomicRanges")

bioc_version_for_r <- function() {
  rver <- paste0(R.version$major, ".", R.version$minor)
  if (utils::compareVersion(rver, "4.5") >= 0) return("3.22")
  if (utils::compareVersion(rver, "4.4") >= 0) return("3.20")
  if (utils::compareVersion(rver, "4.3") >= 0) return("3.18")
  if (utils::compareVersion(rver, "4.2") >= 0) return("3.16")
  "3.15"
}

bioc_ver <- Sys.getenv("EZHU_BIOC_VERSION", unset = "")
bioc_ver <- trimws(bioc_ver)
if (!nzchar(bioc_ver)) bioc_ver <- bioc_version_for_r()
message("Bioconductor version target: ", bioc_ver, " (R ", R.version.string, ")")

install.packages(setdiff(cran_pkgs, rownames(installed.packages())), repos = repos)

# scPagwas is not compatible with SeuratObject v5 (GetAssayData(slot=...) is defunct).
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop(
    "Missing package: SeuratObject\n",
    "If you are using micromamba, install Seurat 4 + SeuratObject 4 via env/ezhu-r-seurat4.yml.\n",
    call. = FALSE
  )
}
so_ver <- as.character(utils::packageVersion("SeuratObject"))
if (utils::compareVersion(so_ver, "5.0.0") >= 0) {
  stop(
    "Incompatible SeuratObject version detected: ", so_ver, "\n",
    "scPagwas requires SeuratObject < 5.0.0 (Seurat v4-era APIs).\n",
    "Fix: recreate the micromamba env from env/ezhu-r-seurat4.yml (clean env; avoid CRAN installing Seurat v5).\n",
    call. = FALSE
  )
}

if (requireNamespace("Seurat", quietly = TRUE)) {
  s_ver <- as.character(utils::packageVersion("Seurat"))
  if (utils::compareVersion(s_ver, "5.0.0") >= 0) {
    stop(
      "Incompatible Seurat version detected: ", s_ver, "\n",
      "scPagwas requires Seurat < 5.0.0.\n",
      "Fix: recreate the micromamba env from env/ezhu-r-seurat4.yml.\n",
      call. = FALSE
    )
  }
}

# Force Bioconductor version to match the current R (or EZHU_BIOC_VERSION) to avoid
# "Bioconductor version X requires R version Y" failures after conda upgrades.
BiocManager::install(version = bioc_ver, ask = FALSE, update = FALSE)
BiocManager::install(setdiff(bioc_pkgs, rownames(installed.packages())), ask = FALSE, update = FALSE)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = repos)
}

scpagwas_ref <- Sys.getenv("EZHU_SCPAGWAS_REF", unset = "sulab-wmu/scPagwas@8dd7975e262a4400928601389142b06cb78367a0")
message("Installing scPagwas from: ", scpagwas_ref)
remotes::install_github(scpagwas_ref, dependencies = FALSE, upgrade = "never")

library(scPagwas)
cat("OK: scPagwas ", as.character(utils::packageVersion("scPagwas")), "\n", sep = "")

if (!disable_renv && requireNamespace("renv", quietly = TRUE)) {
  renv::snapshot(prompt = FALSE)
} else {
  message("Skipping renv::snapshot (EZHU_DISABLE_RENV=1 or renv not available).")
}
