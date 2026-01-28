source("scripts/utils/repro.R")
ezhu_set_repo_root()
ezhu_activate_renv()

repos <- c(CRAN = "https://cloud.r-project.org")

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = repos)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = repos)

BiocManager::install("fgsea", ask = FALSE, update = FALSE)
install.packages(c("Matrix", "Seurat"), repos = repos)
remotes::install_github("sulab-wmu/scPagwas", upgrade = "never", dependencies = TRUE)
renv::snapshot(prompt = FALSE)
cat("OK: scPagwas installed. Version: ", as.character(utils::packageVersion("scPagwas")), "
", sep = "")
