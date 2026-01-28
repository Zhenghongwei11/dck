source("scripts/utils/repro.R")
ezhu_set_repo_root()
ezhu_activate_renv()

options(stringsAsFactors = FALSE, width = 120)

# MR packages are not required for the default local reproduction flow.
# Install them only when running causal MR via OpenGWAS.

repos <- c(
  # IEU packages are distributed via r-universe; keep CRAN as fallback.
  "MRCIEU" = "https://mrcieu.r-universe.dev",
  "CRAN" = "https://cloud.r-project.org"
)

pkgs_cran <- c("ieugwasr", "remotes")
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
}

# TwoSampleMR is typically installed from the MRCIEU r-universe (preferred) or GitHub.
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  tryCatch(
    install.packages("TwoSampleMR", repos = repos),
    error = function(e) {
      message("Failed to install TwoSampleMR from r-universe; falling back to GitHub.")
      remotes::install_github("MRCIEU/TwoSampleMR", upgrade = "never", dependencies = TRUE)
    }
  )
}

cat("OK: MR deps installed. Versions:\n", sep = "")
cat("  ieugwasr: ", as.character(utils::packageVersion("ieugwasr")), "\n", sep = "")
cat("  TwoSampleMR: ", as.character(utils::packageVersion("TwoSampleMR")), "\n", sep = "")

