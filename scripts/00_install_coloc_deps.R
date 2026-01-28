source("scripts/utils/repro.R")
ezhu_set_repo_root()
ezhu_activate_renv()

options(stringsAsFactors = FALSE, width = 120)

repos <- c(
  # IEU packages are distributed via r-universe; keep CRAN as fallback.
  "MRCIEU" = "https://mrcieu.r-universe.dev",
  "CRAN" = "https://cloud.r-project.org"
)

ensure_pkg <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  install.packages(pkg, repos = repos)
  invisible(requireNamespace(pkg, quietly = TRUE))
}

in_conda <- nzchar(Sys.getenv("CONDA_PREFIX", unset = ""))

# Install lightweight CRAN deps from binaries where possible.
pkgs_cran <- c("coloc", "dplyr", "readr", "remotes", "jsonlite")
for (p in pkgs_cran) {
  if (!ensure_pkg(p)) stop("Failed to install package: ", p, call. = FALSE)
}

# OpenGWAS API backend (required): used to query exposure/outcome associations by rsID.
if (!requireNamespace("ieugwasr", quietly = TRUE)) {
  tryCatch(
    install.packages("ieugwasr", repos = repos, dependencies = FALSE),
    error = function(e) {
      # ieugwasr is distributed on MRCIEU r-universe; keep a GitHub fallback.
      remotes::install_github("MRCIEU/ieugwasr", upgrade = "never", dependencies = FALSE)
    }
  )
}

# Optional VCF backend (disabled by default): requires Bioconductor binaries.
use_vcf <- identical(Sys.getenv("EZHU_COLOC_USE_VCF", unset = "0"), "1")
if (isTRUE(use_vcf)) {
  # Bioconductor deps required by gwasvcf (avoid source compilation on fresh VMs).
  bioc_needed <- c("VariantAnnotation", "GenomicRanges")
  bioc_missing <- bioc_needed[!vapply(bioc_needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(bioc_missing) > 0 && isTRUE(in_conda)) {
    msg <- c(
      "Missing Bioconductor packages: ",
      paste(bioc_missing, collapse = ", "),
      "",
      "In micromamba/conda environments, do NOT install these via BiocManager (it may compile from source and fail on modern GCC).",
      "Install prebuilt binaries instead, e.g.:",
      "",
      "  micromamba install -y -n \"$EZHU_MAMBA_ENV\" -c conda-forge -c bioconda \\",
      "    bioconductor-variantannotation bioconductor-genomicranges \\",
      "    bioconductor-rtracklayer bioconductor-rsamtools bioconductor-rhtslib",
      "",
      "Then rerun:",
      "  ./scripts/run_r_mamba.sh scripts/00_install_coloc_deps.R"
    )
    stop(paste(msg, collapse = "\n"), call. = FALSE)
  }

  if (length(bioc_missing) > 0 && !isTRUE(in_conda)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = repos)
    BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
  }

  # Install gwasvcf without pulling Bioc deps (we install those via conda above).
  if (!requireNamespace("gwasvcf", quietly = TRUE)) {
    tryCatch(
      install.packages("gwasvcf", repos = repos, dependencies = FALSE),
      error = function(e) {
        remotes::install_github("MRCIEU/gwasvcf", upgrade = "never", dependencies = FALSE)
      }
    )
  }
}

cat("OK: coloc deps installed. Versions:\n", sep = "")
cat("  coloc: ", as.character(utils::packageVersion("coloc")), "\n", sep = "")
cat("  jsonlite: ", as.character(utils::packageVersion("jsonlite")), "\n", sep = "")
cat("  ieugwasr: ", as.character(utils::packageVersion("ieugwasr")), "\n", sep = "")
if (isTRUE(use_vcf) && requireNamespace("gwasvcf", quietly = TRUE)) {
  cat("  gwasvcf: ", as.character(utils::packageVersion("gwasvcf")), "\n", sep = "")
}
if (isTRUE(use_vcf) && requireNamespace("VariantAnnotation", quietly = TRUE)) {
  cat("  VariantAnnotation: ", as.character(utils::packageVersion("VariantAnnotation")), "\n", sep = "")
}
